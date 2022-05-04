from Bio.PDB import *
import pandas as pd
import numpy as np
import urllib.request
import os.path
import os
import re

def percent(genotype, phenotype,snps2aa):
	geno  = genotype.values
	geno_sum = geno.sum(axis=0)
	pheno = phenotype.values.reshape((1,-1))

	percent = pheno.dot(geno).reshape(-1)
	percent = percent / geno_sum

	percent_df = pd.DataFrame(index=genotype.columns)
	percent_df.loc[:,'es'] = percent
	snp_df = pd.merge(percent_df, snps2aa,left_index=True, right_on='varcode')
	
	return snp_df

def parser_vep(vep_input):
	
	with open(vep_input)as in_file:
		for line in in_file:
			if re.match("^##INFO=<ID=CSQ",line): 
				csq_header_raw = line.split('\"')[1].split(':')[1]
				csq_header = csq_header_raw.split('|')
			if re.match("^#CHROM",line): 
				vcf_header_raw = line.strip('\n').split('\t')
				vcf_header = vcf_header_raw[0:8]
	vep_df = pd.read_csv(vep_input,comment='#',sep='\t',header=None,names=vcf_header, \
				usecols=range(8))
	temp_df = vep_df['INFO'].str.split(',', expand=True).add_prefix('csq')
	vep_df = vep_df.join(temp_df)
	temp_header = temp_df.columns

	csq_df = vep_df.melt(id_vars=vcf_header, value_vars=temp_header,value_name='csq')
	temp_df = csq_df['csq'].str.split('|',expand=True)
	if temp_df.shape[1] < len(csq_header): return pd.DataFrame()
	temp_df.columns = csq_header
	
	csq_df[csq_header] = temp_df[csq_header]
	csq_df_filter = csq_df[ (csq_df['Consequence']=='missense_variant') & \
				(csq_df['CANONICAL']=='YES')]
	pdbentry = csq_df_filter['SWISSPROT'].str.split('_',1, expand=True)
	csq_df_filter.loc[:,'SWISSPROT'] = pdbentry[0]
	csq_df_filter.loc[:,'ID']=csq_df_filter[['#CHROM','POS','REF','ALT']].astype(str).agg(':'.join,axis=1)

	#extend the csq dataframe to all possible amino acid change
	#ref as destinate allele
	csq_df_filter_part1 = csq_df_filter.copy()
	csq_df_filter_part1['des_allele'] = csq_df_filter_part1['REF']
	csq_df_filter_part1['Amino_acids'] = csq_df_filter['Amino_acids'].str[::-1]

	#alt as destinate allele
	csq_df_filter_part2 = csq_df_filter.copy()
	csq_df_filter_part2['des_allele'] = csq_df_filter_part2['ALT']
	csq_df_filter_part2['Amino_acids'] = csq_df_filter['Amino_acids']
	csq_df_filter_ext = pd.concat([csq_df_filter_part1,csq_df_filter_part2],ignore_index=True)

	csq_df_filter_ext['full_id'] = csq_df_filter_ext[['ID','des_allele']].astype(str).agg('_'.join,axis=1)
	return csq_df_filter_ext[['full_id','Feature','SWISSPROT','Protein_position','Amino_acids']]


def filter_snps2aa(snps2aa_noidx, pdb='None'):
	# find the protein with most varaints mapped
	if not pdb:
		uniq_map = snps2aa_noidx.groupby(['structure','chain'])['varcode'].count().reset_index()
		iline = uniq_map['varcode'].argmax()
		pdb = uniq_map.loc[iline,'structure']

	snps2aa = snps2aa_noidx.loc[snps2aa_noidx['structure']==pdb]

	snps_mapped = snps2aa['varcode'].unique().tolist()

	idx_tab = pd.DataFrame()
	idx_tab.loc[:,'id'] = list(range(len(snps_mapped)))
	idx_tab.loc[:,'varcode'] = snps_mapped
	snps2aa = pd.merge(snps2aa, idx_tab, on='varcode')

	return snps2aa,pdb

def snps_to_aa(vep,genotype,map_to_pdb_file,database='pdb'): 

	snps = set(vep['full_id'].values).intersection(genotype.columns.tolist())
	map_to_pdb = pd.read_csv(map_to_pdb_file,usecols=range(3),index_col=False,\
			comment='#',header=0,names=['structure','chain','SWISSPROT'])
	vep_mapping_processing = vep.loc[vep['full_id'].isin(snps)]
	vep_mapping_processing.loc[:,'SWISSPROT'] = vep_mapping_processing['SWISSPROT'].str.split('.').str[0]
	
	if database == "pdb":
		vep_mapping = pd.merge(vep_mapping_processing, map_to_pdb, on='SWISSPROT')
		out_df = map_PDB_structure(vep_mapping)	

	elif database == "alphafold":
		vep_mapping = vep_mapping_processing.copy()
		vep_mapping.loc[:,'structure'] = vep_mapping['SWISSPROT'] 
		out_df = map_alphafold_structure(vep_mapping)

	return out_df
			
def map_PDB_structure(vep_mapping):

	pdb_struct_dir='ref/pdb'
	if not os.path.exists(pdb_struct_dir):os.makedirs(pdb_struct_dir)	

	ori_cols = ['full_id','structure','chain','Protein_position','x','y','z','Amino_acids']
	new_cols = ['varcode','structure','chain','structure_position','x','y','z','aa']

	if len(vep_mapping['structure']) == 0: 
		print("no PDB structure mapped to snps")
		return pd.DataFrame(columns=new_cols)

	vep_mapping.loc[:,'x'] = np.nan
	vep_mapping.loc[:,'y'] = np.nan
	vep_mapping.loc[:,'z'] = np.nan
	
	parser = PDBParser()
	structure_list = {}
	for pdb in set(vep_mapping['structure'].values):
		try:
			pdbl=PDBList()
			pdbl.retrieve_pdb_file(pdb, pdir='ref/pdb', file_format='pdb',obsolete=False)
			structure = parser.get_structure(pdb, "ref/pdb/pdb%s.ent"%pdb)
			structure_list[pdb] = structure
		except: continue

	for irow, row in vep_mapping.iterrows():	
		entry = row['structure']
		chain = row['chain']
		residue = int(row['Protein_position'])
		try:
			structure = structure_list.get(entry)
			# ensure the expression source
			org_sys = structure.header['source']['1']['organism_scientific'].replace(' ','').split(',')
			if len(org_sys) > 1 or org_sys[0] != 'homosapiens': continue 

			atom  = structure[0][chain][residue]["CA"]
			coord = atom.get_coord()
			vep_mapping.loc[irow,['x','y','z']] = coord
		except: continue
	out_df = vep_mapping[ori_cols]
	out_df.columns = new_cols

	return	out_df.dropna().drop_duplicates().reset_index(drop=True)

def map_alphafold_structure(vep_mapping):

	af_struct_dir='ref/alphafold'
	if not os.path.exists(af_struct_dir):os.makedirs(af_struct_dir)	

	ori_cols = ['full_id','structure','chain','Protein_position','x','y','z','Amino_acids']
	new_cols = ['varcode','structure','chain','structure_position','x','y','z','aa']

	if len(vep_mapping['structure']) == 0: 
		print("no alphafold structure mapped to snps")
		return pd.DataFrame(columns=new_cols)

	vep_mapping.loc[:,'x'] = np.nan
	vep_mapping.loc[:,'y'] = np.nan
	vep_mapping.loc[:,'z'] = np.nan

	parser = PDBParser()
	structure_list = {}
	for uniprot in vep_mapping['structure'].unique():
		try:
			chain_id = []
			url='https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v1.pdb'%(uniprot)
			urllib.request.urlretrieve(url, "ref/alphafold/%s.pdb"%uniprot)
			structure = parser.get_structure(uniprot, "ref/alphafold/%s.pdb"%uniprot)
			structure_list[uniprot] = structure
			for chain in structure[0].get_list():
				chain_id += chain.get_id()
			vep_mapping['chain_list'] = ','.join(chain_id)
		except: continue
	vep_mapping.loc[:,'chain'] = vep_mapping.chain_list.str.split(',')
	vep_mapping = vep_mapping.explode('chain')

	for irow, row in vep_mapping.iterrows():	
		entry = row['structure']
		chain = row['chain']
		residue = int(row['Protein_position'])
		try:
			structure = structure_list.get(entry)
			resol = structure.header["resolution"]
			atom  = structure[0][chain][residue]["CA"]
			coord = atom.get_coord()
			vep_mapping.loc[irow,['x','y','z']] = coord
		except: continue
	out_df = vep_mapping[ori_cols]
	out_df.columns = new_cols

	return	out_df.dropna().drop_duplicates().reset_index(drop=True)

def parser_vcf(genotype_file,phenotype_file,cov_file,cov_list,freq):

	# process genotype file
	genotype_raw=pd.read_csv(genotype_file,sep='\s+|\t|,',engine='python',index_col=1)
	genotype = genotype_raw.iloc[:,5:]

	# process covariates
	if cov_file:
		cov_raw = pd.read_csv(cov_file,sep='\s+|\t|,',engine='python',index_col=0)
		cov_raw.dropna(inplace=True)
		cov = cov_raw.loc[:,cov_list]
	else: cov = None

	# process phenotype files
	phenotype = pd.read_csv(phenotype_file,sep='\s+|\t|,',engine='python',usecols=range(2),index_col=0)
	phenotype.dropna(inplace=True)

	if phenotype.shape[1] < 2: phenotype = phenotype.transpose()
	phenotype_ind = phenotype.columns.tolist()

	genotype_ind = genotype.index.tolist()
	individual = set(genotype_ind).intersection(phenotype_ind)
	if cov_file:	
		cov_ind = cov.index.tolist()
		individual = set(individual).intersection(cov_ind)
		cov = cov.loc[individual]

	## calculate freq
	freqs = genotype.sum(axis=0,skipna=True) 
	freqs = freqs /(2 * genotype.count())
	
	freqs_df = freqs.to_frame(name='freq')
	# secure a minor allele count matrix
	
	freqs_df['index'] = freqs_df.index		
	freqs_df[['ID','c_allele']] = freqs_df['index'].str.split('_',expand=True)
	freqs_df[['#CHROM','POS','REF','ALT']] = freqs_df['ID'].str.split(':',expand=True)

	genotype_processing = genotype.T
	flip_allele_bool = freqs > (1-freq)
	genotype_processing.loc[flip_allele_bool] = genotype_processing.loc[flip_allele_bool].replace({2:0, 0:2})
	freqs_df.loc[flip_allele_bool,'index'] = freqs_df.loc[flip_allele_bool,['ID','ALT']].astype(str).agg('_'.join,axis=1)

	genotype_processing.set_index(freqs_df['index'].values, inplace=True)
	ma_genotype = genotype_processing.T

	ma_freqs = ma_genotype.sum(axis=0,skipna=True) 
	ma_freqs = ma_freqs /(2 * ma_genotype.count())
		
	ma_freqs_subset = ma_freqs[(ma_freqs < freq) & (ma_freqs > 0)]
	snps = ma_freqs_subset.index.tolist()

	ma_genotype.fillna(0, inplace=True)
	return ma_genotype.loc[individual,snps],ma_freqs[snps],phenotype[individual], cov

if __name__ == "__main__":
	main()


