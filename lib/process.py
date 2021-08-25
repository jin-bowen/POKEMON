from Bio.PDB import *
import pandas as pd
import numpy as np
import urllib.request
import os.path
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

	return csq_df_filter[['ID','Feature','SWISSPROT','Protein_position','Amino_acids']]

def snps_to_aa(snps,vep,map_to_pdb_file,database='pdb'): 

	map_to_pdb = pd.read_csv(map_to_pdb_file,usecols=range(3),index_col=False,\
			comment='#',header=0,names=['structure','chain','SWISSPROT'])
	vep_mapping_processing = vep.loc[vep['ID'].isin(snps)]
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

	ori_cols = ['ID','structure','chain','Protein_position','x','y','z','Amino_acids']
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
			resol = structure.header["resolution"]
			atom  = structure[0][chain][residue]["CA"]
			coord = atom.get_coord()
			#bfct  = atom.get_bfactor() / (resol*resol)
			vep_mapping.loc[irow,['x','y','z']] = coord
			#vep_mapping.loc[irow,'scaled_bfct'] = bfct
		except: continue
	out_df = vep_mapping[ori_cols]
	out_df.columns = new_cols

	return	out_df.dropna().drop_duplicates().reset_index(drop=True)

def map_alphafold_structure(vep_mapping):

	ori_cols = ['ID','structure','chain','Protein_position','x','y','z','Amino_acids']
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

def parser_vcf(genotype_file,phenotype_file,cov_file,cov_list,freq=None):

	# process genotype file
	genotype_raw=pd.read_csv(genotype_file,sep='\s+|\t|,',engine='python',index_col=0)
	genotype_raw.fillna(0,inplace=True)
	genotype = genotype_raw.iloc[:,5:]
	# process covariates
	if cov_file:
		cov_raw = pd.read_csv(cov_file,sep='\s+|\t|,',engine='python',index_col=0)
		cov_raw.dropna(inplace=True)
		cov = cov_raw.loc[:,cov_list]
	else: cov = None

	# process phenotype files
	phenotype = pd.read_csv(phenotype_file,sep='\s+|\t|,',engine='python',index_col=0)
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
	freqs_df = genotype.sum(axis=0) 
	freqs_df = freqs_df /(2 * genotype_raw.shape[0])
	if freq:
		freqs_subset = freqs_df[(freqs_df < freq) & (freqs_df > 0)]
		snps = freqs_subset.index.tolist()
	else:
		freqs_subset = freqs_df[(freqs_df > 0)]
		snps = freqs_subset.index.tolist()

	return genotype.loc[individual,snps],freqs_subset,phenotype[individual], cov

if __name__ == "__main__":
	main()


