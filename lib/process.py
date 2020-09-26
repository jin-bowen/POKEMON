import dask.dataframe as dd
from Bio.PDB import *
import pandas as pd
import numpy as np
import os.path
import re

def parser_vep(vep_input):
	
	with open(vep_input)as in_file:
		for line in in_file:
			if re.match("^##INFO",line): 
				csq_header_raw = line.split('\"')[1].split(':')[1]
				csq_header = csq_header_raw.split('|')
			if re.match("^#CHROM",line): 
				vcf_header_raw = line.strip('\n').split('\t')
				vcf_header = vcf_header_raw[0:9]

	vep_df = pd.read_csv(vep_input,comment='#',sep='\t',header=None,names=vcf_header, \
				usecols=range(9))
	
	temp_df = vep_df['INFO'].str.split(',', expand=True).add_prefix('csq')
	vep_df = vep_df.join(temp_df)
	temp_header = temp_df.columns

	csq_df = vep_df.melt(id_vars=vcf_header, value_vars=temp_header,value_name='csq')
	temp_df = csq_df['csq'].str.split('|',expand=True)
	temp_df.columns = csq_header

	csq_df[csq_header] = temp_df[csq_header]
	csq_df_filter = csq_df[ (csq_df['Consequence']=='missense_variant') & \
		(csq_df['CANONICAL']=='YES')]

	pdbentry = csq_df_filter['SWISSPROT'].str.split('_',1, expand=True)
	csq_df_filter['SWISSPROT'] = pdbentry[0]
	csq_df_filter['ID']=csq_df_filter[['#CHROM','POS','REF','ALT']].astype(str).agg(':'.join,axis=1)

	return csq_df_filter[['ID','Feature','SWISSPROT','Protein_position','Amino_acids']]

def snps_to_aa(snps,gene_name,vep,map_to_pdb_file): 

	map_to_pdb = pd.read_csv(map_to_pdb_file, header=None,sep=' ',index_col=False,\
			names=['structure','chain','SWISSPROT'])
	vep_mapping_processing = vep.loc[vep['ID'].isin(snps)]
	vep_mapping = pd.merge(vep_mapping_processing, map_to_pdb, on='SWISSPROT')	

	if len(vep_mapping['structure']) == 0: 
		print("no PDB structure mapped to snps")
		return empty_df

	vep_mapping['x'] = np.nan
	vep_mapping['y'] = np.nan
	vep_mapping['z'] = np.nan

	parser = PDBParser()
	structure_list = {}
	for pdb in set(vep_mapping['structure'].values):
		pdbl=PDBList()
		pdbl.retrieve_pdb_file(pdb, pdir='ref', file_format='pdb',obsolete=False)
		structure = parser.get_structure(pdb, "ref/pdb%s.ent"%pdb)
		structure_list[pdb] = structure

	for irow, row in vep_mapping.iterrows():
	
		entry = row['structure']
		chain = row['chain']
		residue = int(row['Protein_position'])

		try:
			structure = structure_list.get(entry)
			atom = structure[0][chain][residue]["CA"]
			coord = atom.get_coord()
			vep_mapping.loc[irow,['x','y','z']] = coord
		except:
			continue
	ori_cols = ['ID','structure','chain','Protein_position','x','y','z']
	new_cols = ['varcode','structure','chain','structure_position','x','y','z']
	out_df = vep_mapping[ori_cols]
	out_df.columns = new_cols

	return	out_df.dropna().drop_duplicates().reset_index(drop=True)

def parser_vcf(genotype,cov_file,cov_list):

	df_raw=pd.read_csv(genotype, sep=' ')
	df_raw.set_index('IID', inplace=True)
	df_raw.fillna(0, inplace=True)

	if cov_file:
		cov_raw = pd.read_csv(cov_file, sep=' ')
		cov_raw.set_index('IID', inplace=True)
	else: 
		cov_raw = None
		cov     = None

	df = df_raw.iloc[:,5:]
	pheno = df_raw.loc[:,'PHENOTYPE'] 

	# for binary phenotype
	if pheno.nunique() == 2:
		pheno = pheno - 1

	if cov_file:
		cov = cov_raw.loc[:,cov_list]
	
	# change plink style 12:56477541:C:T_T  to 12:56477541:C:T 
	df_rename   = list(map(lambda x: x.split('_', 1)[0], df.columns))
	df.columns  = df_rename

	## calculate freq
	freqs = df.sum(axis=0) 
	freqs = freqs /(2 * df_raw.shape[0])

	return df, freqs, pheno, cov

if __name__ == "__main__":
	main()


