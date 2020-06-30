import dask.dataframe as dd
import pandas as pd
import numpy as np
import pickle
import argparse

def snps_to_aa(snps,gene_name,ref_mapping,ref_pdb): 

	ref_mapping_grp = ref_mapping.groupby('transcript')
	ref_pdb_grp = ref_pdb.groupby('structure')

	print(ref_pdb_grp.get_group('2fp0').compute())

	try: gene_ref_mapping = ref_mapping_grp.get_group(gene_name)
	except: 
		print("no coordinate information for input gene")
		return 0
	snp_ref_mapping = gene_ref_mapping.loc[gene_ref_mapping['varcode'].isin(snps)]

	keys = snp_ref_mapping['structure'].unique().compute().tolist()

	fst_key = keys.pop(0)
	sub_ref_pdb = ref_pdb_grp.get_group(fst_key)
	for key in keys:
		temp = ref_pdb_grp.get_group(key)
		sub_ref_pdb = sub_ref_pdb.append(temp)

	mapped_record = dd.merge(sub_ref_pdb, snp_ref_mapping, \
		on=['structure','chain','structure_position'], how='inner')

	cols = ['varcode','structure','chain','structure_position','x','y','z']
	out_df = mapped_record[cols].compute()
	return	out_df.drop_duplicates().reset_index(drop=True)

def generate(gene_name,genetype,cov_file,cov_list,ref_mapping,ref_pdb):

	df_raw=pd.read_csv(genetype, sep=' ')
	df_raw.set_index( 'IID', inplace=True)
	df_raw.fillna(0, inplace=True)

	cov_raw = pd.read_csv(cov_file, sep=' ')
	cov_raw.set_index( 'IID', inplace=True)

	ref_mapping = dd.read_csv(ref_mapping, sep="\t", dtype=str)
	ref_pdb = dd.read_csv(ref_pdb, dtype=str, usecols=list(range(12)))

	# filter individual that carries at least one variant
	filtered_ind = list(map(lambda line: np.sum(line) > 0, df_raw.iloc[:,5:].values))	
	df    = df_raw.iloc[ filtered_ind,5:]
	pheno = df_raw.loc[ filtered_ind,'PHENOTYPE']
	# accomodate plink phenotype: 1 for control and 2 for case
	pheno = pheno - 1

	# change plink style 12:56477541:C:T_T  to 12:56477541:C:T 
	df_rename = map(lambda s: s[:-2], df.columns.values)
	df.columns  = df_rename

	cov  =	cov_raw.loc[ filtered_ind , cov_list ]
	
	## calculate freq
	freqs = df.sum(axis=0) 
	freqs = freqs / df_raw.shape[0]

	## filter snps with freqs > 0.0
	freqs_filtered = freqs[freqs > 0.0]
	df_filtered    = df.loc[:,freqs_filtered.index.values]	

	snps  = df_filtered.columns.tolist()
	snps2aa = snps_to_aa(snps, gene_name, ref_mapping, ref_pdb)
	
	## filter snps with mapped coordinates
	snps_mapped = snps2aa['varcode'].values
	df_clean = df_filtered.loc[:,snps_mapped]
	freqs_clean = freqs_filtered.loc[snps_mapped]

	return df_clean, freqs_clean, pheno, snps2aa, cov

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_name", type=str,help="gene name(capitalized)")
	parser.add_argument("--genotype", type=str,help="vcf file in plink format")
	parser.add_argument("--cov_file", type=str,help="cov_file file in plink format")
	parser.add_argument("--cov_list", type=str,help="individual cov_file name, seperated by comma")

	parser.add_argument("--reference", type=str,help="snp to AA mapping reference")
	parser.add_argument("--out_dir", type=str,help="output path for pickle file")

	args = parser.parse_args()

	gene_name = args.gene_name
	genetype  = args.genotype
	cov_file  = args.cov_file	
	cov_list  = args.cov_list.split(',')	
	reference   = args.reference
	out_dir   = args.out_dir

	df_clean, freqs_clean, pheno, snps2aa, cov = \
		generate(gene_name,genetype,cov_file,cov_list,reference)

	obj = open('%s/%s.pkl'%(out_dir,gene_name),'wb')	
	pickle.dump(df_clean, obj)
	pickle.dump(freqs_clean, obj)
	pickle.dump(pheno,   obj)
	pickle.dump(snps2aa, obj)
	pickle.dump(cov, obj)
	obj.close()

if __name__ == "__main__":
	main()


