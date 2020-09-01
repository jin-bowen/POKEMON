import dask.dataframe as dd
import pandas as pd
import numpy as np
import os.path

def snps_to_aa(snps,gene_name,ref_mapping,ref_pdb_dir): 

	cols = ['varcode','structure','chain','structure_position','x','y','z']
	empty_df = pd.DataFrame(columns=cols)	

	ref_mapping_grp = ref_mapping.groupby('transcript')

	gene_ref_mapping = ref_mapping_grp.get_group(gene_name)
	if len(gene_ref_mapping.index) == 0: 
		print("no PDB structure mapped to snps")
		return empty_df

	snp_ref_mapping = gene_ref_mapping.loc[gene_ref_mapping['varcode'].isin(snps)]	
	keys = snp_ref_mapping['structure'].unique().compute().tolist()

	ref_pdb = pd.DataFrame()	
	for key in keys:	
		pdb_path = ref_pdb_dir + '/' + str(key)
		if not os.path.exists(pdb_path): continue
		sub_ref_pdb = pd.read_csv(pdb_path, dtype=str)
		ref_pdb = ref_pdb.append(sub_ref_pdb, ignore_index=True)

	if ref_pdb.empty:
		print("no PDB structure file available")
		return empty_df

	ref_pdb_dd = dd.from_pandas(ref_pdb, npartitions=3)
	mapped_record = dd.merge(ref_pdb_dd, snp_ref_mapping, \
		on=['structure','chain','structure_position'], how='inner')

	out_df = mapped_record[cols].compute()
	return	out_df.drop_duplicates().reset_index(drop=True)

def generate(gene_name,genetype,cov_file,cov_list,ref_mapping,ref_pdb_dir):

	df_raw=pd.read_csv(genetype, sep=' ')
	df_raw.set_index('IID', inplace=True)
	df_raw.fillna(0, inplace=True)

	if cov_file:
		cov_raw = pd.read_csv(cov_file, sep=' ')
		cov_raw.set_index('IID', inplace=True)
	else: 
		cov_raw = None
		cov     = None

	# filter individual that carries at least one variant
	filtered_ind = list(map(lambda line: np.sum(line) > 0, df_raw.iloc[:,5:].values))	
	ref_mapping = dd.read_csv(ref_mapping, sep="\t", dtype=str)
	
	df    = df_raw.iloc[filtered_ind,5:]
	pheno = df_raw.loc[filtered_ind,'PHENOTYPE']
	if cov_file:
		cov = cov_raw.loc[filtered_ind,cov_list]

	# accomodate plink phenotype: 1 for control and 2 for case
	pheno = pheno - 1

	# change plink style 12:56477541:C:T_T  to 12:56477541:C:T 
	df_rename   = list(map(lambda s: s[:-2], df.columns.values))
	df.columns  = df_rename

	snps    = df.columns.tolist()
	snps2aa = snps_to_aa(snps, gene_name, ref_mapping, ref_pdb_dir)

	## filter snps with mapped coordinates
	#snps_mapped = set(snps2aa['varcode'].values)
	df_mapped   = df

	## calculate freq
	freqs_mapped = df.sum(axis=0) 
	freqs_mapped = freqs_mapped / df_raw.shape[0]

#	# check the number of variants per individual
#	mapped_ind = list(map(lambda line: np.sum(line) > 0, df_mapped.iloc[:,5:].values))
#	num_case   = pheno.loc[ mapped_ind].values.sum()
#
#	print("filtered_individual:%s"%str(len(mapped_ind)))
#	print("original_individual:%s"%str(len(filtered_ind)))
#	print("case_individual:%s"%str(num_case))

	return df_mapped, freqs_mapped, pheno, snps2aa, cov

if __name__ == "__main__":
	main()


