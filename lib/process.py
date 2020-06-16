import dask.dataframe as dd
import pandas as pd
import numpy as np
import pickle
import argparse

## function to match variant with amino acid  #
#def snps_aa(snps, pdbid, conn):
#
#	## for individual variant record
#	## match it with the structural information
#	## default mode is 0
#	rec = pd.DataFrame(columns=['gd_id','structid','chain','seq','x','y','z'])
#	for snp in snps: 
#		chrom = snp.split(':')[0]	
#		pos   = snp.split(':')[1]	
#
#		cur=conn.cursor()
#		## find match variant IndVar in genomicData table: gd_id
#		## use gd_id to find corresponding structural information
#		query  = "select gd.name, a.structid, a.chain, a.seqid, a.x, a.y, a.z "
#		query += "from GenomicData gd "
#		query += "inner join wes_missense_%s_pdb a "%chrom
#		query += "on  gd.gd_id =  a.gd_id "	
#		query += "where gd.chr   = 'chr%s' "%chrom 
#		query += "and   gd.start = %s "%pos
#		query += "and   gd.label = 'wes_missense_%s' "%chrom
#		query += "and   a.structid = '%s'"%pdbid 
#	
#		dstream = cur.execute(query)
#		## one variant can code into several different subnit
#		temp = cur.fetchmany(dstream)
#
#		for line in temp: rec.loc[len(rec)] = list(line)
#
#	rec.set_index('gd_id', inplace=True)
#	return rec

def generate(gene_name,genetype,cov_file,cov_list,reference):

	df_raw=pd.read_csv(genetype, sep=' ')
	df_raw.set_index( 'IID', inplace=True)
	df_raw.fillna(0, inplace=True)

	cov_raw = pd.read_csv(cov_file, sep=' ')
	cov_raw.set_index( 'IID', inplace=True)

	ref_raw = dd.read_csv(reference, sep="\t")
	ref = ref_raw.groupby(['gene'])

	# filter individual that carries at least one variant
	filtered_ind = list(map(lambda line: np.sum(line) > 0, df_raw.iloc[:,5:].values))	
	df    = df_raw.iloc[ filtered_ind,5:]
	pheno = df_raw.loc[ filtered_ind,'PHENOTYPE']
	# accomoddate plink phenotype: 1 for control and 2 for case
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
	snps2aa = snps_to_aa(snps, gene_name, ref)
	
	## filter snps with mapped coordinates
	snps_mapped = snps2aa['snp'].values
	df_clean = df_filtered.loc[:,snps_mapped]
	freqs_clean = freqs_filtered.loc[snps_mapped]

	return df_clean, freqs_clean, pheno, snps2aa, cov

def snps_to_aa(snps, gene_name, ref): 

	try: 
		gene_ref = ref.get_group(gene_name).compute()	
	except: 
		print("no coordinate information for input gene")
		return 0

	snps_modified = list(map(lambda x: ('chr' + x).split(':')[0:2], snps))

	gene_ref['snp'] = None 

	for i, isnp in enumerate(snps_modified):
		bool_row = (gene_ref['chr'] == isnp[0]) & \
			(gene_ref['start'].astype(int) == int(isnp[1])) 
		gene_ref.loc[bool_row, 'snp'] = snps[i] 
			
	snps_intersect = gene_ref.dropna()

	cols = ['snp','structure','chain','structure_position','x','y','z']
	return snps_intersect[cols].reset_index(drop=True)
	
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


