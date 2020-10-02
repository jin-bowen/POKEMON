from lib.process import *
from lib.score import *
from lib.VCT import *
import argparse

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_name", type=str,help="gene name(capitalized)")
	parser.add_argument("--genotype", type=str,help="vcf file in plink format")
	parser.add_argument("--cov_file", type=str,help="cov_file file in plink format", default=None)
	parser.add_argument("--cov_list", type=str,help="individual cov_file name, seperated by comma", default=None)
	parser.add_argument("--annotation", type=str,help="snp to AA mapping reference")
	parser.add_argument("--ref_pdb_dir", type=str,help="AA coordinate reference")
	parser.add_argument("--out_file", type=str, help="output file")
	parser.add_argument("--alpha", type=str, help="proportion of frequency kernel involved")
	parser.add_argument("--use_aa", action='store_true')

	args = parser.parse_args()

	gene_name = args.gene_name
	genotype_file  = args.genotype
	if args.cov_file:
		cov_file  = args.cov_file
		cov_list  = args.cov_list.split(',')
	else:
		cov_file = None
		cov_list = None
	annotation = args.annotation
	ref_pdb_dir = args.ref_pdb_dir
	alpha     = float(args.alpha)
	out_file  = args.out_file
	use_aa_bool = args.use_aa

	map_to_pdb_file = 'ref/pdbsws_chain'
	#map_to_pdb_file = 'trem2_map'
	pwm_file = 'ref/blosum62'

	outf = open(out_file, "a+")

	# generate input file
	pwm = pd.read_csv(pwm_file, index_col=0,delim_whitespace=True)
	vep = parser_vep(annotation)
	genotype,freqs,pheno,cov = parser_vcf(genotype_file,cov_file,cov_list)
	snps = genotype.columns.tolist()
	snps2aa_noidx = snps_to_aa(snps,gene_name,vep,map_to_pdb_file)

	snps = freqs.index.values
	n_snp = snps.shape[0]
	idx_tab = pd.DataFrame()
	idx_tab['id'] = list(range(n_snp))
	idx_tab['varcode'] = snps
	snps2aa = pd.merge(snps2aa_noidx, idx_tab, on='varcode')

	# no structure mapped 
	if snps2aa.empty: 
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None

	chain_stat = snps2aa.groupby(['structure','chain'])['varcode'].count().reset_index()
	line_num = chain_stat['varcode'].argmax()
	pdb = chain_stat.loc[line_num,'structure']

#	obj = open('%s_stat.pkl'%out_file,'wb')
#	pickle.dump(genotype, obj)
#	pickle.dump(snps2aa, obj)
#	pickle.dump(pheno, obj)
#	obj.close()
#	return 0

	# if there is only one element in the kernel
	# do not execute the calculation
	if genotype.shape[1] < 3 or pheno.nunique() < 2:
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None
	
	# get distance matrix
	dist_mat_dict = cal_distance_mat(snps2aa, n_snp)
	distance_mat = dist_mat_dict.get(pdb)

	aa_weight = cal_aa_weight(snps2aa, pwm, n_snp) 
	# generate the score matrix based on frequency and distance
	# alpha=1: freq only; alpha=0: struct only
	freq_w, struct_w, combined_w = \
	 weight_mat(freqs.values,distance_mat,aa_weight,use_aa=use_aa_bool,alpha = alpha)

	# calculate kernel based on the score matrix
	K = cal_Kernel(combined_w, genotype)
	m = genotype.shape[1]

	if args.cov_file:	
		obj = VCT(K, fixed_covariates=cov.values, num_var=m)
	else: obj = VCT(K, num_var=m)

	pval = obj.test(pheno.values, acc=1e-15)	
	record = [gene_name, pdb, str(pval)]
	outf.write('\t'.join(record) + '\n')

if __name__ == "__main__":
	main()

