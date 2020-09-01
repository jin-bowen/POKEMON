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

	parser.add_argument("--ref_mapping", type=str,help="snp to AA mapping reference")
	parser.add_argument("--ref_pdb_dir", type=str,help="AA coordinate reference")
	parser.add_argument("--empirical_pval", type=str,help="choose to use empirical pvalue or not")
	parser.add_argument("--out_file", type=str, help="output file")

	args = parser.parse_args()

	gene_name = args.gene_name
	genotype_file  = args.genotype
	if args.cov_file:
		cov_file  = args.cov_file
		cov_list  = args.cov_list.split(',')
	else:
		cov_file = None
		cov_list = None
	ref_mapping = args.ref_mapping
	ref_pdb_dir = args.ref_pdb_dir
	empirical = args.empirical_pval
	out_file  = args.out_file

	outf = open(out_file, "a+")

	# generate input file
	genotype,freqs,pheno,snps2aa,cov = \
		generate(gene_name,genotype_file,cov_file,cov_list,ref_mapping,ref_pdb_dir)

	# no structure mapped 
	if snps2aa.empty: return None
	# if there is only one element in the kernel
	# do not execute the calculation
	if genotype.shape[1] < 3 or pheno.nunique() < 2:
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None
	
	# get distance matrix
	dist_mat_dict = cal_distance_mat(snps2aa, freqs)
	
	for pdb, distance_mat in dist_mat_dict.items():
		# generate the score matrix based on frequency and distance
		# alpha=1: freq only; alpha=0: struct only

		freq_w, struct_w, combined_w = \
			sim_mat(freqs.values, distance_mat, alpha = 0.0, rho=0.0)

		# calculate kernel based on the score matrix
		K = cal_Kernel(combined_w, genotype)
		#determine if K is sparse or dense matrix
		if sp.sparse.issparse(K): K = K.toarray()

		if args.cov_file:	
			obj = VCT(K, fixed_covariates=cov.values)
		else: obj = VCT(K)
		pval = obj.test(pheno.values, acc=1e-4)	
		if empirical:
			N = np.ceil(1/p).astype(int)
			pvals_w_shuffle = np.zeros(N)
			for i in range(N):
				np.random.shuffle(phenotype_val)
				pval_shuffle = obj.test(phenotype_val, acc=1e-4)
				pvals_w_shuffle[i] = pval_shuffle
			pvals_bool = pvals_w_shuffle < pval
			S = np.sum(pvals_w_shuffle < pval)
			em_p = S/N
			record = [gene_name, pdb, str(em_p)]
			outf.write('\t'.join(record) + '\n')
		else:
			record = [gene_name, pdb, str(pval)]
			outf.write('\t'.join(record) + '\n')

if __name__ == "__main__":
	main()

