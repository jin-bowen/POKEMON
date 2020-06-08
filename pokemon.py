from lib.process import *
from lib.score import *
from lib.skat import *

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_name", type=str,help="gene name(capitalized)")
	parser.add_argument("--genotype", type=str,help="vcf file in plink format")
	parser.add_argument("--cov_file", type=str,help="cov_file file in plink format")
	parser.add_argument("--cov_list", type=str,help="individual cov_file name, seperated by comma")

	parser.add_argument("--reference", type=str,help="snp to AA mapping reference")
	#parser.add_argument("--out_dir", type=str,help="output path for pickle file")

	args = parser.parse_args()

	gene_name = args.gene_name
	genetype  = args.genotype
	cov_file  = args.cov_file
	cov_list  = args.cov_list.split(',')
	reference = args.reference

	# generate input file
	genotype,freqs,pheno,snps2aa,cov = \
		generate(gene_name,genetype,cov_file,cov_list,reference)
	
	# get distance matrix
	distance_mat = cal_distance_mat(snps2aa, freqs)

	# generate the score matrix based on frequency and distance
	freq, score, snp_score_mat = score_mat(freqs.values, distance_mat)
	
	# calculate kernel based on the score matrix
	K = DenseKernel(snp_score_mat, genotype)

	#determine if K is sparse or dense matrix
	if sp.sparse.issparse(K): K = K.toarray()

	# if there is only one element in the kernel
	# do not execute the calculation
	if K.shape[0] < 3 or len(pheno.unique()) < 2:
	        sys.stdout.write('%s\tNA'%pdbid+ '\n')
	        return None

	obj = SKAT(K, fixed_covariates=cov.values)
	pvals = obj.test(pheno.values, acc=1e-7)

	record = [gene_name, str(pvals)]
	sys.stdout.write('\t'.join(record) + '\n')


if __name__ == "__main__":
	main()

