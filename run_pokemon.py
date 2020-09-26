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
	map_to_pdb_file = 'ref/pdbsws_chain'

	outf = open(out_file, "a+")

	# generate input file
	vep = parser_vep(annotation)
	genotype,freqs,pheno,cov = parser_vcf(genotype_file,cov_file,cov_list)
	snps = genotype.columns.tolist()
	snps2aa = snps_to_aa(snps,gene_name,vep,map_to_pdb_file)

	chain_stat = snps2aa.groupby(['structure','chain'])['varcode'].count().reset_index()
	line_num = chain_stat['varcode'].argmax()
	pdb = chain_stat.loc[line_num,'structure']

	# no structure mapped 
	if snps2aa.empty: return None
	# if there is only one element in the kernel
	# do not execute the calculation
	if genotype.shape[1] < 3 or pheno.nunique() < 2:
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None
	
	# get distance matrix
	dist_mat_dict = cal_distance_mat(snps2aa, freqs)
	distance_mat = dist_mat_dict.get(pdb)

	# generate the score matrix based on frequency and distance
	# alpha=1: freq only; alpha=0: struct only
	freq_w, struct_w, combined_w = \
		sim_mat(freqs.values, distance_mat, alpha = alpha, rho=0.0)

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

