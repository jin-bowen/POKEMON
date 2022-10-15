from lib.process import *
from lib.score import *
from lib.VCT import *
import argparse
import os,io 

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_name", type=str,help="gene name(capitalized)")
	parser.add_argument("--genotype", type=str,help="vcf file in plink format")
	parser.add_argument("--phenotype", type=str)	
	parser.add_argument("--cov_file", type=str,default=None)
	parser.add_argument("--cov_list", type=str,default=None)
	parser.add_argument("--annotation", type=str,help="snp to AA mapping reference")
	parser.add_argument("--alpha", type=str, help="proportion of frequency kernel involved")
	parser.add_argument("--use_blosum", action='store_true')
	parser.add_argument("--out_file", type=str, help="output file")
	parser.add_argument("--pdb", type=str, default=None)
	parser.add_argument("--figures", action='store_true')
	parser.add_argument("--out_fig_dir", type=str, default=None)

	parser.add_argument("--maf",type=float,default=0.05)
	parser.add_argument("--database",type=str,default='pdb',\
			help="can be pdb or alphafold, default as pdb")

	args = parser.parse_args()

	gene_name      = args.gene_name
	genotype_file  = args.genotype
	phenotype_file = args.phenotype
	if args.cov_file:
		cov_file  = args.cov_file
		cov_list  = args.cov_list.split(',')
	else:
		cov_file  = None
		cov_list  = None

	annotation    = args.annotation
	alpha         = float(args.alpha)
	out_file      = args.out_file
	use_blosum_bool  = args.use_blosum
	freq_filter   = args.maf
	database      = args.database

	if args.pdb: pdb = args.pdb
	else: pdb = None

	draw_figures = args.figures
	figures_dir  = args.out_fig_dir
	
	# default files
	cwd = os.getcwd()
	map_to_pdb_file = cwd + '/ref/pdb_chain_uniprot.csv.gz'
	map_to_pdb = pd.read_csv(map_to_pdb_file, encoding='utf8',\
				usecols=range(3),index_col=False,compression='gzip',\
				comment='#',header=0,names=['structure','chain','SWISSPROT'])
	
	pwm_file = cwd + '/ref/blosum62.dat'
	pwm = pd.read_csv(pwm_file,encoding='utf8',index_col=0,delim_whitespace=True)
	
	# generate input file
	genotype,freqs,phenotype,cov = \
		parser_vcf(genotype_file,phenotype_file,cov_file,cov_list,freq_filter)
	vep = parser_vep(annotation)
	snps2aa_noidx = snps_to_aa(vep,genotype,map_to_pdb,database=database)

	outf = open(out_file, "a+")
	# no structure mapped 
	if snps2aa_noidx.empty: 
		outf.write('%s\tNA\tNA\n'%gene_name)
		print("no structure mapped")
		return None

	snps2aa,pdb = filter_snps2aa(snps2aa_noidx, pdb=pdb)

	# restrict to mapped variant only
	snps_mapped = snps2aa['varcode'].unique()
	genotype = genotype[snps_mapped]
	freqs    = freqs[snps_mapped]
	# if there is only one element in the kernel, do not execute the calculation
	snps_sum = genotype.sum(axis=0)
	snps_sum = snps_sum.loc[snps_sum>0]
	snps2aa_subset = snps2aa.merge(snps_sum.to_frame(),left_on='varcode',right_index=True)	
	if snps2aa_subset['varcode'].nunique() < 5:
		print("# variance mapped  < 5")
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None
	
	# get distance matrix
	distance_mat = cal_distance_mat(snps2aa)

	# variants weight induced by aa change
	aa_weight = cal_aa_weight(snps2aa,pwm,use_pwm=use_blosum_bool)

	# generate the score matrix based on frequency and distance
	# alpha=1: freq only; alpha=0: struct only
	freq_w, struct_w, combined_w = \
		weight_mat(freqs.values,distance_mat,aa_weight,use_aa=use_blosum_bool,alpha=float(alpha))

	# calculate kernel based on the score matrix
	K = cal_Kernel(combined_w, genotype)
	m = genotype.shape[1]

	if args.cov_file:	
		obj = VCT(K, fixed_covariates=cov.values, num_var=m)
	else: obj = VCT(K, num_var=m)

	for lab, ipheno in phenotype.iterrows():
		temp = ipheno.values.astype(np.float64)
		ipheno_T = ipheno.to_frame(lab).T
		number_pheno_val = len(unique(temp))

		if (number_pheno_val < 3):
			percent_df = percent(genotype,ipheno_T,snps2aa)

		if np.all(percent_df['es']<0.5):
			print("all case or control varaints")
			continue

		pval = obj.test(temp, acc=1e-8)	
		record = [gene_name, pdb,lab, str(pval)]
		outf.write('\t'.join(record) + '\n')
		if pval > 0.05: continue
		if draw_figures and (number_pheno_val < 3): 
			from lib.cluster import cluster
			out_fig_prefix = figures_dir + '/' + gene_name + '_' + lab  +'_' + pdb
			cls = cluster(genotype,snps2aa,ipheno_T,distance_mat,pdb)
			cls.plot(out_fig_prefix)
			cls.plot_cluster(out_fig_prefix)

if __name__ == "__main__":
	main()


