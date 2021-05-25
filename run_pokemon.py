from lib.process import *
from lib.score import *
from lib.VCT import *
import argparse
import os 

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_name", type=str,help="gene name(capitalized)")
	parser.add_argument("--genotype", type=str,help="vcf file in plink format")
	parser.add_argument("--phenotype", type=str)	
	parser.add_argument("--cov_file", type=str,default=None)
	parser.add_argument("--cov_list", type=str,default=None)
	parser.add_argument("--annotation", type=str,help="snp to AA mapping reference")
	parser.add_argument("--ref_pdb_dir", type=str,help="AA coordinate reference")
	parser.add_argument("--alpha", type=str, help="proportion of frequency kernel involved")
	parser.add_argument("--use_blosum", action='store_true')
	parser.add_argument("--use_bfactor", action='store_true')
	parser.add_argument("--out_file", type=str, help="output file")
	parser.add_argument("--pdb", type=str, default=None)
	parser.add_argument("--figures", action='store_true')
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
	ref_pdb_dir   = args.ref_pdb_dir
	alpha         = args.alpha
	out_file      = args.out_file
	use_blosum_bool  = args.use_blosum
	use_bfactor_bool = args.use_bfactor
	use_aa_bool   = args.use_blosum & args.use_bfactor

	if args.pdb: pdb = args.pdb
	else: pdb = None
	draw_figures = args.figures
	
	# default files
	cwd = os.getcwd()
	map_to_pdb_file = cwd + '/ref/pdb_chain_uniprot.csv'
	pwm_file = cwd + '/ref/blosum62'
	pwm = pd.read_csv(pwm_file,index_col=0,delim_whitespace=True)
	
	# generate input file
	genotype,freqs,phenotype,cov = \
	  parser_vcf(genotype_file,phenotype_file,cov_file,cov_list)
	vep = parser_vep(annotation)
	snps = genotype.columns.tolist()
	n_snp = len(snps)
	idx_tab = pd.DataFrame()
	idx_tab['id'] = list(range(n_snp))
	idx_tab['varcode'] = snps
	snps2aa_noidx = snps_to_aa(snps,vep,map_to_pdb_file)
	snps2aa = pd.merge(snps2aa_noidx, idx_tab, on='varcode')

	outf = open(out_file, "a+")
	# no structure mapped 
	if snps2aa.empty: 
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None
	# if there is only one element in the kernel
	# do not execute the calculation
	if genotype.shape[1] < 3:
		outf.write('%s\tNA\tNA\n'%gene_name)
		return None

	# find the protein with most varaints mapped
	if not pdb:
		uniq_map = snps2aa.groupby(['structure','chain'])['varcode'].count().reset_index()
		iline = uniq_map['varcode'].argmax()
		pdb = uniq_map.loc[iline,'structure']

	fig_to_dir = out_file.split('/')[0:-1]
	fig_to_dir = "/".join(fig_to_dir)
	if draw_figures: 
		from lib.plot import score_on_var
		score_on_var(genotype,snps2aa,phenotype,pdb,fig_to_dir)

	#obj = open('%s_%s.pkl'%(gene_name,pdb),'wb')
	#pickle.dump(genotype, obj)
	#pickle.dump(phenotype, obj)
	#pickle.dump(snps2aa, obj)
	#print(len(genotype.columns.tolist()))
	#print(n_snp)
	#print(snps2aa.loc[snps2aa['structure']==pdb,'varcode'].nunique())
	#return 0

	# get distance matrix
	dist_mat_dict = cal_distance_mat(snps2aa, n_snp)
	distance_mat = dist_mat_dict.get(pdb)

	# variants weight induced by aa change
	aa_weight = cal_aa_weight(snps2aa,pwm,n_snp,use_pwm=use_blosum_bool,use_bfct=use_bfactor_bool) 

	# generate the score matrix based on frequency and distance
	# alpha=1: freq only; alpha=0: struct only
	freq_w, struct_w, combined_w = \
	 weight_mat(freqs.values,distance_mat,aa_weight,use_aa=use_aa_bool,alpha=float(alpha))

	# calculate kernel based on the score matrix
	K = cal_Kernel(combined_w, genotype)
	m = genotype.shape[1]

	if args.cov_file:	
		obj = VCT(K, fixed_covariates=cov.values, num_var=m)
	else: obj = VCT(K, num_var=m)

	for lab, ipheno in phenotype.iterrows():
		temp = ipheno.values.astype(np.float)
		pval = obj.test(temp, acc=1e-15)	
		record = [gene_name, pdb,lab, str(pval)]
		outf.write('\t'.join(record) + '\n')
	
if __name__ == "__main__":
	main()

