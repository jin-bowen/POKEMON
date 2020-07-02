import statsmodels.discrete.discrete_model as sm
import numpy as np 
import pandas as pd
import pickle
import math
import sys
sys.path.insert(1, '/home/bxj139/POKEMON/lib')
from process import *
from score import *
from VCT import *

def main():
	pdb_combo_file = sys.argv[1]
	
	pdb_combo = open(pdb_combo_file,'rb')
	genotype_df = pickle.load(pdb_combo)
	phenotype_df= pickle.load(pdb_combo) 
	snps2aa_df  = pickle.load(pdb_combo) 
	var_es_df   = pickle.load(pdb_combo)

	case_df = phenotype_df[phenotype_df['phenotype'] < 0.5].head(2000)
	ctrl_df = phenotype_df[phenotype_df['phenotype'] > 0.5].head(2000)
	case_ind = case_df.index.tolist()
	ctrl_ind = ctrl_df.index.tolist()
	tot_ind  = case_ind + ctrl_ind

	## filter snps with mapped coordinates
	snps_mapped = snps2aa_df['varcode'].values
	phenotype = phenotype_df.loc[tot_ind]
	genotype  = genotype_df.loc[tot_ind, snps_mapped]

	print(var_es_df.loc[snps_mapped].sort_values(by=['es']))

	phenotype.to_csv('phenotype')
	genotype.to_csv('genotype')

	phenotype_val = phenotype.values.flatten()	
	genotype_val = genotype.values
	freqs = genotype.sum(axis=0)/len(tot_ind)

	clf = sm.Logit(phenotype_val, genotype_val)
	clf_es = clf.fit()

	#print(np.all(genotype.index == phenotype.index))
	var_es_df.loc[ snps_mapped, 'lr_coef'] = clf_es.params
	filename=pdb_combo_file.split('.')[0]
	var_es_df.to_csv( filename + '.csv')
	
	# get distance matrix
	dist_mat_dict = cal_distance_mat(snps2aa_df, freqs)
	
	for pdb, distance_mat in dist_mat_dict.items():
		# generate the score matrix based on frequency and distance
		freq_w, struct_w, combined_w = \
		sim_mat(freqs.values, distance_mat, alpha = 1.0, rho=0)
		
		# calculate kernel based on the score matrix
		K = cal_Kernel(combined_w, genotype)
		#determine if K is sparse or dense matrix
		if sp.sparse.issparse(K): K = K.toarray()

		obj = VCT(K)
		pvals = obj.test(phenotype.values, acc=1e-7)
		
		record = [pdb, str(pvals)]
		sys.stdout.write('\t'.join(record) + '\n')

# main body
if __name__ == "__main__":
	main()


