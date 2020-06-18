from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from scipy import stats

#import pymol
#from pymol import cmd

import dask.dataframe as dd
import numpy as np 
import pandas as pd
import itertools as it
import pickle
import math
import sys

def generate_phenotype(genotype, snp_corr_es_df):
	"""

	Args:

	Returns: phenotype - list of phenotype for number_ind, 1D array of size number_ind
	"""

	snp_corr_es = snp_corr_es_df.values
	sub_snps    = snp_corr_es_df.index.tolist()

	snp_es = np.sqrt(np.diagonal(snp_corr_es))	
	np.fill_diagonal(snp_corr_es,0)

	ind_genotype = genotype[sub_snps].values

	ind_log_odds =  np.dot(ind_genotype, snp_es)
	
	ind_ind_corr = np.linalg.multi_dot([ind_genotype, snp_corr_es, ind_genotype.T])
	ind_log_odds += np.diagonal(ind_ind_corr)
	ind_log_odds += 1

	ind_prob = np.exp(ind_log_odds) / (1 + np.exp(ind_log_odds))

	phenotype = list(map(lambda prob: np.random.binomial(1, p=prob), list(ind_prob)))

	phenotype_df = pd.DataFrame(phenotype, columns=['phenotype'])

	return phenotype_df

def exponential_es(snp_es, dist_mat_df, t=14):
	"""

	Args:

	Returns:
	"""
	sub_snps = dist_mat_df.index.values
	dist_mat = dist_mat_df.values

	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	if min(dist_list) < 0:
		raise ValueError('invalid distance exist')
	else:
		struct_w_list = list(map(lambda r: np.exp(-r * r/(2 * t * t)), dist_list))

	struct_w = np.array(struct_w_list).reshape(num_row, num_col)
	sub_snp_es = snp_es.loc[sub_snps].values

	snp_snp_corr = np.dot(sub_snp_es, sub_snp_es.T)
	sub_snp_es_combined = snp_snp_corr * struct_w
	
	sub_snp_es_combined_df = pd.DataFrame(sub_snp_es_combined, \
					index=sub_snps,\
					columns=sub_snps)
	return sub_snp_es_combined_df

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
	return snps_intersect[cols].drop_duplicates().reset_index(drop=True)

#######################################
def cal_distance_mat(snps2aa_tot):
	"""
	Args:

	Returns:
	"""

	snps2aa_grp = snps2aa_tot.groupby(['structure'])

	dist_mat_dict = {}
	for key, snps2aa in snps2aa_grp:

		snps2aa['id'] = list(range(len(snps2aa)))
	
		snp_coord = snps2aa[['x','y','z']].values
		distance_vec = squareform(pdist(snp_coord)).flatten()
	
		snp_idx = snps2aa['id'].tolist()
		snp_pair_idx = list(map(list, it.product(snp_idx, repeat=2)))

		snp_pair_distance = pd.DataFrame(data=snp_pair_idx, columns=['i','j'], dtype=int)
		snp_pair_distance['r'] = distance_vec
	
		snp_pair_distance_uniq = snp_pair_distance.groupby(['i','j']).min().reset_index()
	
		row  = snp_pair_distance_uniq['i'].values
		col  = snp_pair_distance_uniq['j'].values
		data = snp_pair_distance_uniq['r'].values
	
		distance_mat = sparse.coo_matrix((data,(row,col))).toarray()
		distance_mat_df = pd.DataFrame(distance_mat, \
					columns=snps2aa['snp'].values,\
					index=snps2aa['snp'].values)	
	
		dist_mat_dict[key] = distance_mat_df
		# check if the distance matrix if symmetrical	
		#print((distance_mat.T == distance_mat).all())
	return dist_mat_dict

def plot(var, pdb_id, chain=None):

	pymol.finish_launching()

	if chain:
		pdb_id = pdb_id + chain

	cmd.fetch(pdb_id)
	cmd.alter(pdb_id, 'b = 0.5')
	cmd.show_as('cartoon',pdb_id)
	cmd.color('white',pdb_id)

	for i, row in var.iterrows():
		resi  = row['seq']
		chain = row['chain']
		pheno = row['es']
		cmd.alter('resi %s and chain %s'%(resi, chain), 'b=%s'%pheno)

	cmd.spectrum("b", "blue_white_red", "%s"%pdb_id, maximum=1.0, minimum=0)
	cmd.zoom()

def main():

	case_variant_file = sys.argv[1]	
	ctrl_variant_file = sys.argv[2]	
	reference = sys.argv[3]

	gene_name = case_variant_file.split('_')[0]

	cols=['chr','pos','if','ref','alt', 'qual', 'filter', 'info', 'format']
	case_variant = pd.read_csv(case_variant_file, sep="\t", \
			header=0, names=cols, usecols=list(range(9))) 
	ctrl_variant = pd.read_csv(ctrl_variant_file, sep="\t", \
			header=0, names=cols, usecols=list(range(9))) 

	case_variant['snp'] = case_variant['chr'].astype(str) + ":" + \
				case_variant['pos'].astype(str) + ":" + \
				case_variant['ref'].astype(str) + ":" + \
				case_variant['alt'].astype(str)

	ctrl_variant['snp'] = ctrl_variant['chr'].astype(str) + ":" + \
				ctrl_variant['pos'].astype(str) + ":" + \
				ctrl_variant['ref'].astype(str) + ":" + \
				ctrl_variant['alt'].astype(str) 


	case_variant_es = pd.DataFrame(columns=['es'], index=case_variant['snp'].values)
	case_variant_es['es'] = 0.03
	ctrl_variant_es = pd.DataFrame(columns=['es'], index=ctrl_variant['snp'].values)
	ctrl_variant_es['es'] = 0.00

	num_ind = 2000
	num_var = case_variant['snp'].values.shape[0] + ctrl_variant['snp'].values.shape[0]
	var =  case_variant['snp'].tolist() + ctrl_variant['snp'].tolist()
	ind_genotype = np.random.randint(low=0, high=2, size=(num_ind, num_var))
	ind_genotype_df = pd.DataFrame(ind_genotype, columns=var)

	ref_raw = dd.read_csv(reference, sep="\t")
	ref = ref_raw.groupby(['transcript'])
	case_snp = case_variant['snp'].values
	snps2aa = snps_to_aa(case_snp, gene_name, ref)

	dist_mat_dict = cal_distance_mat(snps2aa)

	for pdb, dist_mat in dist_mat_dict.items():
		case_variant_corr_es = exponential_es(case_variant_es, dist_mat, t=14)
	
		case_var_es = case_variant_corr_es.sum(axis=1).to_frame(name='es')
			
		phenotype = generate_phenotype(ind_genotype_df,case_variant_corr_es)
		var = pd.concat([case_var_es, ctrl_variant_es])

		plot(var, pdb_id, chain=None)

# main body
if __name__ == "__main__":
	main()







