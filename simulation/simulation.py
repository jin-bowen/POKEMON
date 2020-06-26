from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from scipy import stats

import pymol
from pymol import cmd

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

	Returns:
	"""
	snp_corr_es = snp_corr_es_df.values

	snp_es = np.sqrt(np.diagonal(snp_corr_es))	
	np.fill_diagonal(snp_corr_es,0)

	ind_genotype = genotype.values
	ind_log_odds =  np.dot(ind_genotype, snp_es)
	
	ind_ind_corr = np.linalg.multi_dot([ind_genotype, snp_corr_es, ind_genotype.T])
	ind_log_odds += np.diagonal(ind_ind_corr)

	ind_prob = np.exp(ind_log_odds) / (1 + np.exp(ind_log_odds))
	phenotype = list(map(lambda prob: np.random.binomial(1, p=prob), list(ind_prob)))

	phenotype_df = pd.DataFrame(phenotype, columns=['phenotype'])

	return phenotype_df

def exponential_es(var_val, dist_mat_df, t=14):
	"""
	Args:

	Returns:
	"""
	snps = var_val.index
	dist_mat = dist_mat_df.values

	#print(np.all(var_val.index == dist_mat_df.index))

	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	if min(dist_list) < 0:
		raise ValueError('invalid distance exist')
	else:
		struct_w_list = list(map(lambda r: np.exp(-r * r/(2 * t * t)), dist_list))

	struct_w = np.array(struct_w_list).reshape(num_row, num_col)
	snp_es = var_val['es'].values

	snp_corr = np.dot(snp_es, snp_es.T)
	snp_es_combined = snp_corr * struct_w
	
	snp_es_combined_df = pd.DataFrame(snp_es_combined, \
					index=snps,\
					columns=snps)
	return snp_es_combined_df

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
def cal_distance_mat(snps2aa_tot, freqs):
	"""
	Args:
	
	Returns:
	"""	
	snps = freqs.index.values
	n_snp = snps.shape[0]
	
	idx_tab = pd.DataFrame()
	idx_tab['id'] = list(range(n_snp))
	idx_tab['snp'] = snps
	
	snps2aa_tot_idx = pd.merge(snps2aa_tot, idx_tab, on='snp')
	snps2aa_grp = snps2aa_tot_idx.groupby(['structure'])
	dist_mat_dict = {}

	for key, snps2aa in snps2aa_grp:

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
	
		distance_mat = sparse.coo_matrix((data,(row,col)),shape=(n_snp,n_snp)).toarray()
		distance_mat[distance_mat == 0.0] = np.inf
		distance_mat_df = pd.DataFrame(distance_mat,columns=snps,index=snps)

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
		resi  = row['structure_position']
		chain = row['chain']
		pheno = row['es']
		cmd.alter('resi %s and chain %s'%(resi, chain), 'b=%s'%pheno)
	max_es = max(var['es'].values)
	cmd.spectrum("b", "white_red", "%s"%pdb_id, maximum=max_es, minimum=0)
	cmd.zoom()

def main():

	case_var_file = sys.argv[1]	
	ctrl_var_file = sys.argv[2]	
	reference = sys.argv[3]

	gene_name = case_var_file.split('_')[0]

	cols=['chr','pos','if','ref','alt', 'qual', 'filter', 'info', 'format']
	case_var = pd.read_csv(case_var_file, sep="\t", \
			header=0, names=cols, usecols=list(range(9))) 
	ctrl_var = pd.read_csv(ctrl_var_file, sep="\t", \
			header=0, names=cols, usecols=list(range(9))) 

	case_var['snp'] = case_var['chr'].astype(str) + ":" + \
				case_var['pos'].astype(str) + ":" + \
				case_var['ref'].astype(str) + ":" + \
				case_var['alt'].astype(str)

	ctrl_var['snp'] = ctrl_var['chr'].astype(str) + ":" + \
				ctrl_var['pos'].astype(str) + ":" + \
				ctrl_var['ref'].astype(str) + ":" + \
				ctrl_var['alt'].astype(str) 

	case_var_val = pd.DataFrame(index=case_var['snp'].values)
	case_var_val['es'] = 0.1
	case_var_val['freq'] = 0.01
	ctrl_var_val = pd.DataFrame(index=ctrl_var['snp'].values)
	ctrl_var_val['es'] = 0.00
	ctrl_var_val['freq'] = 0.01
	var_val = pd.concat([case_var_val, ctrl_var_val])

	num_ind = 2000
	num_var = var_val.shape[0]

	snp      =  var_val.index
	var_es   = var_val['es'].values
	var_freq = var_val['freq'].values

	# generate genotype
	ind_genotype = np.zeros((num_ind, num_var))	

	for i, ifreq in enumerate(var_freq):
		num_case = np.random.binomial(num_ind, ifreq)
		case_id = np.random.choice(num_ind, num_case)
		ind_genotype[case_id, i] = 1
	
	ind_genotype_df = pd.DataFrame(ind_genotype, columns=snp)

	# mapped snp
	ref_raw = dd.read_csv(reference, sep="\t")
	ref = ref_raw.groupby(['transcript'])
	snps2aa = snps_to_aa(snp, gene_name, ref)
	dist_mat_dict = cal_distance_mat(snps2aa, var_val)
	
	for pdb, dist_mat in dist_mat_dict.items():

		var_corr_es = exponential_es(var_val, dist_mat, t=7)
	
		phenotype = generate_phenotype(ind_genotype_df, var_corr_es)
		pdb_snps2aa = snps2aa.loc[snps2aa['structure'] == pdb]
		var_aa = pd.merge(var_val, pdb_snps2aa, left_index=True, right_on='snp')
		plot(var_aa, pdb, chain=None)

# main body
if __name__ == "__main__":
	main()

