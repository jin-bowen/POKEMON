from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from scipy import stats
import numpy as np 
import pandas as pd
import itertools as it
import pickle
import math
import sys

#######################################
### The weight of variant
def beta(freq, a=1, b=25):
	"""
	Given a set of minor allele frequency, beta function parameter a and b, 
		compute the weight of variants
	Args:
	freq - minor allele frequency, float or a 1D array
	a - beta function param1, float
	b - beta function param2, float
	
	Returns:
	the weight 
	"""
	beta_freq = stats.beta.pdf(freq, a, b)
	
	return np.array(beta_freq)

#######################################
### The score function(exponential form)
def score_exp(dist_mat, t = 14):
	"""
	Given the euclidean distance of (variant_i, variant_j), 
	compute the similarity score between them.

	Args:
	dist_mat - Euclidean distance of variant_i, variant_j, 2D array

	Returns:
	score_mat - given the score function with high and low bound
	"""
	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	if min(dist_list) < 0:
		raise ValueError('invalid distance exist')
	else:
		score_list = list(map(lambda r: np.exp(-r * r/(2 * t * t)), dist_list))

	score_mat = np.array(score_list).reshape(num_row, num_col)
	return score_mat

#######################################
### The score function(with high and low bound)
def score_cut(dist_mat, l=8, h=20):
	"""
	Given the euclidean distance of (variant_i, variant_j),
	compute the similarity score between them.

	Args:
	dist_mat - Euclidean distance of variant_i, variant_j, 2D array

	Returns:
	score_mat - given the score function with high and low bound
	"""

	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	score_list = map(lambda r: ValueError('invalid distance') if r < 0 else
	(1.0                                       if r <= l and r >= 0 else
	(0.5*math.cos((r - l)*np.pi/(h - l)) + 0.5 if r <= h and r > l
	else 0.0)),dist_list)

	score_mat = np.array(score_list).reshape(num_row, num_col)

	return score_mat

#######################################
### The weight of variant
def score_mat(freqs, distance_mat, score_fun='exp', optimize='False',alpha = 0.0):
	"""
	Given a vector of minor allele frequency and matrix of cartisen coordinates
	compute the individual pair matrix with similarity score in the cells representing
	the similarity between two 

	Args:
	freqs - minor allele frequency vector, a 1D array of size (num_var)
	coordinates - a matrix with variants coordinates in a protein
			a matrix of size (num_var, dimension) 

	Returns:
	freq_score_mat - a score matrix for protein 
	"""

	# calculate frequency based kernel
	beta_freqs = beta(freqs)
	freqs_mat  = np.dot(beta_freqs, beta_freqs)

	# calculate protein structure based kernal
	score_mat = score_exp(distance_mat)

	# calculate kernel for burden test
	#burden_mat = 

	freq_score_mat = freqs_mat * score_mat
	# normalize score matrix
	norm = np.diagonal(freq_score_mat).max()

	return freqs_mat, score_mat, freq_score_mat * 1.0/norm

#######################################
### Construct non redundant distance matrix
def cal_distance_mat(snps2aa, freqs):
	"""
	Deal with situation when a single snp are mapped to more than two different residues 
	  in a protein structure. Maintain the smallest pairwise distance

	Args:
	snps2aa - a pandas dataframe with header:
	  'snp','structid','chain','chain_seqid','x','y','z'

	Returns:
	distance_mat - a matrix containing the minimum pairwise distance between snps
	"""

	snps = freqs.index.values
	idx_tab = pd.DataFrame()
	idx_tab['id'] = list(range(snps.shape[0]))
	idx_tab['snp'] = snps

	snps2aa_idx = pd.merge(snps2aa, idx_tab, on='snp')

	snp_coord = snps2aa_idx[['x','y','z']].values
	distance_vec = squareform(pdist(snp_coord)).flatten()

	snp_idx = snps2aa_idx['id'].tolist()
	snp_pair_idx = list(map(list, it.product(snp_idx, repeat=2)))

	snp_pair_distance = pd.DataFrame(data=snp_pair_idx, columns=['i','j'], dtype=int)
	snp_pair_distance['r'] = distance_vec

	snp_pair_distance_uniq = snp_pair_distance.groupby(['i','j']).min().reset_index()

	row  = snp_pair_distance_uniq['i'].values
	col  = snp_pair_distance_uniq['j'].values
	data = snp_pair_distance_uniq['r'].values

	distance_mat = sparse.coo_matrix((data,(row,col))).toarray()

	# check if the distance matrix if symmetrical	
	#print((distance_mat.T == distance_mat).all())
	
	return distance_mat

#def DenseKernel(snp_score_mat, freqs, IndVec, pheno):
#	"""
#	seperate common and rare variants
#	for common variants: a dense matrix
#	for rare variants: a sparse matrix
#	output kernel is the sum of above matrix
#
#	Args:
#	maf - minor allele frequency threshold. float 
#
#	Returns:
#	
#	base_K - kernel tha describe the similarity score between individuals
#
#	"""
#	maf = 0.1
#	dim = len(pheno)
#
#	base_K = np.array([[0 for i in range(dim)] for j in range(dim)])
#
#	rarevar_list = set(freqs.index.values)
#	rarevar_pdb = IndVec[ IndVec['gd_id'].isin(rarevar_list)].copy()
#	ind_rare = rarevar_pdb['ind'].unique()	
#
#	rare_mat = []
#
#	for iid, jid in it.combinations_with_replacement(ind_rare, 2):
#		
#		idx_iid = pheno.query('IID == "%s"'%iid).index.values[0]
#		idx_jid = pheno.query('IID == "%s"'%jid).index.values[0]
#		
#		igd = rarevar_pdb.query('ind == "%s"'%iid)['gd_id'].unique().tolist()
#		jgd = rarevar_pdb.query('ind == "%s"'%jid)['gd_id'].unique().tolist()
#
#		idx_igd = map(lambda x: freqs.index.tolist().index(x), igd)
#		idx_jgd = map(lambda x: freqs.index.tolist().index(x), jgd)
#
#		score_ij = 0
#		for si, sj in it.product(idx_igd, idx_jgd): 
#			score_ij += snp_score_mat[si, sj] 
#
#		rare_mat.append([score_ij,idx_iid, idx_jid])
#
#	rare_mat = np.array(rare_mat)
#	rare_K = csr_matrix((rare_mat[:,0], (rare_mat[:,1],rare_mat[:,2])),shape=(dim, dim)).toarray()
#	rare_K = rare_K + rare_K.T - np.diag(rare_K.diagonal())
#
#	K = rare_K + base_K
#	return K

def DenseKernel(snp_score_mat, genotype):
	"""

	calculate kernal from frequency and distance matrix

	Args:
	snp_score_mat: scored pairwise distance between snps	
	genotype: genotype dataframe with row for ind and col for snps

	Returns:
	
	K - kernel tha describe the similarity score between individuals

	"""
	genotype_mat = sparse.csr_matrix(genotype.values)
	snp_score_sp = sparse.csr_matrix(snp_score_mat)

	K = genotype_mat.dot(snp_score_sp)
	K = K.dot(genotype_mat.transpose())

	return K

def main():

	pdb_combo_file = sys.argv[1]
	# outfile path
	out  = sys.argv[2]

	pdb_combo = open(pdb_combo_file,'rb')
	genotype = pickle.load(pdb_combo)
	freqs    = pickle.load(pdb_combo) 
	pheno    = pickle.load(pdb_combo) 
	snps2aa  = pickle.load(pdb_combo) 
	cov      = pickle.load(pdb_combo)

	snps = freqs.index.values
	idx_tab = pd.DataFrame()
	idx_tab['id'] = list(range(snps.shape[0]))
	idx_tab['snp'] = snps

	# get distance matrix
	distance_mat = cal_distance_mat(snps2aa, idx_tab)

	# generate the score matrix based on frequency and distance	
	freq, score, snp_score_mat = score_mat(freqs.values, distance_mat)	
	K = DenseKernel(snp_score_mat, genotype)

	obj = open('%s_K.pkl'%out,'wb')
	pickle.dump(K, obj)
	pickle.dump(pheno, obj)
	pickle.dump(cov, obj)
	obj.close()

	# only for ploting
#	obj = open('%s_mat.pkl'%out,'wb')
#       pickle.dump(freq, obj)
#       pickle.dump(score, obj)
#	pickle.dump(snp_score_mat, obj)
#       obj.close()

# main body
if __name__ == "__main__":
	main()


