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
	the beta weight: float or a 1D array
	"""
	beta_freq = stats.beta.pdf(freq, a, b)
	
	return np.array(beta_freq)

#######################################
### The score function(exponential form)
def exponential(dist_mat, t = 14):
	"""
	Given the euclidean distance of (variant_i, variant_j), 
	compute the similarity score between them.

	Args:
	dist_mat - Euclidean distance of variant_i, variant_j, 2D array

	Returns:
	struct_w - given the similarity function with high and low bound
	"""
	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	if min(dist_list) < 0:
		raise ValueError('invalid distance exist')
	else:
		struct_w_list = list(map(lambda r: np.exp(-r * r/(2 * t * t)), dist_list))

	struct_w = np.array(struct_w_list).reshape(num_row, num_col)
	return struct_w

#######################################
### The score function(with high and low bound)
def bounded(dist_mat, l=8, h=20):
	"""
	Given the euclidean distance of (variant_i, variant_j),
	compute the similarity score between them.

	Args:
	dist_mat - Euclidean distance of variant_i, variant_j, 2D array

	Returns:
	struct_w - given the similarity function with high and low bound
	"""

	num_row = dist_mat.shape[0]
	num_col = dist_mat.shape[1]

	dist_list = dist_mat.flatten()

	struct_w_list = map(lambda r: ValueError('invalid distance') if r < 0 else
	(1.0                                       if r <= l and r >= 0 else
	(0.5*math.cos((r - l)*np.pi/(h - l)) + 0.5 if r <= h and r > l
	else 0.0)),dist_list)

	struct_w = np.array(struct_w_list).reshape(num_row, num_col)

	return struct_w

#######################################
### The weight of variant
def weight_mat(freqs, distance_mat, aa_weight,use_aa=False,sim_fun='exponential',alpha = 0.5):
	"""
	Given a vector of minor allele frequency and matrix of cartisen coordinates
	compute the similarity matrix with score in the each cell representing
	the similarity between two individual 

	Args:
	freqs - minor allele frequency vector, a 1D array of size (num of variants)
	coordinates - a matrix with variants coordinates in a protein
			array of size ( num of variants , dimension) 
	sim_fun - 'exponential' or 'bounded': function to calculate pair-wise similarity
	alpha - parameter that determines the weights between structure kernal and frequency kernal
	rho - paramter that determines the weights between skat and burden test for frequency kernal	

	Returns:
	freq_struct_w - weights matrix for all variants
	"""

	# calculate frequency based kernel
	beta_freqs = beta(freqs)
	beta_freqs_diag = np.diag(beta_freqs)
	p = beta_freqs.shape[0]
	R  = np.diag(np.full(p,1)) 

	freq_w = np.linalg.multi_dot([beta_freqs_diag, R, beta_freqs_diag])
	beta_freqs_vec = beta_freqs.reshape((-1,1))
	freq_scale = beta_freqs_vec.dot(beta_freqs_vec.T)

	# calculate protein structure based kernal
	if sim_fun == 'exponential': struct_w = exponential(distance_mat)
	elif sim_fun == 'bounded': struct_w = bounded(distance_mat) 

	# scale by AA weight
	aa_scale = aa_weight.dot(aa_weight.T)
	if use_aa: struct_w *= aa_scale

	struct_norm = np.sum(struct_w)
	freq_norm = np.sum(freq_w)

	# combine two weights kernel with a parameter alpha
	if alpha == 0:
		freq_struct_w = struct_w / struct_norm
	elif alpha == 0.5:
		freq_struct_w = freq_scale * struct_w / (struct_norm * freq_norm)
	elif alpha == 1:
		freq_struct_w = freq_w / freq_norm	
			
	return freq_w, struct_w, freq_struct_w  

#######################################
def cal_aa_weight(snps2aa,pwm,use_pwm=False):
	"""
	Args:
	snps2aa - a pandas dataframe with header:
	  'snp','structid','chain','chain_seqid','x','y','z','aa'

	Returns:
	"""	
	snps2aa_unique = snps2aa[['varcode','aa','id']].drop_duplicates()
	n_snp = snps2aa['varcode'].nunique()

	if n_snp != snps2aa_unique.shape[0]:
		raise Exception('non-unique projection found')

	snps2aa_unique = snps2aa_unique.set_index('id')
	snps2aa_unique[['ref_aa','alt_aa']] = snps2aa_unique['aa'].str.split('/', expand=True)

	if use_pwm:
		log_score = snps2aa_unique.apply(lambda x: pwm.loc[x['ref_aa'],x['alt_aa']], axis=1)
		weight = log_score.apply(lambda x: np.exp(-x))
	else:
		weight = pd.Series(data=np.ones(len(snps2aa_unique)),
			index=snps2aa_unique.index.tolist())	

	aa_weight_mat = weight.values.reshape((n_snp,1))
	return np.sqrt(aa_weight_mat)

#######################################
### Construct non redundant distance matrix
def cal_distance_mat(snps2aa):
	"""
	Deal with situation when a single snp are mapped to more than two different residues 
	  in a protein structure. Maintain the smallest pairwise distance

	Args:
	snps2aa - a pandas dataframe with header:
	  'snp','structid','chain','chain_seqid','x','y','z'
	freqs 

	Returns:
	distance_mat - a matrix containing the minimum pairwise distance between snps
	"""
	
	n_snp = snps2aa['varcode'].nunique()
	snp_coord = snps2aa[['x','y','z']].values
	distance_vec = squareform(pdist(snp_coord)).flatten()
	
	snp_idx = snps2aa['id'].tolist()
	snp_pair_idx = list(map(list, it.product(snp_idx, repeat=2)))
	
	snp_pair_distance = pd.DataFrame(data=snp_pair_idx, columns=['i','j'], dtype=int)
	snp_pair_distance['r'] = distance_vec
	
	snp_pair_distance_uniq = snp_pair_distance.groupby(['i','j']).min().reset_index()
	
	row  = snp_pair_distance_uniq['i'].values
	col  = snp_pair_distance_uniq['j'].values
	data = snp_pair_distance_uniq['r'].values + 1.0
	
	distance_mat = sparse.coo_matrix((data,(row,col)),shape=(n_snp,n_snp)).toarray()
	np.fill_diagonal(distance_mat, 1.0)
	distance_mat[(distance_mat == 0.0)] = np.inf
	distance_mat -= 1.0

	# check if the distance matrix if symmetrical	
	#print((distance_mat.T == distance_mat).all())
	return distance_mat

def cal_Kernel(combined_w, genotype):
	"""
	calculate kernal from frequency and distance matrix

	Args:
	combined_w: weight matrix, describing pairwise distance/freqs between snps	
	genotype: genotype dataframe with row for ind and col for snps

	Returns:
	
	K - kernel that wraps the similarity score between individuals

	"""
	genotype_mat = sparse.csr_matrix(genotype.values)
	combined_w_sp = sparse.csr_matrix(combined_w)

	K = genotype_mat.dot(combined_w_sp)
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
	idx_tab.loc[:,'id'] = list(range(snps.shape[0]))
	idx_tab.loc[:,'snp'] = snps

	# get distance matrix
	distance_mat = cal_distance_mat(snps2aa, idx_tab)

	# generate the score matrix based on frequency and distance	
	freq, score, combined_w = score_mat(freqs.values, distance_mat)	
	K = DenseKernel(combined_w, genotype)

	obj = open('%s_K.pkl'%out,'wb')
	pickle.dump(K, obj)
	pickle.dump(pheno, obj)
	pickle.dump(cov, obj)
	obj.close()

	# only for ploting
#	obj = open('%s_mat.pkl'%out,'wb')
#       pickle.dump(freq, obj)
#       pickle.dump(score, obj)
#	pickle.dump(combined_w, obj)
#       obj.close()

# main body
if __name__ == "__main__":
	main()


