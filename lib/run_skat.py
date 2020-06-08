import pandas as pd
import scipy as sp
from skat import *
import random
import pickle
import sys

def subsampling(scoreset, pheno, fold=10):

	# the original sample
	xy = np.unique(list(scoreset['ind1']) + list(scoreset['ind2'])) 
	mpheno = np.mean(pheno['status'])

	# individual number
	obs_pheno = pheno[ pheno['ind'].isin(xy) ]
	n_obs = len(xy)
	n_sample = fold * len(xy)
 
	n_case  = mpheno * n_sample - np.sum(obs_pheno['status'])
	n_case = int(n_case) if n_case > 0 else 0
	n_ctrl = n_sample - n_case - n_obs
	n_ctrl = int(n_ctrl) if n_ctrl > 0 else 0

		
	# for case not in observation
	re_case = pheno[ ~pheno['ind'].isin(xy) & pheno['status']==1 ]
	if re_case.shape[0] > n_case:
		case_sample = random.sample(re_case['ind'], n_case)
	else:
		case_sample = re_case['ind'].values.tolist()

	# for ctrl not in observation
	re_ctrl = pheno[ ~pheno['ind'].isin(xy) & pheno['status']==0 ] 
	if re_ctrl.shape[0] > n_ctrl:
		ctrl_sample = random.sample(re_ctrl['ind'], n_ctrl)
	else:
		ctrl_sample = re_ctrl['ind'].values.tolist()

	xy_new = np.unique(list(xy) + case_sample + ctrl_sample)	 
	dim = len(xy_new)

	xy_re = {x:i for i, x in enumerate(xy_new)}
	sc = scoreset['score'].values.tolist()
	x_re = map(lambda x: xy_re[x], scoreset['ind1'])	
	y_re = map(lambda x: xy_re[x], scoreset['ind2'])

	mat = kernel(x_re, y_re, sc, dim)

	subpheno = pheno[ pheno['ind'].isin(xy_new)].loc[:,'status'].values
	return mat, subpheno

def main():

        input_pk = sys.argv[1]
        pdbid    = sys.argv[2]

        input_combo = open(input_pk, 'rb')
        K     = pickle.load(input_combo)
        pheno = pickle.load(input_combo)
	cov   = pickle.load(input_combo)

	#determine if K is sparse or dense matrix
	if sp.sparse.issparse(K): K = K.toarray()

        # if there is only one element in the kernel
        # do not execute the calculation
        if K.shape[0] < 3 or len(pheno.unique()) < 2:
                sys.stdout.write('%s\tNA'%pdbid+ '\n')
                return None

        obj = SKAT(K, fixed_covariates=cov.values)
        pvals = obj.test(pheno.values, acc=1e-7)

        record = [pdbid, pvals]
        sys.stdout.write('\t'.join(map(str,record)) + '\n')

if __name__ == "__main__":

        main()

