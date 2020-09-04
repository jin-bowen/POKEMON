from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
import fastlmm.util.stats.quadform as qf
import numpy as np
import scipy as sp
from scipy import sparse
import sys

class VCT:
	def process_covariates(self, fixed_covariates=None):
		self.X = fixed_covariates

		# add bias term
		if self.X is None: self.X = np.ones((self.n, 1))
		else: self.X = np.hstack([self.X, np.ones((self.n, 1))])
		self.p = self.X.shape[1]
		
		# Calculate X+
		self.Xdagger = np.linalg.pinv(self.X)

	def compute_scores(self, phenotypes):
		pheno_vec = phenotypes.reshape(self.n, -1)

		# determine whether it is binary or continuous variable
		# calculate the number of unique values in the pheno_vectype
		uniq = np.unique(phenotypes).shape[0]

		# for binary pheno_vectype
		if uniq == 2:
			log_fit = LogisticRegression(fit_intercept=False).fit(self.X,phenotypes)
			coef    = log_fit.coef_	
			logit_p = np.dot(self.X, coef.T) 
			p_case_null = np.exp(logit_p)/(1 + np.exp(logit_p))
			pheno_vec = pheno_vec - p_case_null
		else:
			log_fit = LinearRegression(fit_intercept=False).fit(self.X,phenotypes)
			coef    = log_fit.coef_
			pred_p = np.dot(self.X, coef.T) 
			pheno_vec = pheno_vec - pred_p

		# estimated environmental variance under null hypothesis
		denonimators = np.dot(pheno_vec.T,pheno_vec)
		denonimators /= self.n - self.p

		# calculate part of the score
		nominators = self.K.dot(pheno_vec).T.dot(pheno_vec)
		return nominators / denonimators 

	def compute_p_value(self, r, acc):
		return self.davies(r, self.phis[np.where(self.phis > 1e-20)], acc)

	def davies(self, squaredform, eigvals, acc):
		#print(squaredform, eigvals)
		#print(qf.qf(squaredform, eigvals, acc=acc,lim=10000))
		return qf.qf(squaredform, eigvals, acc=acc,lim=10000)[0]

	def test(self, phenotypes, acc=1e-6):

		scores = self.compute_scores(phenotypes)
		pvals = self.compute_p_value(scores, acc)
		return pvals

	def __init__(self, kernel_matrix=None, fixed_covariates=None, num_var=None,\
			zero_threshold=1e-8, phis=None):
		self.K = kernel_matrix
		self.n = np.shape(self.K)[0]
		self.m = int(num_var)
		self.process_covariates(fixed_covariates=fixed_covariates)

		self.S =  np.identity(self.n)
		self.S -= np.linalg.multi_dot([self.X, self.Xdagger,self.Xdagger.T, self.X.T])
		self.SKS = sparse.csr_matrix(self.K.dot(self.S).T)	
		self.SKS = self.SKS.dot(self.S)

		rank = min(self.n, self.m) - 1
		self.phis = sparse.linalg.eigsh(self.SKS, k=rank, return_eigenvectors=False)	
		self.phis = np.sort(self.phis)[::-1]	

