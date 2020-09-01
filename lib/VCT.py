from sklearn.linear_model import LogisticRegression
import fastlmm.util.stats.quadform as qf
import numpy as np
import scipy as sp
import sys

class VCT:
	def process_covariates(self, fixed_covariates=None):
		self.X = fixed_covariates

		# add bias term
		if self.X is None: self.X = np.ones((self.n, 1))
		else: self.X = np.hstack([self.X, np.ones((self.n, 1))])

		self.p = self.X.shape[1]
		
		# Calculate X+
		if self.X is not None: self.Xdagger = np.linalg.pinv(self.X)

	def compute_scores(self, phenotypes):
		pheno = phenotypes.reshape(self.n, -1)

		# determine whether it is binary or continuous variable
		# calculate the number of unique values in the phenotype
		uniq = np.unique(phenotypes).shape[0]
		flag = 1 if uniq > 2 else 0

		# for binary phenotype
		if flag:
			log_fit = LogisticRegression(fit_intercept=False).fit(self.X, pheno)
			coef    = log_fit.coef_
	
			pheno = np.dot(self.X, coef.T) 
		
		# estimated environmental variance under null hypothesis
		denonimators = np.linalg.multi_dot([pheno.T, self.S, pheno])
		denonimators /= self.n - self.p
 
		# calculate part of the score
		nominators = np.linalg.multi_dot([pheno.T, self.SKS, pheno])

		return nominators / denonimators 

	def compute_p_value(self, r, acc):
		return self.davies(r, self.phis[np.where(self.phis > 1e-8)], acc)

	def davies(self, squaredform, eigvals, acc):
		return qf.qf(squaredform, eigvals, acc=acc,lim=10000)[0]

	def test(self, phenotypes, acc=1e-6):
		scores = self.compute_scores(phenotypes)
		pvals = self.compute_p_value(scores, acc)
		return pvals

	def __init__(self, kernel_matrix=None, fixed_covariates=None,\
			zero_threshold=1e-8, phis=None):
		self.K = kernel_matrix
		self.n = np.shape(self.K)[0]
		self.process_covariates(fixed_covariates=fixed_covariates)

		self.S =  np.identity(self.n)
		self.S -= np.linalg.multi_dot([self.X, self.Xdagger,self.Xdagger.T, self.X.T])
		self.SKS = np.linalg.multi_dot([self.S, self.K, self.S])

		# Perform the eigendecomposition
		if phis is not None:self.phis = phis 
		else:
			self.phis = np.linalg.eigvalsh(self.SKS)
		
			# Round to zero
			self.phis[self.phis < zero_threshold] = 0.0
			self.phis = np.sort(self.phis)[::-1]
	
			eigen_mean = np.mean(self.phis)
			eigen_var  = np.var(self.phis)

