from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from fastlmmclib.quadform.qfc_src import wrap_qfc
import numpy as np
import scipy as sp
from scipy import sparse
import sys

class VCT:
	def process_covariates(self, fixed_covariates=None):	
		
		# add bias term
		if fixed_covariates is None: self.X = np.ones((self.n, 1), dtype='float32')
		else:
			#scaler = StandardScaler() 
			#scaler.fit(fixed_covariates)
			self.X = fixed_covariates
			self.X = np.hstack([self.X, np.ones((self.n, 1))])

		self.X = self.X.astype('float32')
		self.p = self.X.shape[1]
		
		# Calculate X+
		self.Xdagger = np.linalg.pinv(self.X)
		self.Xdagger = self.Xdagger.astype('float32')

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
			coef = coef.reshape((1,-1))
			pred_p = np.dot(self.X, coef.T) 
			pheno_vec = pheno_vec - pred_p
		# estimated environmental variance under null hypothesis
		denonimators = np.dot(pheno_vec.T,pheno_vec)
		denonimators /= self.n - self.p

		# calculate part of the score
		nominators = self.K.dot(pheno_vec).T.dot(pheno_vec)
		return nominators / denonimators 

	def compute_p_value(self, r, acc):
		return self.davies(r, self.phis[np.where(self.phis > 1e-8)], acc)

	def davies(self, squaredform, eigvals, acc):
		sigma = 0.0
		lim=np.floor(1/acc)
		coeffs=eigvals
		chi2val=squaredform
	
		size = coeffs.shape[0]
		dof = np.ones(size,dtype = 'int32')
		noncentrality = np.zeros(size)
		ifault=np.zeros(1,dtype = 'int32')
		trace = np.zeros(7)
		ipval = 1.0-wrap_qfc.qf(coeffs,noncentrality,dof,sigma,chi2val,lim,acc,trace,ifault)

		if ifault != 0:
			return np.nan
		if ifault == 0 and ipval <= 0.0:
			return acc
		else:
			return ipval

	def test(self, phenotypes, acc=1e-8):

		scores = self.compute_scores(phenotypes)
		pvals = self.compute_p_value(scores, acc)
		return pvals

	def __init__(self, kernel_matrix=None, fixed_covariates=None, num_var=None,\
			zero_threshold=1e-8, phis=None):
		self.K = kernel_matrix.astype('float32')
		self.n = np.shape(self.K)[0]
		self.m = int(num_var)

		self.process_covariates(fixed_covariates=fixed_covariates)

		self.S =  np.identity(self.n, dtype='float32')
		
		XX_dagger = sp.linalg.blas.dgemm(alpha=1.0,a=self.X,b=self.Xdagger)
		XX_dagger = XX_dagger.astype('float32')
		XX_daggetT = sp.linalg.blas.dgemm(alpha=1.0,a=XX_dagger,b=XX_dagger,trans_b=True)
		self.S -= XX_daggetT
	
		self.SKS = self.K.dot(self.S).T
		self.SKS = self.SKS.astype('float32')
		self.SKS = sp.linalg.blas.dgemm(alpha=1.0,a=self.SKS,b=self.S)
			
		rank = min(self.n, self.m) - 1
		self.phis = sparse.linalg.eigsh(self.SKS, k=rank, return_eigenvectors=False)	
		self.phis = np.sort(self.phis)[::-1]
