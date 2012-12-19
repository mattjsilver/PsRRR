# PATHWAYS SPARSE REDUCED-RANK REGRESSION
#
# Distributed under GNU general public licence (see COPYRIGHT.txt) 
# Copyright (C) 2012 Matt Silver
#
# WARRANTY DISCLAIMER
# Please note there is no warranty for the program, to the extent permitted by applicable law.
# Except when otherwise stated in writing the copyright holders and/or other parties provide the program
# "as is" without warranty of any kind, either expressed or implied, including, but not limited to, 
# the implied warranties of merchantability and fitness for a particular purpose. The entire risk
# as to the quality and performance of the program is with you. Should the program
# prove defective, you assume the cost of all necessary servicing, repair or correction.
#
# Please send comments or bugs to g.montana@imperial.ac.uk
#
#**************************************************************************************************
'''
select variables using lasso penalty
'''

import numpy as np

def lasso(predType=None, pred=None, resp=None, N=0, nPreds=0, params=None, computePath=False, nSNPsToSelect=10):
	'''
	lasso soft thresholding
	
	predType =  'b': sparse genotype selection (pred = X, resp = univariate Y)
			'a': sparse phenotype selection (pred = Y, resp = univariate X)
	pred = predictors
	resp = response
	(N,P) size of predictor matrix (X or Y)	
	
	computePath - compute full regularisation path? (currently uses scikit.learn)
	'''
	
	if computePath:
		pass
		return lasso_fullPath(predType=predType, pred=pred, resp=resp, N=N, nPreds=nPreds, 
						params=params, computePath=computePath, nSNPsToSelect=nSNPsToSelect)
	else:
		return lasso_pointWise(predType=predType, pred=pred, resp=resp, N=N, nPreds=nPreds, 
						params=params, computePath=computePath)


def lasso_fullPath(predType=None, pred=None, resp=None, N=0, nPreds=0, params=None, computePath=False, nSNPsToSelect=10):
	'''
	compute lasso over regularisation path using sklearn module
	** currently only implemented for simulation study
	nSNPsToSelect is min number of SNPs to select **
	'''
	
	from sklearn.linear_model import lasso_path
	class results:
		pass
	
	print 'running lasso estimation over regularisation path to select', nSNPsToSelect, 'SNPs'
	
	X = pred
	y = resp
	y.shape = (y.shape[0],)

	# Compute lasso path using sklearn.linear_model.lasso_path
	# produce (n_alphas x P) matrix of estimated coefficients
	# n_alphas is number of lambdas to compute (sklearn uses 'alpha' as reg param)
	# ensure that at least nSNPsToSelect SNPs are selected at the minimum lambda value
	maxNselSNPs = nSNPsToSelect - 1 
	eps = 0.5 # eps = lambda/maxLambda
	firstRun = True

	while maxNselSNPs < nSNPsToSelect:
		if not firstRun:
			print maxNselSNPs, 'SNPs selected, recomputing path with smaller minimum lambda value'
			firstRun = False
		eps /= 2.0 # decrease minimum lambda value
		models = lasso_path(X, y, eps=eps, n_alphas=100) # eps is min lambda as proportion of lambdaMax 
		coefs_lasso = np.array([model.coef_ for model in models])
		# compute number of selected SNPs at each lambda
		path_nSelSNPs = np.sum(coefs_lasso!=0,axis=1)
		# max no of SNPs selected
		maxNselSNPs = path_nSelSNPs[-1] # maximum number of SNPs selected (at minimum lambda)

	# find first point on reg path where at least maxNselSNPs are selected
	idx = np.where(path_nSelSNPs >= nSNPsToSelect)[0][0]
	# find nonzero coeffs at this point
	results.selPreds = np.where(coefs_lasso[idx,:]!=0)[0]

	return results


def lasso_pointWise(predType=None, pred=None, resp=None, N=0, nPreds=0, params=None, computePath=False, assume_orthogonalX=False):
	'''run lasso at single lambda - i.e. don't compute full path
		if assume_orthogonalX = True, then coordinate descent is skipped'''
	

	print 'running lasso estimation with lambda = ', params.Lambda

	if assume_orthogonalX:
		print 'assuming orthogonal design (no coordinate descent)'
	
	else:
		print 'assuming non-orthogonal design (run coordinate descent)'
		CDtol = 1e-5 # tolerance for CD convergence
	
	class results:
		pass
	
	
	if predType == 'b': # sparse genotype selection
		X = pred
		y = resp
		
		# get smallest lambda value for which no groups selected - initial upper bound for lambda
		lambdaMax = getLambdaMax(X, y, N, nPreds)
		Lambda = lambdaMax * params.Lambda
		
		if assume_orthogonalX:
			results.b = softThresh(np.dot(X.T,y), Lambda)
			results.selPreds = np.where(results.b != 0)[0]
		else:
			# run coordinate descent
			P = X.shape[1]
			b = np.zeros(P,) # initialise b
			r_full = y; r_full.shape = (y.shape[0],) # full residual

			CDconverged = False
			CDidx = 0
						
			while not CDconverged:
				
				b_last = b.copy()
				
				for j in range(P):
					Xjbj_init = np.dot(X[:,j].T,b[j]) # for fast update of full residual
					b[j] = softThresh(np.dot(X[:,j].T,r_full) + b[j],Lambda)
					r_full += Xjbj_init - np.dot(X[:,j].T,b[j]) # update full residual
					
				print 'CDidx:', CDidx,
				print 'selPreds:', b.nonzero()[0]
				
				# test convergence
				if CDidx: 
					bdotb = np.dot(b_last.T/np.linalg.norm(b_last),b/np.linalg.norm(b)) # change in b
					print 'bdotb',bdotb
					CDconvTest = 1 - bdotb
				
					print 'b convergence',  CDconvTest
					print
					
					if CDconvTest < CDtol:
						CDconverged = True
						print 'CD converged in', CDidx, 'iterations'
						print
						
				CDidx += 1
			
			results.selPreds = b.nonzero()[0]

		

	elif predType == 'a': # sparse phenotype selection
		
		Y = pred
		x = resp

		# get smallest lambda value for which no voxels selected - initial upper bound for lambda
		lambdaMax = getLambdaMax(Y, x, N, nPreds)
		Lambda = lambdaMax * 0.25 #params.Lambda
		
		results.a = softThresh(np.dot(Y.T,x), Lambda).T
		results.selPreds = np.where(results.a != 0)[0]
	
	return results
	

	
def getLambdaMax(pred,resp,N,nPreds):
	'''
	determine smallest lambda value at which no predictors are selected
	'''
	lambdaMax = abs(np.dot(pred[:,0].T,resp))
	for j in range(1,nPreds):
		Lambda = abs(np.dot(pred[:,j].T,resp))
		if Lambda > lambdaMax:
			lambdaMax = Lambda
	return lambdaMax



def softThresh(z,Lambda):
	'''
	lasso soft tresholding function
	z can be scalar or vector
	'''
	return np.sign(z) * np.maximum(0,abs(z) - Lambda)