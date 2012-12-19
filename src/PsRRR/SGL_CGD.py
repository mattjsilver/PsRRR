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
Run sparse group lasso estimation
using COORDINATE GRADIENT DESCENT
'''

from scipy.linalg import norm
from scipy.optimize import brentq
import numpy as np
from lasso import softThresh


def run_SGL_CGD(X=None,y=None,bPenalty=None,groups=None,N=None,P=None,gWeights=None,params=None):
	
	## RUN SGL ESTIMATION

	# get smallest lambda value for which no groups selected initial upper bound for lambda
	print 'running SGL estimation with lambda = ', params.Lambda, 'and alpha = ', params.alpha
	lambdaMax = getLambdaMax(X, y, groups, gWeights, params.alpha)

	Lambda = lambdaMax * params.Lambda

	# run estimation with specified penalty
	res = sgl(X=X, y=y, groups=groups, gWeights=gWeights, Lambda=Lambda, params=params)
	
	return res


def sgl(X=None, y=None, b=None, groups = None, gWeights=None, Lambda=None, params=None):
	'''
	solve SGL using coordinate gradient descent using Newton's method
	assume group independence
	
	INPUTS
	data      - data structure with expanded X,Y,groups, weights etc
	params    - various parameters, file locations etc
	Lambda    - regularization param
	params.alpha	  - distribution of regularisation between SNPs and pathways
	'''
	alpha = params.alpha
	nGroups = len(groups)

	# CD (coordinate descent within block) 
	showCDconvergence = False # verbose output for CD convergence
	CDtol = 1e-4 # tolerance for BCD convergence test
	maxCDits = 25 # maximum number of CD iterations

	# --------------------------------------------------------------------------
	# --------------------------------------------------------------------------
	## MAIN LOOP

	lSel = [] # selected groups in this BCD cycle
	CDits = [] # record CD iterations for each selected group
	allSelSNPs = [] # indices of all selected SNPs
	

	for l in range(nGroups):
		
		Xl = X[:,groups[l]] # group genotype matrix
		
		# check if group selected
		if np.linalg.norm(softThresh(np.dot(Xl.T,y), alpha*Lambda)) > (1-alpha)*gWeights[l]*Lambda:
			
			lSel.append(l) # all selected groups
			gSize = len(groups[l]) # size of current group
			theta = np.zeros((gSize,1)) # initialise group coeff vector
			weight_l = gWeights[l] # weight for this group

			if showCDconvergence:
				print '\nestimating group: ', l
			
			CDconverged = False
			CDidx = 0 # CD iteration counter for SNP j

			if showCDconvergence:
				blockObjF = SGLblockObjF(bl=theta, rl=y, Xl=Xl, weight=weight_l, Lambda=Lambda, alpha=alpha)
				print 'pre CGD block ObjF:', blockObjF

			
			# CD LOOP
			while not CDconverged:
				
# 				if showCDconvergence:
# 					print '\nCD pass:', CDidx
				
				lastTheta = theta.copy() 
				
				for j in range(gSize):
	
					
					if theta[j] == 0:
						# theta[j] = zero update:
						# test 1st deriv either side of zero to establish if theta[j] is at minimum
						# if not shift theta[j] accordingly
						theta[j] = sgl_j_zero_update(j=j,Xl=Xl,rl=y,theta=theta,weight=weight_l,Lambda=Lambda,alpha=alpha)
						
					
					else:
						
						theta[j] = sgl_j_nonzero_update(j=j,Xl=Xl,rl=y,theta=theta,weight=weight_l,alpha=alpha,Lambda=Lambda)
	
				# end j loop
				
				# test for CD convergence
				if not CDconverged and CDidx:
	
					# test for convergence of theta[j]
					CDconvTest = 1 - abs( np.dot(theta.T/norm(theta), lastTheta/norm(lastTheta)) )
					if showCDconvergence:
						print 'test CD convergence: 1-|bl.T last_bl| =', CDconvTest
	
					if CDconvTest < CDtol:
						CDconverged = True
						if showCDconvergence:
							print 'CD converged\n'
	
						
					if CDidx == maxCDits:
						CDconverged = True
						if showCDconvergence:
							print 'max number of CD iterations passed'
				
				CDidx += 1
		

			# CD converged for this block - record selected SNPs
			CDits.append(CDidx)
			selSNPsThisGroup = [groups[l][idx] for idx in theta.nonzero()[0]]
			allSelSNPs.extend(selSNPsThisGroup)

			if showCDconvergence:
				blockObjF = SGLblockObjF(bl=theta, rl=y, Xl=Xl, weight=weight_l, Lambda=Lambda, alpha=alpha)
				print 'post CGD block ObjF:', blockObjF
				print '# sel SNPs this group:', len(selSNPsThisGroup)
			
	
	
	# return relevant results
	class sglRes:
		pass
	sglRes.selGroups = lSel
	sglRes.selSNPs = list(set(allSelSNPs)) # count overlapping SNPs select multiple times once only

	# print some convergence info
	print '\nCGD iterations for each selected group:',CDits
	print 'selected groups:', sglRes.selGroups
	print len(sglRes.selSNPs),'SNPs selected'

	return sglRes



def sgl_j_zero_update(j=None,Xl=None,rl=None,theta=None,weight=None,Lambda=None,alpha=None):
	'''theta[j] update for selected group in the case where theta[j] = 0
	test 1st deriv either side of theta[j]=0
	if directional derivative in increasing both sides, theta[j] remains at zero
	otherwise do Newton update for theta[j] in specified direction
	return updated value of theta[j] (or 0 if no change)
	
	since group is selected, we assume norm(theta) != 0
	'''
	# no group penalty at theta[j] = 0
	# deriv in +ve direction
	d_pos = -np.dot(Xl[:,j].T,(rl - np.dot(Xl,theta))) + alpha * Lambda
	# deriv in -ve direction
	d_neg = -np.dot(Xl[:,j].T,(rl - np.dot(Xl,theta))) - alpha * Lambda
	
	# if d_pos and d_neg both < 0; shift theta[j] right => step = + 1
	# if d_pos and d_neg both > 0; shift theta[j] left  => step = - 1
	# if d_pos > 0 and d_neg both < 0 cannot decrease objF in either directio => step = 0
	step = -np.sign(np.sign(d_pos) + np.sign(d_neg))

	# update theta[j] depending on step; if step == 0, theta[j] remains at zero	
	if step:
		theta_j = 0
		theta_norm = norm(theta)
		if theta_norm == 0: # second deriv = 1
			nextTheta_j = np.dot(Xl[:,j].T,rl) - alpha * Lambda * step
		else: # second deriv = 1 + (1-alpha) * weight * Lambda / theta_norm
			nextTheta_j = (np.dot(Xl[:,j].T,(rl - np.dot(Xl,theta))) - alpha * Lambda * step)  \
				/ (1 + (1-alpha) * weight * Lambda / theta_norm)
		# check obj function is decreasing; if not halve step size
		newTheta = theta.copy()
		newTheta[j] = nextTheta_j
		newBlockObjF = SGLblockObjF(bl=newTheta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
		oldBlockObjF = SGLblockObjF(bl=theta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
		counter = 0 # keep track of number of times step size is halved
		while newBlockObjF > oldBlockObjF: # block obj f increasing, halve step size
			nextTheta_j = (nextTheta_j + theta_j)/2.
			newTheta[j] = nextTheta_j
			newBlockObjF = SGLblockObjF(bl=newTheta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
			if counter > 5:
			# first selected SNP often slow to converge; quicker on next CD pass
				oldBlockObjF = newBlockObjF
			counter += 1
		
		return nextTheta_j

	
	else: # objF increasing either side of zero, theta[j] remains at zero

		return 0
	


def sgl_j_nonzero_update(j=None,Xl=None,rl=None,theta=None,weight=None,alpha=None,Lambda=None):
	''' Newton update with theta[j] and norm(theta) nonzero'''
	theta_j = theta[j]
	theta_norm = norm(theta)
	
	# 1st deriv
	d1 = -np.dot(Xl[:,j].T,(rl - np.dot(Xl,theta))) \
		+ alpha * Lambda * np.sign(theta_j) \
		+ (1-alpha) * weight * Lambda * theta_j / theta_norm
	# 2nd deriv
	d2 =  1 + (1-alpha) * weight * Lambda / theta_norm * (1 - theta_j**2 / theta_norm**2)
	
	# Newton update
	nextTheta_j = theta_j - d1 / d2
	
	# check obj function is decreasing; if not halve step size
	newTheta = theta.copy()
	newTheta[j] = nextTheta_j
	newBlockObjF = SGLblockObjF(bl=newTheta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
	oldBlockObjF = SGLblockObjF(bl=theta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
	counter = 0 # keep track of number of times step size is halved
	while newBlockObjF > oldBlockObjF: # block obj f increasing, halve step size
		nextTheta_j = (nextTheta_j + theta_j)/2.
		newTheta[j] = nextTheta_j
		newBlockObjF = SGLblockObjF(bl=newTheta, rl=rl, Xl=Xl, weight=weight, Lambda=Lambda, alpha=alpha)
		if counter > 5:
			# first selected SNP often slow to converge; quicker on next CD pass
			oldBlockObjF = newBlockObjF
		counter += 1
	
	return nextTheta_j



def SGLglobalObjF(y=None, X=None, b=None, groups=None, gWeights=None, Lambda=None, alpha=None):
	'''
	compute objective function for single block only
	'''
	norms = map(norm,[b[groups[l]] for l in range(len(groups))])
	sumWeightedNorms = sum([gWeights[l] * norms[l] for l in range(len(groups))])
	return 0.5 * norm(y - np.dot(X, b))**2 + (1-alpha) * Lambda * sumWeightedNorms + Lambda * alpha * norm(b,1)


def SGLblockObjF(bl=None, rl=None, Xl=None, weight=None, Lambda=None, alpha=None):
	'''
	compute objective function for single block only
	'''
	return 0.5 * norm(rl - np.dot(Xl, bl))**2 + (1-alpha) * Lambda * weight * norm(bl) + Lambda * alpha * sum(abs(bl))


def getLambdaMax(X, y, groups, weights, alpha):
	'''
	determine minimum Lambda value at which no groups are selected
	'''
	
	lambdaMax = np.zeros((len(groups),))
	
	# need to find Lambda (>0) at which
	# norm(S(Xl'y,alpha*Lambda)) = (1-alpha)*Lambda
	# this is piecewise quadratic
	# to solve must first find upper bound for Lambda
	# at which soft thresholding function produces zero vector
	# this upper bound is given by: 
	# max(abs(np.dot(X[:,groups[l]].T,y)))
	# we can then find Lambda as the value between 0 and this upper bound

	for l in range(len(groups)):
		# determine upper bound for Lambda (see above)
		Xly = np.dot(X[:,groups[l]].T,y) # 2 x faster
		if alpha:
			maxLambda_l = max(abs(Xly)/alpha)
			lambdaMax[l] = brentq(lambdaMaxFunc,0,maxLambda_l,args=(Xly,weights[l],alpha))
		else:
			lambdaMax[l] = norm(Xly)/weights[l]
	
	return max(lambdaMax)

def adaptWeights_selGroup(X=None, y=None, groups=None, weights=None, alpha=None):
	'''
	UNIVARIATE PHENOTYPE ONLY
	determine first group to be selected when lambda is reduced from lambdaMax
	this is equivalent to 
	max_l: lambdaMax
	
	for multivariate phenotypes, full SGL estimation must be carried out
	in order to determine a vector during PsRRR
	'''
	
	# check univariate phenotype
	if y.shape[1] > 1:
		print 'fast weight tuning algorithm works only with univariate phenotype'
		print 'SGL weight tuning for multivariate phenotype not yet implemented!'
		print '...exiting'
		import sys; sys.exit()
	
	lambdaMax = np.zeros((len(groups),))
	
	# need to find Lambda (>0) at which
	# norm(S(Xl'y,alpha*Lambda)) = (1-alpha)*Lambda
	# this is piecewise quadratic
	# to solve must first find upper bound for Lambda
	# at which soft thresholding function produces zero vector
	# this upper bound is given by: 
	# max(abs(np.dot(X[:,groups[l]].T,y)))
	# we can then find Lambda as the value between 0 and this upper bound

	for l in range(len(groups)):
		Xly = np.dot(X[:,groups[l]].T,y) # 2 x faster
		# determine upper bound for Lambda (see above)
		maxLambda_l = max(abs(Xly)/alpha)
		lambdaMax[l] = brentq(lambdaMaxFunc,0,maxLambda_l,args=(Xly,weights[l],alpha))

	selGroup = np.argmax(lambdaMax)
	print '[SGL-CGD.adaptWeights_selGroup] selected group:',selGroup
	
	return selGroup

		
		
def lambdaMaxFunc(x, Xly, weight, alpha):
	'''piecewise quadratic function used to determine lambdaMax'''
	return norm(softThresh(Xly, alpha*x))**2 - ((1-alpha)*x*weight)**2