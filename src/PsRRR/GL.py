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
Group lasso estimation, using block coordinate descent
'''

import activeSet
import numpy as np

def run_GL(X=None,y=None,bPenalty=None,groups=None,N=None,P=None,gWeights=None,params=None):
	
	# get smallest lambda value for which no groups selected initial upper bound for lambda
	if not params.adaptWeights:
		print 'running GL estimation with lambda = ', params.Lambda		
	lambdaMax = getLambdaMax(X, y, groups, gWeights)
	Lambda = lambdaMax * params.Lambda
#
#	# compute AS (note ASdata.AS will always contain one group only when adapting weights)
	ASdata = activeSet.ASData(X=X, y=y, groups=groups, Lambda=Lambda, 
					gWeights=gWeights, params=params, penalty=bPenalty,lambdaMax=lambdaMax)
	if params.verbose:
		print 'lambda/lambda_max', Lambda/lambdaMax, '; AS size:', len(ASdata.AS)

			
	# run estimation
	kktOK = False # flag to indicate no kkt violations

	while not kktOK: # continue until no kkt violations

		# run estimation with specified penalty
		glRes = glBCD(data=ASdata, Lambda=Lambda, params=params)
		
		selGroups = [ASdata.AS[g] for g in glRes.selGroupsAS] # current selected groups
		if params.verbose:
			print '\nselected groups:', selGroups

		# generate coefficient vector, b, in original parameter space
		b = np.zeros((P,1)) # initialise coeff vector
		for g in range(len(selGroups)):
			# update complete coeff vector from coeffs in sel groups
			b[groups[selGroups[g]]] = glRes.bAS[ASdata.groups[glRes.selGroupsAS[g]]]
		
		if params.adaptWeights:
			kktOK = True # save time when adapting weights
		else:
			# check if any groups not currently in AS violate KKT
			# if so run estimation again with expanded set
			AScheck = activeSet.checkAS(X=X, y=y, b=b, groups=groups, Lambda=Lambda, gWeights=gWeights, penalty=bPenalty, params=params)
			# check for extra groups in AS
			extraGroups = [g for g in AScheck if g not in ASdata.AS]
			if not extraGroups:
				kktOK = True # leave loop
				if params.verbose:
					print 'no KKT violations\n'
			else:
				if params.verbose:
					print '****KKT violations****'
					print 'expanding active set'
				# append extraGroups to AS
				[ASdata.AS.append(g) for g in extraGroups]
				# update ASdata
				activeSet.generateAS(ASdata,X, groups, gWeights)
		
		
	# return results
	class results:
		pass
	results.selGroups = selGroups # pathway selection frequencies
	results.b = b # final coefficient vector
	return results


def getLambdaMax(X, y, groups, weights):
	'''
	determine minimum Lambda value at which no groups are selected
	'''

	lambdaMax = np.zeros((len(groups),))
	
	for g in range(len(groups)):
		lambdaMax[g] = np.linalg.norm(np.dot(X[:,groups[g]].T, y)) / weights[g]
			
	return max(lambdaMax)


def glBCD(data, Lambda, params):
	'''
	GL estimation using block coordinate descent with 
	local Taylor approximation
		
	INPUTS
	data      - data structure with expanded X,y,groups, weights etc
	params    - various parameters
	Lambda    - regularization param
	'''
	
	# throughout this function, 'CD' refers to coordinate descent within a SINGLE
	#	block, whereas 'BCD' refers to coordinate descent across multiple blocks
	
	## set up data
	X = data.X
	y = data.Y
	b = data.b.copy()
	groups = data.groups
	gWeights = data.gWeights
	nGroups = len(groups)
	
	## CONVERGENCE TESTING
	# BCD (block coordinate descent)
	BCDconverged = False # test for BCD convergence
	showBCDconvergence = False # verbose output for BCD convergence
	BCDtol = 1e-3 # tolerance for BCD convergence test
	maxBCDits = 40 # maximum number of BCD iterations	

	# CD (coordinate descent within block) 
	showCDconvergence = False # verbose output for CD convergence
	maxCDits = 500 # maxiumum number of CD iterations


	# --------------------------------------------------------------------------
	# --------------------------------------------------------------------------
	## MAIN LOOP

	BCD_idx = 0 # BCD iteration counter
	r = y # initial total residual (b = 0)

	if showBCDconvergence:
		print 'initial global obj F:', glObjF(y, X, b, groups, gWeights, Lambda)

	while not BCDconverged: # BCD

		lSel = [] # selected groups in this BCD cycle
		lastb = b.copy() # used for BCD convergence test

		if showBCDconvergence:
			CDits = [] # keep record of number of CD iterations for ref
			
		if showBCDconvergence:
			print '\nBCD pass: ', BCD_idx
		
		for l in range(nGroups):

			if showCDconvergence:
				print '\nestimating group: ', l
			
			gSize = len(groups[l]) # size of current group
			Xl = X[:,groups[l]] 
			rl = r + np.dot(Xl, b[groups[l]]) # compute partial residual, rl from full residual, r
			
			if np.linalg.norm(np.dot(Xl.T, rl)) <= Lambda * gWeights[l]:

				b[groups[l]] = 0
				
			else: # beta_l not zero - CD within block with Taylor approximation

				lSel.append(l)

				if np.linalg.norm(b[groups[l]]) == 0:
					bl = np.ones((gSize,1)) * 1e-6 # avoid divide by zero
				else:
					bl = b[groups[l]].copy()

				CDconverged = False # check for CD convergence within block
				CDidx = 0 # CD iteration counter

				# block objF before CD
				blockF = blockObjF(bl, gWeights[l], rl, Xl, Lambda)
				if showCDconvergence:
					print 'pre-CD block obj f:', blockF
					
					
				while not CDconverged: # CD
					
					preF_block = blockF.copy()
					
					# vectorised version
					bl = (np.dot(Xl.T, (rl - np.dot(Xl, bl))) + bl)/(1 + (Lambda * gWeights[l] / np.linalg.norm(bl)))
					
					# TEST CD CONVERGENCE
					if CDidx: # don't test on first iteration

						blockF = blockObjF(bl, gWeights[l], rl, Xl, Lambda)
						delta_f = preF_block - blockF
						if showCDconvergence:
							print 'CD idx:', CDidx, 'post-update block obj f:', blockF, 'decrease in block objF:', delta_f 
								
						if delta_f < 0: # good measure of CD convergence
							CDconverged = True
							if showCDconvergence:
								print 'CD converged\n'
							
						if CDidx == maxCDits:
							# CD convergence slow
							CDconverged = True
							if showCDconvergence:
								print 'max number of CD iterations passed'

					CDidx += 1

				b[groups[l]] = bl.copy() # update b for this group

				if showBCDconvergence: # keep a record of number of CD iterations to convergence
					CDits.append(CDidx)
				
				
				
			# update full residual
			r = rl - np.dot(Xl, b[groups[l]])
		
			# END FOR l in range(nGroups):
			
	


		# TEST GLOBAL CONVERGENCE
		if BCD_idx: # skip first BCD iteration

			BCDconvTest = 1-abs( np.dot (b.T/np.linalg.norm(b, 2), lastb/np.linalg.norm(lastb,2) ) )
			
			if BCDconvTest < BCDtol or BCD_idx > maxBCDits:
				BCDconverged = True 
			
			if showBCDconvergence:
				print 'BCD CONVERGENCE TEST:'
				print 'sel groups (AS):', lSel
				print 'No CD iterations (all selected groups):', CDits
				print 'global obj F:', glObjF(y, X, b, groups, gWeights, Lambda)
				print 'test global convergence: 1-|b.T lastb| =', BCDconvTest

				
			if BCDconvTest < BCDtol:
				BCDconverged = True
			
				if showBCDconvergence:
					print '\n**BCD CONVERGED**'
					print
	
			if BCD_idx > maxBCDits: # max number of BCD iterations exceeded
				BCDconverged = True
				if showBCDconvergence:
					print 'max number of CD iterations passed'
					print
	
	
		BCD_idx +=1
	# END WHILE NOT BCD CONVERGED
				
	
	# RETURN RESULTS
	class glRes:
		pass
	glRes.selGroupsAS = lSel
	glRes.bAS = b
	return glRes


def glObjF(y, X, b, groups, gWeights, Lambda):
	'''
	compute objective function (all blocks)
	'''
	#   objective function for complete model with all groups
	
	# first determine weighted sum of group l2 norms
	sumGpL2Norms = sum([np.linalg.norm(b[groups[g]]) * gWeights[g] for g in range(len(groups))])
	f = 0.5 * np.linalg.norm(y - np.dot(X, b))**2 + Lambda * sumGpL2Norms
	return f


def blockObjF(b_l, gweight, r_l, Xl, Lambda):
	'''
	compute objective function for single block only
	'''
	fblock = 0.5 * np.linalg.norm(r_l - np.dot(Xl, b_l))**2 + Lambda * gweight * np.linalg.norm(b_l)
	return fblock


def blockObjF_1d(theta_j, gweight, theta_j_, w_j, Z_j, Lambda):
	'''
	objective function for single group
	used to find theta_j using 1d optimisation
	
	INPUTS
	theta_j     - what we want to estimate
	gweight     - group weight
	theta_j_    - current estimate of bl, excl theta_j
	w_j         - current single predictor partial residual 
	                 (i.e. residual using current estimate for b, excluding j)
	Z_j         - predictor j genotype vector
	Lambda  	- regularization param
	'''
	# create group parameter vector
	bl = np.append(theta_j_, theta_j) # order doesn't matter

	f = 0.5 * np.linalg.norm(w_j - np.dot(Z_j, theta_j))**2 + Lambda * gweight * np.linalg.norm(bl)
	return f			

#  	import interactiveShell as ipy; ipy.ipshell() # start ipython shell