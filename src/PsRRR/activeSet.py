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
classes and methods related to creating and checking active set
'''
import numpy as np

class ASData(object):
	'''
	data object corresponding to current AS only
	ASdata.AS			- groups in active set
	ASdata.X			- reduced genotype matrix corresponding to groups in active set only
	ASdata.b			- associated SNP coefficient vector
	ASdata.groups		- group -> SNP idx mappings corresponding to groups in AS only
	ASdata.gWeights		- AS group weights
	ASdata.skipEstimation (boolean)	- True: skip GL estimation (AS too large or too small)
	ASdata.computeAS (boolean)		- True: recompute AS (AS too small) 
	'''

	def __init__(self,X=None,y=None,groups=None,Lambda=None,gWeights=None,params=None,penalty=None,lambdaMax=None):

		self.Y = y
		
		if params.verbose:
			print(' computing active set\n')

		self.AS = findAS(X=X, y=y, groups=groups, Lambda=Lambda, weights=gWeights, penalty=penalty, params=params)
		
		# when adapting weights, recompute AS with larger lambda value if more than one group in AS
		if params.adaptWeights:
			while len(self.AS) > 1:
				print 'len(AS) =', len(self.AS), 'recomputing AS with smaller lambda'
				Lambda = (lambdaMax + Lambda) * 0.5
				self.AS = findAS(X=X, y=y, groups=groups, Lambda=Lambda, weights=gWeights, penalty=penalty, params=params)
		
		generateAS(self, X, groups, gWeights)
			


def findAS(X=None, y=None, groups=None, Lambda=None, weights=None, penalty=None, params=None):
	'''
	determine active set - i.e. selected groups at current Lambda
	'''
	AS = [g for g in xrange(len(groups)) if np.linalg.norm(np.dot(X[:,groups[g]].T, y)) >= Lambda * weights[g]]
		
	return AS



def checkAS(X=None,y=None,b=None,groups=None,Lambda=None,gWeights=None,penalty=None,params=None):
	'''
	check for violations of kkt conditions to see if current AS needs updating
	'''
	
	AS = []
	fullResidual = y - np.dot(X, b)

	for l in range(len(groups)):
		if np.linalg.norm(np.dot(X[:,groups[l]].T, (fullResidual + np.dot(X[:,groups[l]], b[groups[l]])))) >= Lambda * gWeights[l]:
			AS.append(l)

	return AS



def generateAS(self,X,groups,gWeights):
	'''
	generate reduced data set (X,groups,weights) corresponding to groups in current AS
	(self.AS) only	
	'''
	# 1. design matrix
	snpAS = [snp for group in self.AS for snp in groups[group]] # SNPs in AS
	self.X = X[:,snpAS] # submatrix containing only SNPs in AS				

	# 2. corresponding group mapping
	self.groups = []				
	gSizes = [len(groups[group]) for group in self.AS] # AS group sizes
	groupStart = 0
	for g in range(len(self.AS)):
		self.groups.append(range(groupStart,groupStart + gSizes[g]))
		groupStart += gSizes[g]
		
	# 3. weights
	self.gWeights = [gWeights[i] for i in self.AS]
	
	# 4. b (coeff vector)
	self.b = np.zeros((len(snpAS),1)) # initialise coeff vector
	self.computeAS = False # no need to compute AS again, unless lambda is reduced