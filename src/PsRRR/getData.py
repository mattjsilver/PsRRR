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
Generate data instance - X, Y, groups etc
with or without overlap expansion, scaling, subsampling etc
(specified by args)
'''

import cPickle,sys
from random import sample
from math import floor
import numpy as np
import getPhenotypes

class Data(object):
	'''
	Data class containing genotype, phenotype and group information
	With overlap expansion (default):
	eX			- (N x eP) expanded genotype design matrix
	eGroups 	- dictionary of expanded SNP indices
	Y			- (N x Q) phenotype matrix
	Without overlap expansion (lasso and SGL penalties):
	X			- (N x P) unexpanded genotype design matrix
	Y, Q as above
	'''

	def __init__(self,params=None,subsample=False,standardise=False,useExpanded=True):
		'''
		logical args:
		subsample 	- create N/2 subsample?
		standardise - centre and scale data?
		useExpanded - use expanded genotype matrix to account for overlaps? (default is true)
		'''
		
		self.loadData(params,subsample,useExpanded)
		
		# center and scale data
		if standardise:

			print 'scaling and centering X and Y....',
			N = self.Y.shape[0] # sample size

			# GENOTYPES
			# here we convert dtype to float64 *for subsample data only* (save memory)
			if useExpanded:
				if self.eX.dtype != 'float64':
					self.eX = self.eX.astype(np.float64) # convert to float64 if not already
					self.eX -=  np.mean(self.eX,axis=0) # center predictors
					self.eX *=  1/np.std(self.eX,axis=0)/np.sqrt(N) # standardise predictors sum(Xj**2)=1
			else:
				if self.X.dtype != 'float64':
					self.X = self.X.astype(np.float64) # convert to float64 if not already
					self.X -=  np.mean(self.X,axis=0) # center predictors
					self.X *=  1/np.std(self.X,axis=0)/np.sqrt(N) # standardise predictors sum(Xj**2)=1
			
			# PHENOTYPES
			if self.Y.dtype != 'float64':
				self.Y = self.Y.astype(np.float64) # convert to float64 if not already
			self.Y -= np.mean(self.Y,axis=0) # center response					
			print ' [done]'
			
		

	def loadData(self,params,subsample,useExpanded):
		
		# load existing genotype, and pathway mapping data

		print 'importing data (genotypes, phenotypes and groups) from', params.DATA_FILE
		
		dataFh = open(params.DATA_FILE)
		data = cPickle.load(dataFh)
		dataFh.close()
		self.N = data.N
		self.P = data.P
		self.eGroups = data.eGroups
		self.groups = data.groups
		print 'SNP mappings for ' + str(len(self.eGroups)) + ' pathways loaded'
				

		# load phenotype
		
		phenotypes = getPhenotypes.Phenotypes(params)
		self.Q = phenotypes.Q
			
		
		# check sample size same for genotype and phenotype
		if phenotypes.N != self.N:
			print 'sample size for phenotypes and genotypes don''t match!  Exiting...'
			sys.exit()
			
		# generate subsample if required
		if subsample: # create N/2 bootstrap samples

			# genotype QC filtering ensures that there are no monomorphic SNPs
			# i.e. sum(X_j) = 0 for all j; therefore we can be certain that 
			# sum(X[:,j]>=1) for all j.  However, it is still possible that 
			# we create monomorphic SNPs in a single N/2 subsample.  
			# This will cause divide by zero problems later on
			# (e.g. with standardising predictors).  We therefore check for
			# monomorphic SNPs and resample where necessary
			
			monomorphicSNPs = True
			
			while monomorphicSNPs:

				sampIdx = createSubsamp(params, data.N)
				
				# load genotypes
				if useExpanded: # (default)
					print 'loading EXPANDED design matrix'
					self.eX = data.eX[sampIdx,]

					# check for monomorphic SNPs - resample if there are any
					if not np.any(np.sum(self.eX,axis=0)==0):
						monomorphicSNPs = False
					else:
						print 'monomoprhic SNPs detected - resampling!'

					self.eP = data.eP

				else:
					print 'loading UNEXPANDED design matrix'
					self.X = data.X[sampIdx,]

					# check for monomorphic SNPs - resample if there are any
					if not np.any(np.sum(self.X,axis=0)==0):
						monomorphicSNPs = False
					else:
						print 'monomoprhic SNPs detected - resampling!'

					self.P = data.P
	
				
				
			self.Y = phenotypes.Y[sampIdx,]
			self.N = len(sampIdx)

			if useExpanded:
				print '('+str(self.eX.shape[0]) + ',' + str(self.eX.shape[1]) + ') genotype ',
			else:
				print '('+str(self.X.shape[0]) + ',' + str(self.X.shape[1]) + ') genotype ',
				
			print 'and ' + '('+str(self.Y.shape[0]) + ',' + str(self.Y.shape[1]) + ') phenotype matrices generated'
						
		else: # no subsample (adaptWeights)

			# load genotypes
			if useExpanded: # (default)
				self.eX = data.eX
				self.eP = data.eP
			else:
				print 'loading UNEXPANDED design matrix'
				self.X = data.X
				self.P = data.P
	
			self.Y = phenotypes.Y
			print ' [done]'			
	

		
def createSubsamp(params, N=[]):
	'''
	create N/2 subsamp - balanced or unbalanced
	'''
	if params.balanced: # create N/2 bootstrap samples 
		# ensure subjects are balanced w.r.t. disease status
		# *** this is hard coded for the N=464 ADNI dataset at the moment ***
		# *** will need to be adapted for another dataset ***
		# load subject statuses

		try:
			status = np.genfromtxt(params.PsRRR_ROOT_DIR + 'data/' + params.STATUS_FILE,dtype=np.int8)
		except:
			print 'generating balanced subsample (params.balanced = True)'
			print 'failed to load status file (see params.STATUS_FILE)'
			sys.exit()
			
		# get number of cases and controls
		CN = np.where(status==0)[0] # controls have status '0'
		CS = np.where(status==1)[0] # cases have status '1'
		nCN = len(CN); nCS = len(CS)
		print 'generating N/2 balanced subsample for',nCS,'cases and',nCN,'controls'
		return np.concatenate((sample(CN,int(floor(N/4.))),sample(CS,int(floor(N/4.)))))			

		
	else:
		print 'generating N/2 subsample...',
		# create N/2 subsample
		return sample(xrange(N),int(floor(N/2)))