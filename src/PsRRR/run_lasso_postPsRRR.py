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
Run lasso SNP selection on selected pathways
This is used for highlighting SNPs and genes that are driving
the selection of top ranked pathways.  We do this by using the 
existing PsRRR code, but with a lasso, rather than a group lasso
penalty on b.  The estimation is performed on a modified X, containing
only SNPs belonging to the selected pathways
'''

import getPhenotypes, PsRRR, setPathAndParams
import cPickle, sys
import numpy as np


def runLasso(paramFile):

	print '\n\n------------------------------------------'
	print 'sRRR for SNP/gene selection'
	print '------------------------------------------\n\n'


	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile)

	# get program parameters, file locations etc
	params.adaptWeights = False
	params.simulateResponse = False
	params.mode = 'std' # not yet implemented on cluster	
	params.Lambda = 0.8

	# run single PsRRR with lasso selection on set of pathways selected at each subsample in previous analysis
	# load all selected pathways
	resultsFn = "" # pathway results text file
	fn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + resultsFn
	print 'processing pathway selection results from: ' + fn
	fh = open(fn)
	allSelPaths = []
	for line in fh.readlines():
		allSelPaths.append(eval(line.rstrip()))
	fh.close()
	
	subsampsProcessed = len(allSelPaths)
	print 'selPaths for', subsampsProcessed, 'loaded'
	
	results = []
	params.B = 1 # run PsRRR once only - no subsampling
	nSNPs = [] # number of SNPs included in model at each model
	nSelSNPs = [] # number of SNPs selected at each model fit


	for subsamp in range(subsampsProcessed):
		print '\n[run_PsRRR_Lasso] running subsamp', subsamp+1, 'of', subsampsProcessed
		selPaths = allSelPaths[subsamp]
	
		# get data, including reduced design Matrix, X
		data = ReducedData(selPaths,params)
		nSNPs.append(data.P) # number of cols in X
		res = PsRRR.setup_PsRRR(params=params,data=data,bPenalty='lasso_post_PsRRR',subsample = False)
		nSelSNPs.append(len(res[0]['selSNPs'][0]))
		results.append(res)

	# save lasso SNP selection results to file for processing
	print 'all subsamps completed...'
	
	fn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + resultsFn + '_SNPselFreqs_' + params.time + '.pickle'
	fh = open(fn, 'w')
	cPickle.dump(results,fh,-1)
	fh.close()
	
	print 'SNP selection results saved to file:', fn
	
	report_SNP_stats(nSNPs,nSelSNPs) # report SNP selection statistics

		

class ReducedData(object):
	
	def __init__(self,selPaths,params):
	
		self.selPaths = selPaths

		# get original, unexpanded X
		fh = open('_'.join(params.DATA_FILE.split('_')[0:-1]) + '_genotypes.pickle')
		data = cPickle.load(fh)
		fh.close()
		self.X = data.X # not expanded!
		self.N = data.N
		
		# get unexpanded pathway map - this maps SNP indices (cols of X) to groups
		fh = open('_'.join(params.DATA_FILE.split('_')[0:-1]) + '_pathwayMapping.pickle')
		data = cPickle.load(fh)
		fh.close()
		self.eGroups = data.groups # not expanded!
		
		
		# get phenotypes
		phenotypes = getPhenotypes.Phenotypes(params)
		self.Y = phenotypes.Y
		self.Q = phenotypes.Q
		
		# generate list of selected pathways !!!!! for simulation only !!!!!
		if params.simulateResponse:
			from random import sample
			print 'simulating 20 selected pathways above threshold, including actual supported pathway'
			print '!!!!!!!! FOR SIMULATION ONLY !!!!!!!!'
			selPaths = (sample(self.eGroups,19))
			selPaths.append(phenotypes.suppPath[0][0])
	
		# create reduced X
		# self.selSNPs is list of original, unexpanded SNP indices to be included
		self.selSNPs = [snp for group in selPaths for snp in self.eGroups[group]]
		self.selSNPs = list(set(self.selSNPs)) # remove overlaps
		self.X = self.X[:,self.selSNPs]
		self.P = len(self.selSNPs)
		
		print '(', self.N, 'x', self.P, ')', 'reduced X created with SNPs from selected pathways:', selPaths
	
	
		print 'scaling and centering X and Y....',
		N = self.Y.shape[0] # sample size

		# GENOTYPES
		# here we convert dtype to float64 *for subsample data only*
		if self.X.dtype != 'float64':
			self.X = self.X.astype(np.float64) # convert to float64 if not already
			self.X -=  np.mean(self.X,axis=0) # center predictors
			self.X *=  1/np.std(self.X,axis=0)/np.sqrt(N) # standardise predictors sum(Xj**2)=1
		
		# PHENOTYPES
		if self.Y.dtype != 'float64':
			self.Y = self.Y.astype(np.float64) # convert to float64 if not already
		self.Y -= np.mean(self.Y,axis=0) # center response					
		print ' [done]'
	

def report_SNP_stats(nSNPs,nSelSNPs):
	'''REPORT MODEL SIZE STATISTICS ACROSS ALL SUBSAMPS'''
	import numpy as np
	print 'mean SNPs persubsamp:', np.mean(nSNPs)
	print 'min SNPs persubsamp:', min(nSNPs)
	print 'max SNPs persubsamp:', max(nSNPs)
	print 'SD SNPs persubsamp:', np.std(nSNPs)			

	print 'mean sel SNPs persubsamp:', np.mean(nSelSNPs)
	print 'min sel SNPs persubsamp:', min(nSelSNPs)
	print 'max sel SNPs persubsamp:', max(nSelSNPs)
	print 'SD sel SNPs persubsamp:', np.std(nSelSNPs)			



if __name__ == '__main__':
	
	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()
	
	runLasso(paramFile)
	
#  	import interactiveShell as ipy; ipy.ipshell() # start ipython shell