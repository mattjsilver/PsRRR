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
Run PsRRR with specified penalty:
bPenalty = GL/SGL/lasso - sparse regression penalty on genotype coefficient vector
aPenalty = currently lasso only - sparse regression penalty on phenotype coefficient vector

Behaviour varies according to the dimensionality (Q) of the phenotype:
Q = 1 - (univariate phenotype) there is no sparse reduced rank regression loop
[see Silver and Montana, Fast Identification of Biological Pathways Associated with
# a Quantitative Trait Using Group Lasso with Overlaps. Statistical Applications in
Genetics and Molecular Biology 11, no. 1 (2012)].
Q > 1 - (multivariate phenotype) full sparse reduced rank regression
[see # Silver et al., Identification of gene pathways implicated in Alzheimer's disease
using longitudinal imaging phenotypes with sparse regression, NeuroImage (in press)]

There are two use cases:
A. Weight tuning prior to full analysis
PsRRR.run_PsRRR() called direct from adaptWeights.py, with data passed as arg
B. Full PsRRR analysis:
PsRRR.setup_PsRRR() is called from runPsRRR.py, data not passed as arg
The main reason for the difference is to do with computational efficiency and memory management

All methods run in one of two possible 'modes':
standard mode (params.mode = 'std)) - for running on single workstation
parallel python mode (params.mode = 'pp') - for parallel operation on cluster
'''

import os, sys, time
import numpy as np
import getData
import groupWeights


def setup_PsRRR(params,data=None,bPenalty='GL',aPenalty=None,subsample=True,nSNPsToSelect=None):
	'''
	Primary purpose is to set up PsRRR for parallel operation, although it can handle 
		the case of non-parallel operation as well
	'''
	
	# LOAD GROUPS WEIGHTS
	if bPenalty == 'GL' or bPenalty == 'SGL':
		if params.weights_iteration == 0: # standard size weighting
			if params.verbose:
				print 'generating group weights using standard size weighting'
			gWeights = groupWeights.createWeights(params, stdWeights=True)
		else: # use previously generated weights file
			weightsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/' + params.WEIGHTS_FILE + '_iteration' + str(params.weights_iteration) 
			print 'generating weights using', weightsFn
			gWeights = groupWeights.createWeights(params, stdWeights=False, weightsFn=weightsFn)
	else:
		gWeights = None

	allResults = []

	######################################
	####### PARALLEL OPERATION ###########
	# SET UP PARALLEL PYTHON OPERATION ON CLUSTER
	if params.mode == 'pp':
	
		# launch ppserver across all remote servers
		import getPP_job_server
		job_server = getPP_job_server.Server(params) # parallel python job_server instance
		
		params.nWorkers = job_server.nWorkers # number of available workers (i.e. CPUs)

		# check that number of desired subsamples to run is divisible by number of workers
		if params.B % params.nWorkers:
			print 'total number of subsamples:', params.B, 'not divisible by number of available worker:', params.nWorkers
			print 'please ammend pglawParams.py file'
			print 'exiting...'
			sys.exit()
	
		# submit jobs
		subsampsPerWorker = params.B/params.nWorkers
		print 'launching', params.nWorkers, 'workers', 'to run', params.B, 'subsamples'
		print 'each worker will fit', subsampsPerWorker, 'subsamples'
		sys.stdout.flush()

		print 'launching ppserver jobs...'
		print 'expect to wait a long time until completion!\n'
		
		for job in xrange(params.nWorkers):
		
			allResults.append(job_server.submit(
				func=run_PsRRR,
				args=(params,subsampsPerWorker,gWeights,data,job,bPenalty,aPenalty,subsample),
				modules=("sys","os","time","cPickle","import numpy as np","getData","GL")))
			
		job_server.wait() # wait till all jobs completed		
		
		
		# print some job_server stats
		print
		job_server.print_stats()

		# kill job_server
		job_server.destroy()

		print 'all jobs completed'
		print 'completion time:', time.strftime("%d_%b_%Y_%H_%M") 		# execution time


		return allResults
	
	######################################
	####### NON-PARALLEL OPERATION #######
	else: # no parallel processing - run single subsamp at a time

		for subsamp in range(params.B):
			print '\n*****************'
			print 'running subsamp:', subsamp + 1
			print '*****************'				
			allResults.append(run_PsRRR(params=params,nSubsamps=1,weights=gWeights,data=data,job=1,
									bPenalty=bPenalty,aPenalty=aPenalty,subsample=subsample, nSNPsToSelect=nSNPsToSelect))
		
		return allResults			
			
		
		

def run_PsRRR(params=None,nSubsamps=None,weights=None,data=None,job=None,
			bPenalty='GL',aPenalty=None,subsample=True,nSNPsToSelect=None):
	'''
	run PsRRR in one of two modes:
	1. adapting weights
		this function called directly from adaptWeights.py
		run PsRRR over X, Y:
			no subsampling
			no centering or scaling as this is already done by calling function
	2. not adapting weights
		create standardised and centered N/2 subsamp of X and Y before
		performing PsRRR 
	'''
	
	# containers for results
	if bPenalty == 'GL' or bPenalty == 'SGL':	
		allSelGroups = [] # selected groups at each subsample
	if bPenalty == 'lasso' or bPenalty == 'lasso_post_PsRRR' or bPenalty == 'SGL':
		allSelSNPs = [] # selected SNPs at each subsample
	if  aPenalty == 'lasso': # save selected voxels too
		selVoxels = []
	
	
	if params.mode == 'pp' and not params.adaptWeights: # log output to local file

		# make separate directory labelled with launch time for results files
		resultsDir = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + params.bPenalty + '_' + params.time
		if not os.path.exists(resultsDir):
			os.makedirs(resultsDir)
		# divert output to logfile
		stdout = sys.stdout
		sys.stdout = open(resultsDir + '/results_' + str(params.time) + '_job' + str(job) + '.log','w')

		print 'ALL PARAMS:'
		print params.__dict__
		
		print 'total subsamps:', params.B, '; total threads:', params.nWorkers
		print '\n[run_PsRRR] running job:', job
		sys.stdout.flush()


	for subsamp in xrange(nSubsamps):
		
		if not params.adaptWeights:
			print '\nrunning subsamp', subsamp + 1, 'of', nSubsamps
		
		if bPenalty == 'GL':
			print '\n\nUSING GROUP LASSO PENALTY'
		elif bPenalty == 'SGL':
			print '\n\nUSING SPARSE GROUP LASSO PENALTY'
		elif bPenalty == 'lasso':
			print '\n\nUSING LASSO PENALTY'
		elif bPenalty == 'lasso_post_PsRRR':
			print '\n\nUSING LASSO PENALTY ON PsRRR SELECTED GROUPS ONLY'
			

		######################################################################################################################
		# get all data
		
		if params.adaptWeights: 	# data passed in arg; just need to permute Y

			print 'WEIGHT TUNING: permuting (' + str(data.N) + ' x ' + str(data.Q) + ') phenotype matrix'
			X = data.eX
			P = data.eP
			Y = np.random.permutation(data.Y)

		else: 						# data not passed as arg

			if bPenalty == 'GL':
				data = getData.Data(params=params, subsample=subsample, standardise=True, useExpanded=True)
				X = data.eX
				P = data.eP

			elif bPenalty == 'SGL' or bPenalty == 'lasso':
									# lasso and SGL selection performed on the unexpanded genotype matrix
				data = getData.Data(params=params, subsample=subsample, standardise=True, useExpanded=False)
				X = data.X
				P = data.P
			
			elif bPenalty == 'lasso_post_PsRRR': 
									# doing post PsRRR analysis SNP selection with a reduced X corresponding to selected pathways
									# full standardised dataset passed in arg
				XsnpIDs = data.selSNPs # IDs of SNPs in X
				X = data.X
				P = data.P
			
			else:
				print 'no recognised genotype sparse regression penalty specified... exiting'
				sys.exit()

			if params.generateNull: # do pathway selection under the null, with correspondence between X and Y subjects broken
				print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
				print '!!!!!!!		WARNING		params.generateNull = True....			!!!!!!!!!!'
				print '!!!!!!!    	PERMUTING PHENOTYPE SUBJECTS TO GENERATE NULL DATA      !!!!!!'
				print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'				
				Y = np.random.permutation(data.Y)
			else:
				Y = data.Y			


		N = data.N		
		Q = data.Q

	
		######################################################################################################################
	
		######################################################################################################################
		# PsRRR
		# initialise a - phenotype coeff vector and b - genotype coeff vector
		a = np.ones((1,Q)); a = a/np.linalg.norm(a) # normalise
		b = np.zeros((P,1))
		PsRRR_idx = 0 # PsRRR loop iterations to convergence
	
		PsRRRconverged = False # flag to indicate convergence of PsRRR loop
		
		while not PsRRRconverged:
			
			if Q == 1: # univariate phenotype
				
				print 'univariate phenotype - not running PsRRR'
				PsRRRconverged = True
				
			else: # multivariate phenotype
				
				if bPenalty == 'SGL':
					print 'SGL not yet implemented for multivariate phenotype... exiting'
					sys.exit()
				if not params.adaptWeights:
					print '\n\n**********\nPsRRR iteration:', PsRRR_idx
			sys.stdout.flush()
			
			# current values of a and b to test for convergence
			last_a = a; last_b = b
			
			
			########################
			# UPDATE B
			if bPenalty == 'GL':
				import GL
				res_b = GL.run_GL(X=X, y=np.dot(Y,a.T), bPenalty=bPenalty, groups = data.eGroups, N=N, P=P, gWeights=weights, params=params)
			
			elif bPenalty =='SGL':
				import SGL_CGD
				res_b = SGL_CGD.run_SGL_CGD(X=X, y=np.dot(Y,a.T), bPenalty=bPenalty, groups=data.groups, N=N, P=P, gWeights=weights, params=params)
				selSNPs = res_b.selSNPs
				print 'SGL estimation completed'
			
			elif bPenalty == 'lasso_post_PsRRR':
				# lasso - update b using Y*a' in place of Y to give sparse estimate
				import lasso
				res_b = lasso.lasso(predType='b',pred=X,resp=np.dot(Y,a.T), N=N, nPreds=P, params=params)
				# map selected SNPs back to original SNP IDs
				selSNPs = [XsnpIDs[idx] for idx in res_b.selPreds]
				print 'selected SNPs:',selSNPs
				
			elif bPenalty == 'lasso':
				# lasso - update b using Y*a' in place of Y to give sparse estimate
				import lasso
				res_b = lasso.lasso(predType='b',pred=X,resp=np.dot(Y,a.T), N=N, nPreds=P, params=params, computePath=False, nSNPsToSelect=nSNPsToSelect)
				selSNPs = res_b.selPreds
				
			else:
				print 'no recognised genotype sparse regression penalty specified... exiting'
				sys.exit()
	
	
			########################
			# UPDATE A if multivariate phenotype, otherwise not required
			if Q > 1:
				# normalise b
				b = res_b.b / np.linalg.norm(res_b.b)
				if aPenalty == None:
					a = np.dot(np.dot(b.T, X.T), Y) / (np.dot(np.dot(b.T, X.T), np.dot(X, b)))
				else:
					# lasso - update a using b'X in place of X to give sparse estimate
					import lasso
					res_a = lasso.lasso(predType='a',pred=Y,resp=np.dot(X,b), N=N, nPreds=Q, params=params)
					# record selected SNPs - note these are unexpanded SNP ids!
					print 'n selected voxels:', len(res_a.selPreds)
					a = res_a.sparseCoeffvec
	
				a = a/np.linalg.norm(a) # normalise
	
	
				########################
				# TEST FOR CONVERGENCE
				if PsRRR_idx: # don't test on first iteration
					PsRRRconvTest = 1 - abs( np.dot (b.T/np.linalg.norm(b, 2), last_b/np.linalg.norm(last_b,2) ) )
					if params.verbose:
						print '\nPsRRR CONVERGENCE:'
						print 'b vector convergence test (1-|b\'b''|):', PsRRRconvTest
						print 'a vector convergence:  (1-|a\'a''|):', 1 - abs( np.dot (a/np.linalg.norm(a, 2), last_a.T/np.linalg.norm(last_a,2) ) ), '\n\n'

					if PsRRR_idx > params.max_PsRRR_iterations:
						print 'Maximum PsRRR cd iterations exceeded!'
						print '... terminating cd using current selected pathways'
						PsRRRconverged = True
					elif PsRRRconvTest < params.PsRRRconvergenceThresh:
						PsRRRconverged = True
						print '\nPsRRR converged in', PsRRR_idx, 'iterations'

				else:

					last_b = 1
				
					
				PsRRR_idx += 1

		
		# record selected groups/SNPs etc for this subsamp
		if bPenalty == 'GL' or bPenalty == 'SGL':
			allSelGroups.append(res_b.selGroups)
		if bPenalty == 'lasso' or bPenalty == 'lasso_post_PsRRR' or bPenalty == 'SGL':
			allSelSNPs.append(selSNPs)
		if  aPenalty == 'lasso': # save selected voxels too
			selVoxels.append(res_a.selPreds)
		
	# end for subsamp in range(nSubsamps)
	
	sys.stdout.flush()
	if params.mode == 'pp' and not params.adaptWeights:
		sys.stdout = stdout
		
				
	########################
	#### SAVE/RETURN RESULTS
	if params.adaptWeights:
		return {'selGroups': allSelGroups, 'nGroups': len(weights)}
	elif bPenalty == 'SGL':
		return {'selGroups': allSelGroups, 'selSNPs': allSelSNPs}
	elif bPenalty == 'GL':
		return {'selGroups': allSelGroups}
	elif bPenalty == 'lasso' or bPenalty == 'lasso_post_PsRRR':
		return {'selSNPs': allSelSNPs}
	elif  aPenalty == 'lasso': # save selected voxels too
		return {'selGroups': allSelGroups, 'selVoxels': selVoxels}



		
#  	import interactiveShell as ipy; ipy.ipshell() # start ipython shell for debugging