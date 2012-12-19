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
# If you use PsRRR in a published paper, please cite at least one of the following 2 papers:
#
# PATHWAYS ANALYSIS WITH UNIVARIATE Y
# Silver, Matt, and Montana, G. "Fast Identification of Biological Pathways Associated with
# a Quantitative Trait Using Group Lasso with Overlaps." Statistical Applications in
# Genetics and Molecular Biology 11, no. 1 (2012).
#
# PATHWAYS ANALYSIS WITH MULTIVARIATE Y
# (In Press) Silver, M., et al., Identification of gene pathways implicated in
# Alzheimer's disease using longitudinal imaging phenotypes with sparse regression, NeuroImage (2012). 
# Also available at http://arxiv.org/abs/1204.1937
#**************************************************************************************************
'''
pGLAW adapt weights algorithm
N/2 subsamples version
'''

import os, sys, datetime, cPickle
import numpy as np
import PsRRR, groupWeights, getData, setPathAndParams
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def adaptWeights(paramFile=None):

	print '------------------------------------------'
	print 'PsRRR weight tuning'
	print '------------------------------------------\n\n'
	
	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile)
	params.adaptWeights = True

	if params.mode == 'pp':
	
		# launch ppserver across all remote servers
		import getPP_job_server
		job_server = getPP_job_server.Server(params) # parallel python job_server instance
		
		params.nWorkers = job_server.nWorkers # number of available workers (i.e. CPUs)
		
		# check that number of desired iterations at every successive weight adaptation
		# will be divisible by the number of processors available
		if any([params.R[r] % params.nWorkers for r in range(len(params.R))]):
			print 'number of iterations:', params.R, 'not divisible by number of available processors:', params.nWorkers
			print 'please ammend pglawParams_xxx.py file'
			print 'exiting...'
			sys.exit()



	for iteration in xrange(params.firstIteration,len(params.R)):
		print '\n------------------------------------------'
		print 'adaptive weights iteration', iteration, '; eta = ', params.eta[iteration]
		print 'start time:', datetime.datetime.now() # current date and time
		if iteration == 0: # first run - use standard size weighting 
			weights = groupWeights.createWeights(params, stdWeights=True)
		else:
			# load current group weights (i.e. from previous iteration)
			print 'loading current group weights'
			weightsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/' + params.WEIGHTS_FILE + '_iteration' + str(iteration)
			weights = groupWeights.createWeights(params, stdWeights=False, weightsFn=weightsFn)
		
		# submit jobs
		selFreqThreads = []

		if params.mode == 'pp':
			# distribute generation of selection frequencies across workers
			print 'calculating empirical selection frequencies over', params.R[iteration], \
				'simulations distributed across', params.nWorkers, 'processes'
			print 'this may take some time...'
			sys.stdout.flush()

			for job in xrange(params.nWorkers):
				selFreqThreads.append(job_server.submit(
					func=fit_model,
					args=(iteration, job, params, weights),
					modules=("import numpy as np","getData","PsRRR","sys","datetime")))
			job_server.wait() # wait till all jobs completed

			# print some job_server stats
			print
			job_server.print_stats()

			print 'Iteration', str(iteration), 'completed'
			sys.stdout.flush()
		
		
		else: # no parallel processing - single job for testing

			selFreqThreads.append(fit_model(iteration, 1, params, weights))


		# generate new set of weights using these selfreqs		
		processWeights(iteration=iteration, selFreqThreads=selFreqThreads, oldWeightsName='iteration' 
			+ str(iteration), newWeightsName='iteration' + str(iteration + 1), params=params)

		# save plot of selFreq distribution for iteration just completed
		visualiseKLdivConvergence(iteration, iteration, params)		
		
		print 'iteration', iteration, 'completed'
		sys.stdout.flush()
		
	
	print 'all iterations completed'

	return
	
		

def fit_model(iteration, job, params, weights):
	
	# run GL over R model fits, each with different permuted (null) response
	# select single group at each fit; override various params as follows:
	params.B = 1 					# single estimation for each permuted phenotype
	params.Lambda = 0.999 			# lambda very close to lambdaMax - require 1 selected group
	
	if params.mode == 'pp':
		# stdout to file
		stdout = sys.stdout
		sys.stdout = open(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/selFreqs_' + params.WEIGHTS_FILE + '_iteration' + str(iteration) + '_thread' + str(job) + '.log','w')
		print 'execution time:', datetime.datetime.now() # current date and time
		print 'total iterations:', params.R[iteration], '; total threads:', params.nWorkers
		its_this_thread = params.R[iteration]/params.nWorkers
	
	else:
		params.nWorkers = 1
		its_this_thread = params.R[iteration]

	print '\n[fit_model] running iteration:', iteration, '; thread:', job
	print 'total null fits this thread:', its_this_thread
	
	# load complete dataset once and create subsample at each model fit
	# GL uses expanded dataset to account for overlaps
	if params.bPenalty == 'GL':
		data = getData.Data(params=params, subsample=False, standardise=True, useExpanded=True)
	elif params.bPenalty == 'SGL':
		data = getData.Data(params=params, subsample=False, standardise=True, useExpanded=False)
	else:
		print 'unrecognised genotype penalty... exiting'
		sys.exit()
	
	r = 0
	while r < (params.R[iteration] / params.nWorkers): # each thread performs R/nWorkers model fits
		print '\nnull fit', r + 1, 'of', its_this_thread
		
		# if using SGL penalty with unvariate phenotype, no need to fit model as we don't need to estimate b
		if params.bPenalty == 'SGL' and data.Y.shape[1] == 1:
			import SGL_CGD
			print 'SGL adaptWeights: permuting phenotype'
			y = np.random.permutation(data.Y)
			selGroup = SGL_CGD.adaptWeights_selGroup(X=data.X, y=y, groups=data.groups, weights=weights, alpha=params.alpha)
		else:
			# (Y permuted in PsRRR.run_PsRRR when params.adaptWeights = True)
			res = PsRRR.run_PsRRR(params=params, nSubsamps=1, weights=weights, data=data, 
				job=job, bPenalty=params.bPenalty, aPenalty=params.aPenalty)
			selGroup = res['selGroups'][0][0]
									
		
		# collate results
		if r == 0:
			selFreqs = np.zeros((len(data.groups), 1))
			selFreqs[selGroup] = 1
		else:
			selFreqs[selGroup] += 1
			
		r += 1

	selFreqs /= params.R[iteration] / params.nWorkers
	print 'group selection frequencies:'
	print selFreqs.T
	sys.stdout.flush()

	if params.mode == 'pp':
		sys.stdout = stdout
	
	return selFreqs



def processWeights(iteration, selFreqThreads, oldWeightsName, newWeightsName, params):
	'''
	load selFreq information generated using last set of weights
	and generate new set of weights
	'''	
	
	# load selFreqs for each thread
	if params.mode == 'pp':
		selFreqs = selFreqThreads[0]() # first thread selFreqs
	else:
		selFreqs = selFreqThreads[0] # first thread selFreqs

	for thread in range(1, params.nWorkers): # add remaining thread results
		if params.mode == 'pp':
			selFreqs += selFreqThreads[thread]()
		else:
			selFreqs += selFreqThreads[thread]

	selFreqs /= params.nWorkers # normalise


	# compute KL divergence
	Pi_l = len(selFreqs)**-1 # unbiased pathway selection frequency
	klDiv = computeKLdiv(selFreqs, Pi_l)
	print '\nKL divergence:', klDiv

	# compute number of pathways with zero selection frequency
	print 'pathways with zero sel freq:', sum(selFreqs==0)
	print

	# save cumulative selFreqs and klDiv
	resultsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/selFreqs_' + params.WEIGHTS_FILE + '_' + oldWeightsName	
	resultsFh = open(resultsFn + '_all.pickle', 'w')
	cPickle.dump(selFreqs, resultsFh)
	cPickle.dump(klDiv, resultsFh)
	resultsFh.close()	

	# --------------
	# ADAPT WEIGHTS
	
	# load old weights
	print 'adapting weights'
	oldWeightsFh = open(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/' + params.WEIGHTS_FILE + '_' + oldWeightsName + '.pickle')
	oldWeights = cPickle.load(oldWeightsFh)
	oldWeightsFh.close()
	d = selFreqs - Pi_l # diff between empirical and unbiased selFreq
	print 'calculating new weight distribution with eta = ', params.eta[iteration]
	# apply standard weight correction
	newWeights = oldWeights * (1-np.sign(d)*(params.eta[iteration]-1) * (len(selFreqs)**2 * d**2)); print 'square update'
	# newWeights = oldWeights * (1-(params.eta[iteration]-1) * (len(selFreqs) * d)); print 'linear update'


	# apply extra weight correction to groups with zero selection frequency over past 4 iterations
	# create nPaths x 4 matrix with columns representing selFreqs over previous 4 iterations (including current)
	if iteration > 3:
		
		allSelFreqs = selFreqs
		for idx in range(iteration-3,iteration):
			# load prev selFreqs
			resultsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/selFreqs_' + params.WEIGHTS_FILE + '_' + 'iteration' + str(idx)
			resultsFh = open(resultsFn + '_all.pickle')
			newSelFreqs = cPickle.load(resultsFh)
			resultsFh.close()
			allSelFreqs = np.concatenate((allSelFreqs,newSelFreqs),axis=1) # columnwise concatenation
		# compute paths with zero selFreq across all 4 previous iterations
		allSelFreqsZero = np.sum(allSelFreqs==0,axis=1) == 4
		# apply extra weight reduction to these pathways only
		newWeights[allSelFreqsZero] = oldWeights[allSelFreqsZero] * 0.95
		print 'extra weight adjustment 0.95 applied to', sum(allSelFreqsZero), 'pathways' 
	
	
	newWeightsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/' + params.WEIGHTS_FILE + '_' + newWeightsName + '.pickle'
	newWeightsFh = open(newWeightsFn,'w')
	cPickle.dump(newWeights, newWeightsFh, -1)
	newWeightsFh.close()
	print 'new weights saved to:', newWeightsFn
	

def computeKLdiv(selFreqs, Pi_l):
	'''
	compute KL divergence
	# selFreqs is observed pathway selection frequency distribution.  This must sum to 1
	# Pi_l is unbiased (uniform) distribution of selection frequencies; Pi must also sum to 1 and be same length as selFreqs
	'''
	
	nPaths = len(selFreqs)
	
	kl = 0
	for i in range(nPaths):
		if selFreqs[i]: # 0 log 0 = 0
			kl = kl + (selFreqs[i] * np.log(selFreqs[i]/Pi_l))
	
	return kl




def visualiseKLdivConvergence(firstIteration, lastIteration, params):
	'''
	visualise selection frequencies and KLdiv from firstIteration to lastIteration
	'''
	import cPickle
	
	for iteration in xrange(firstIteration, lastIteration + 1):
		# load results file
		resultsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/selFreqs_' + params.WEIGHTS_FILE + '_iteration' + str(iteration) + '_all.pickle'
		resultsFh = open(resultsFn)
		selFreqs = cPickle.load(resultsFh)
		klDiv = cPickle.load(resultsFh)
		resultsFh.close()
		
		# create fig
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.bar(range(len(selFreqs)), selFreqs, color='k')
		plt.xlim([0,len(selFreqs)])
		plt.ylim([0, min(1,4./len(selFreqs))]) # may want to change this - ylim is 4 x expected selFreq under the null 
		plt.xlabel('pathway #')
		plt.ylabel('selection frequency')
		# add some text
		txtLabel = 'R: ' + str(params.R[iteration]) + '\n' + \
					'KLdiv: ' + str(np.around(klDiv,decimals=2)) + '\n' + \
					'n(Pi_l = 0): ' + str(sum(selFreqs==0))
		ax.text(0.75,0.85,txtLabel,transform = ax.transAxes) # transform to axis coords
		if iteration == 0:
			plt.title('pathway selection frequencies: adaptive weights - standard size weighting')
		else:
			plt.title('pathway selection frequencies: adaptive weights - iteration ' + str(iteration))
#		plt.show()
		# save as pdf
		if not os.path.exists(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/plots/'):
			print 'making plots directory'
			os.makedirs(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/plots/')
		figFN = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/plots/selFreqs_iteration' + str(iteration) + '.pdf'
		plt.savefig(figFN,format='pdf')

	
if __name__ == "__main__":

	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()
	
	adaptWeights(paramFile=paramFile)
