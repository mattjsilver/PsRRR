'''
Created on Mar 14, 2012
simulate data: X,Y and groups
@author: mattsilver
'''
import numpy as np, cPickle

class SimXandGroups(object):
	
	def __init__(self,params):
		# generate random allele counts from multinomial distribution
		# random (n x p) design matrix generated from multinomial distribution with allele
		# counts following HW distribution (A is major allele): 
		#   maf=0:pAA = pA^2
		#   maf=1:pAa = 2*pA*(1-pA)
		#   maf=2:paa = (1-pA)^2
		# where pA is frequency of major allele
		# n = sample size
		# p = number of SNPs
		# X = (n x p) design matrix

		from scipy.stats import rv_discrete
		
		# simulate X
		print 'generating X from multinomial with HW frequencies'
		if params.constMAF:
			pA = np.array([0.2] * params.P) # maf = const, 0.2
			print 'constant MAF = 0.2'
		else:
			pA = np.random.uniform(0.1,0.5,params.P) # 0.1 <= maf <= .5
			print '0.1 <= MAF <= .5'
			
		pVec = np.array([(1-pA)**2, 2*pA*(1-pA), pA**2])
		self.X = np.zeros((params.N,params.P),dtype='int')

		for j in range(params.P):
			vals = [range(3),pVec[:,j]]
			z=rv_discrete(name='adhoc',values=vals)
			self.X[:,j] = z.rvs(size=params.N)
		
		self.N = params.N; self.P = params.P
		
		if not params.overlap:
			self.eX = self.X; self.eP = params.P

		# simulate groups
		gpSize = params.groupSize # nSNPs per group
		self.groups = {}
		if params.overlap:
			print 'simulating', params.nGroups, '**overlapping** groups of size', gpSize
			for group in range(params.nGroups): # note this assumes groups overlap 50%
				self.groups[group] = range(group * (gpSize - params.overlapSize),gpSize + group*(gpSize - params.overlapSize))
#				print group,self.groups[group]
			# generate eX and eGroups
			import pathwaysSNPs,getDesignMatrix
			pathwaysSNPs.PathSNPs.expand_Pathways.__func__(self) # generates self.eGroups, self.eP and self.orig2exp
			getDesignMatrix.DesignMatrix.expandDesignMatrix.__func__(self,self.orig2exp) # generates self.eX
		
		else:

			print 'simulating', params.nGroups, 'non-overlapping groups of size', gpSize
			for group in range(params.nGroups):
				self.groups[group] = range(gpSize * group, group*gpSize + gpSize)
			self.eGroups = self.groups # assuming no overlaps
		
		# save X and groups
		dataFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/' + 'simData.pickle'
		dataFh = open(dataFn,'w')
		cPickle.dump(self, dataFh, protocol=-1)
		print '(' + str(self.N) + ',' + str(self.eP) + ') genotype matrix and groups saved to:',\
		dataFn + '\n'
		dataFh.close()
		
		params.DATA_FILE = dataFn



class SimY(object):

	def __init__(self,params=None,data=None):
		'''
		generate (N x Q) phenotype multivariate normal vector/matrix
		'''
		
		from random import sample

		# load X and groups
		dataFh = open(params.DATA_FILE)
		data = cPickle.load(dataFh)
		print '(' + str(data.N) + ',' + str(data.eP) + ') genotype matrix and groups loaded from:',\
		params.DATA_FILE + '\n'
		dataFh.close()


		print 'SIMULATING (' + str(params.N) + ',' + str(params.Q) + ') PHENOTYPE sampled from multivariate normal distbn'
		self.N = params.N; 
		self.Q = params.Q
		self.P = params.P
		
		# create multivariate normal baseline response
		if params.verbose:
			print 'generating multivariate normal baseline response'

		baseY = np.random.multivariate_normal(np.random.normal(10,1,self.Q),np.eye(self.Q),self.N)
		self.Y = baseY.copy() # for later calculation of variance explained
		
		
		# GENERATE PHENOTYPE (from unexpanded X)
		
		# determine causal SNPs / affected phenotypes
		self.nSNPs = params.nSNPs
		self.effSize = params.effSize # overall effect size
		# DELETE AS APPROPRIATE!
		if params.overlap:
			
			# CAUSAL SNPS DISTRIBUTED FROM ADJACENT (OVERLAPPING) PATHWAYS:		
			self.suppPath = sample(range(len(data.groups)-1),1) # first causal pathway
			self.suppPath.append(self.suppPath[0]+1) # add adjacent overlapping pathway
			suppPathSNPs = [snp for group in self.suppPath for snp in data.groups[group]]
			self.suppSNPs = sample(suppPathSNPs,params.nSNPs) # causal SNPs
			# draw SNPs so that they overlap two pathways only
#			self.suppSNPs = sample(suppPathSNPs[params.overlapSize:len(suppPathSNPs)-params.overlapSize],params.nSNPs) # causal SNPs
			if params.overlap:
				# supported paths
				self.suppPathsExp = []
				for group in range(len(data.groups)):
					for snp in self.suppSNPs:
						if snp in data.groups[group]:
							self.suppPathsExp.append(group)
				self.suppPathsExp = list(set(self.suppPathsExp))
				# supported SNPs in expanded var space
				self.suppSNPsExp = data.orig2exp[self.suppSNPs,:].nonzero()[1]
		else:
			# CAUSAL SNPS DISTRIBUTED IN SINGLE PATHWAY:		
			self.suppPath = sample(range(len(data.groups)), 1) # causal pathway
			self.suppSNPs = sample(data.groups[self.suppPath[0]],params.nSNPs) # causal SNPs
				
			# CAUSAL SNPS RANDOMLY DISTRIBUTED:
			# self.suppSNPs = sample(range(data.P),params.nSNPs) # causal SNPs

		# print MAF
		print 'ave MAF:', np.sum(data.X[:,self.suppSNPs])/float(2*self.N)/params.nSNPs
		self.suppPheno = sample(range(params.Q), params.nCausalPheno); # causal phenotypes
		self.v = np.zeros((1,params.nCausalPheno)) # variance explained by genetic effects for each affected phenotype

		w_q = np.zeros((self.N,params.nCausalPheno)) # used to calculate SNR
		for q_idx in range(params.nCausalPheno):

			q = self.suppPheno[q_idx] # affected phenotype idx
			
			# allelic effect size
			delta_q = self.effSize * self.N * np.mean(self.Y[:,q]) /   np.sum(data.X[:,self.suppSNPs])      
			
			# change in y_q due to additive genetic effects
			w_q[:,q_idx] = delta_q * np.sum(data.X[:,self.suppSNPs],axis=1)
			
			self.Y[:,q] += w_q[:,q_idx]

			# percentage of variance explained			
			self.v[0,q_idx] = np.var(w_q[:,q_idx])/(np.var(baseY[:,q]) + np.var(w_q[:,q_idx]))
		
		print '--------'
		print 'mean effect size [mean(w_q)/mean(baseY)]:', np.sum(np.mean(w_q,axis=0))/np.sum(np.mean(baseY,axis=0))
		print 'mean SNR:', (1./self.Q) * np.sum( np.var(w_q,axis=0) / np.var(baseY[:,self.suppPheno],axis=0) )
		print 'mean % var explained:', np.mean(self.v)
		print '--------'
				
		print 'suppPath:', self.suppPath
		if params.overlap:
			print 'overlapping paths:', self.suppPathsExp

		self.suppSNPs.sort()
		print 'suppSNPs:', self.suppSNPs
		if params.overlap:
			print 'SNPs in expanded var space:', self.suppSNPsExp

	
		self.PHENOTYPE_FILE = params.DATA_FILE.split('.pickle')[0] + '_pheno.pickle'
		phenoFh = open(self.PHENOTYPE_FILE,'w')
		cPickle.dump(self, phenoFh, protocol=-1)
		if params.verbose:
			print '(' + str(self.N) + ',' + str(self.Q) + ') phenotype matrix generated and saved to:',\
			self.PHENOTYPE_FILE + '\n'
		phenoFh.close()
		# ensure this simulated phenotype is used from now on
		params.PHENOTYPE_FILE = self.PHENOTYPE_FILE

		
def findMatchPaths(groups, p, suppSNPs):
	## GENERATE SNP TO PATHWAY LOOKUP TABLE
	# Need this to check for pathways with same support
	# snpPathMap is p x nPaths matrix 
	# with a 1 denoting SNP i occurs in pathway j
	nGroups = len(groups)
	snpPathMap = np.zeros((p,nGroups));
	for path in range(nGroups):
		snpPathMap[groups[path],path] = 1
	reducedMap = snpPathMap[suppSNPs,:]
	sumReducedMap = np.sum(reducedMap,0)
	# exact matches (all SNPs the same)
	matchPaths = [i for i in range(len(sumReducedMap)) if sumReducedMap[i]==len(suppSNPs)]
	print 'exact matching paths:', matchPaths
	# paths with at least one matching SNP
	matchPaths = [i for i in range(len(sumReducedMap)) if sumReducedMap[i]>=2]
	print 'paths with > 2 matching SNP:', matchPaths

