import sys, os, datetime

class Params:
	def __init__(self,printParams=True,simulate=False):

		# print params?
		if printParams == False:
			stdout = sys.stdout
			sys.stdout = open('/dev/null','w')
			
		print '[params] loading program parameters'


		#***************************************************************************************
		# FILE LOCATIONS
		self.PsRRR_ROOT_DIR			= os.path.abspath('../..') + '/' # program root directory, containing src, data, docs etc.  

		# SNP to pathway mapping
		self.MAPPING_OUTPUT_DIR = '/sample_output/' # all files created during the pathway mapping process will be placed here
		self.PATHWAYS_DB_FILE = self.PsRRR_ROOT_DIR + 'data/genePathways.map' # pathway mapping file
		self.GENE_LOCATIONS_FILE = self.PsRRR_ROOT_DIR + 'data/genePos_GRCh37p3.csv' # gene locations - must match genome build of genotypes
		self.PLINK_FILE_NAME = self.PsRRR_ROOT_DIR + 'data/sample' 	# plink filename without bed/bim/fam suffix
		self.PLINK_PROG_FOLDER	= '/home/silverm/plink/' # directory containing PLINK executable - include trailing '/'
		
		# PsRRR
		self.PsRRR_OUTPUT_DIR = 'sample_output/' # all files created during runtime will be placed here - include trailing '/'
		self.DATA_FILE 				= self.PsRRR_ROOT_DIR + 'sample_output/pathwayMapping/sample_dataFile.pickle'
		# phenotype file is white space delimited text file with N x Q array of values, where N = sample size, and Q = number of phenotypes
		# (Q=1 for univariate phenotype).  First column can be list of IDs if desired (set self.PHENOTYPE_FIRST_COL_IDs = True)
		self.PHENOTYPE_FILE 		= self.PsRRR_ROOT_DIR + 'data/phenotypes.txt' # multivariate pheno for GL
		# self.PHENOTYPE_FILE 		= self.PsRRR_ROOT_DIR + 'data/univariate_phenotypes.txt' # univariate pheno for SGL
		self.PHENOTYPE_FIRST_COL_IDs = True # 1st column of phenotype file is subject IDs if 'True'
		self.STATUS_FILE			= 'status.txt'  # for subject disease status when using balanced subsampling (see getData.createSubsamp())
		self.WEIGHTS_FILE 			= 'sample_weights'	# weights file rootname
		self.weights_iteration		= 0	# adaptive weights iteration to use for analysis
										# '0' specifies standard group size weighting ( w_l = sqrt(S_l) )
		
		print '\nSOURCE DATA FILES:'
		for printItem in ['PATHWAYS_DB_FILE','GENE_LOCATIONS_FILE','PLINK_FILE_NAME','PLINK_PROG_FOLDER','DATA_FILE',
						'PHENOTYPE_FILE','PHENOTYPE_FIRST_COL_IDs','STATUS_FILE',	'WEIGHTS_FILE','weights_iteration']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************
		
		#***************************************************************************************
		# SNP TO PATHWAY MAPPING
		self.SNP_TO_GENE_MAPPING_RANGE = 	10000 	# (bp) Any SNPs +/- SNP_TO_GENE_MAPPING_RANGE bp will be mapped to corresponding gene
		self.MIN_SNPS_PER_PATHWAY = 		10 		# pathways with < MIN_SNPS_PER_PATHWAY mapped SNPs will be excluded 
		self.MIN_GENES_PER_PATHWAY = 		1		# pathways with < MIN_GENES_PER_PATHWAY mapped genes will be excluded

		print '\nSNP TO PATHWAY MAPPING PARAMS:'
		for printItem in ['SNP_TO_GENE_MAPPING_RANGE', 'MIN_SNPS_PER_PATHWAY', 'MIN_GENES_PER_PATHWAY']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		
		#***************************************************************************************
		
		#***************************************************************************************
		# REGULARISATION
		self.bPenalty				= 'GL' # sparsity constraint on 'b' (genotype) coeff vector
											# 'GL' = group lasso (standard setting)
											# 'SGL'= sparse group lasso
											# 'lasso' = lasso
											# 'lasso_post_PsRRR' = lasso on groups selected in previous PsRRR analysis
		self.aPenalty				= None 	# sparsity constraint on 'a' (phenotype) coeff vector 
											# None or 'lasso' only
		self.Lambda 				= 0.85  # lambda value to use in estimation 
											# self.lambda = lambda/lambda_max where:
											# lambda - is actual value used in estimation
											# lambda_max is smallest lambda where no groups are selected
		self.alpha					= 0.85   # SGL only; distribution of l1 / l2 penalties (0 < alpha < 1)
											# alpha = 0 is pure GL; alpha = 1 is pure lasso 

		print '\nGROUP LASSO PARAMS:'
		for printItem in ['bPenalty', 'aPenalty', 'Lambda','alpha']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************			
		
		#***************************************************************************************			
		# PATHWAYS SPARSE REDUCED RANK REGRESSION
		self.PsRRRconvergenceThresh	= 1e-3 # threshold for PsRRR convergence, epsilon where 1-|b'b| < epsilon
		self.max_PsRRR_iterations	= 20   # maximum number of allowed iterations in PsRRR descent
		print '\nPsRRR PARAMS:'
		for printItem in ['PsRRRconvergenceThresh','max_PsRRR_iterations']:
			print "%s = '%s'" % (printItem, getattr(self, printItem)) 			
		
		# Number of N/2 subsamples for ranking pathways, SNPs etc
		self.B						= 400

		print '\nBootstrap PARAMS:'
		for printItem in ['B']:
			print "%s = '%s'" % (printItem, getattr(self, printItem)) 			
		#***************************************************************************************
		
		#***************************************************************************************				
		# WEIGHT TUNING
		self.R 						= [5,5,5,5,5]	# list of number of model fits at each 
											# weight tuning iteration.		
											# weights are adapted a total of n = len(R) times
											# with R(i) increasing for i = 1, ..., n
											# we begin with 'course' tuning, over R(1) iterations, where R(1) ~ 2 x nPathways
											# since the weights are expected to be heavily biased at the beginning, meaning
											# that certain pathways will be chosen many more than 1/R(1) times
		self.firstIteration			= 0 	# where to start iterative adaptation process:
											# '0' indicates start with standard size weighting
											# '1' indicates start with params.WEIGHTS_FILE_iteration1 etc
		self.eta					= [0.95,0.95,0.99,0.99,0.99]
											# maximum proportionate weight decrease to be used at each iteration
											# thus eta[iteration] is the value used at [iteration] to calculate the
											# weights to be used at [iteration+1]
											# this should decrease slowly as weights become better-tuned
		if len(self.R) != len(self.eta):
			print 'number of eta values different from number of iterations, please amend params file\nexiting...'
			sys.exit()

									
		print '\nADAPTIVE WEIGHTS PARAMS:'
		for printItem in ['R','eta', 'firstIteration']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************									
									
		#***************************************************************************************		
		# SIMULATION PARAMETERS
		if simulate:
			self.N						= 400	# sample size
			self.P						= 1010	# number of SNPs
			self.nGroups				= 50	# number of pathways; each pathway contains P/nGroups
			self.groupSize				= 30	# pathway size = P/nGroups for non-overlapping pathways ONLY
			self.overlap				= True  # allow groups to overlap?
			self.overlapSize			= 10	# n overlapping SNPs between adjacent groups
			self.constMAF				= True 	# simulate genotypes with constant MAF?
			self.Q						= 1		# phenotype dimensionality (set to 1 for univariate phenotype
			self.nCausalPheno			= 1		# number of affected phenotypes for each effect (set to 1 for univariate phenotype)
			self.nSNPs					= 10	# number of causal SNPs
			self.effSize				= 0.5	# constant effect size per SNP for each rank
			# sanity checks:
			# check that number of SNPs is divisible by number of pathways
			if not self.overlap and self.P % self.nGroups:
				if printParams == False:
					sys.stdout = stdout
				print 'total number of SNPs:', self.P, 'not divisible by number of pathways:', self.nGroups
				print 'ammend params file'
				print 'exiting...'
				sys.exit()
			# check that number of causal phenotypes < number of phenotypes
			if self.nCausalPheno > self.Q:
				print 'number of causal phenotypes > phenotype dimensionality'
				print 'ammend params file'
				print 'exiting now'
				sys.exit()

			print '\nSIMULATED RESPONSE PARAMS:'
			for printItem in ['N','P','nGroups','Q','nCausalPheno','nSNPs','effSize']:
				print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************		
		
		#***************************************************************************************
		# OPERATING MODE (PARALLELISATION)
		self.mode					= 'pp' # parallel 	processing mode:
												# 'std' 	= no parallel processing
												# 'pp' 		= use parallel python across cluster

		print '\nPARALLEL PROCESSING MODE:'
		for printItem in ['mode']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************

		#***************************************************************************************		
		# PARALLEL PYTHON
		# tuple of remote server names:
		self.ppServers 				= ("Montana-02","Montana-03","Montana-04","Montana-05","Montana-06","Montana-07","Montana-08","Montana-09")
		self.ppPort					= 60000 	# port to communicate between client (master) and remote server nodes
		self.workersPerServer		= 5	# number of CPUs to use on each remote server
		self.nLocalCpus				= 0		# number of *LOCAL* cpus to use as workers (used for testing PP on single multi-cpu workstation)
		
		print '\nPARALLEL PYTHON PARAMS:'
		for printItem in ['ppServers','ppPort','workersPerServer','nLocalCpus']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************

		#***************************************************************************************
		# MISC
		self.balanced				= False # ensure subsamples are balanced w.r.t. disease status
		self.verbose 				= True  # verbose output
		self.time 					= datetime.datetime.now().strftime("%d_%b_%Y_%H_%M") # for time stamping logfiles
		self.generateNull			= False # permute subjects in Y only when True
		
		print '\nMISC PARAMS:'
		for printItem in ['verbose','time']:
			print "%s = '%s'" % (printItem, getattr(self, printItem))
		#***************************************************************************************
			

		if printParams == False:
			sys.stdout = stdout