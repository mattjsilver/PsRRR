#PATHWAYS SPARSE REDUCED-RANK REGRESSION
#
#Distributed under GNU general public licence (see COPYRIGHT.txt) 
#Copyright (C) 2012 Matt Silver
#
#WARRANTY DISCLAIMER
#Please note there is no warranty for the program, to the extent permitted by applicable law.
#Except when otherwise stated in writing the copyright holders and/or other parties provide the program
#"as is" without warranty of any kind, either expressed or implied, including, but not limited to, 
#the implied warranties of merchantability and fitness for a particular purpose. The entire risk
#as to the quality and performance of the program is with you. Should the program
#prove defective, you assume the cost of all necessary servicing, repair or correction.
#
#Please send comments or bugs to g.montana@imperial.ac.uk
#
#******************************************************************************************
'''
Pathways Sparse Reduced-rank Regression (PsRRR)

create pathway to SNP and SNP to pathway mappings using info already
gathered on pathway to gene and gene to SNP maps

Following filtering is applied:
1. remove pathways with less than params.MIN_SNPS_PER_PATHWAY SNPs
2. remove pathways with less than MIN_GENES_PER_PATHWAY genes ** NOT IMPLEMENTED **
3. remove identical pathways, i.e. pathways with the SAME SNPs; 
	for identical pathways we retain one 'representative' pathway
	and keep a list of matching pathways for ref	

Post-filtering, SNP variable expansion is performed, following the method of
		   Jacob, Laurent, Guillaume Obozinski, and Jean-philippe Vert
		   "Group Lasso with Overlap and Graph Lasso
		   Proc 26th Intl Conf on Machine Learning (2009).
	
Copyright (C) 2012 Matt Silver.  See docs/COPYRIGHT.txt
'''

from scipy import sparse
import numpy as np


class PathSNPs(object):
	'''
	All maps are post-filtering - i.e. non-mapped SNPs, genes and pathways are excluded
	
	PATHWAY TO SNP MAPPING AND VICE VERSA
	snpPathMap		- snpPathMap[SNP_RS_id] = [mappedPath1 mappedPath2...]
	pathSNPmap		- pathSNPmap[pathway] = [SNP_RS_id1 SNP_RS_id2...]
	
	P 				- number of SNPs that map to a pathway
	NsnpsNotMapped	- number of SNPs that DON'T map to pathway
	nMappedPathways - number of pathways that map to a SNP
	duplicatePathways - dictionary mapping retained pathways to removed identical pathways
	mappedSNPlist	- final list of SNPs mapped to pathways after filtering

	CONVERSION TO PATHWAY / SNP INDICES
	groups					- list of SNP indices mapping to each pathway	
	pathNameToPathIdxMap	- mappings of pathway names to pathway indices for future ref
	
	SNP	EXPANSION TO ACCOUNT FOR OVERLAPS
	eGroups  - as groups, but with unique SNP indices for SNPs in each group
	orig2exp - (P x eP) matrix mapping expanded SNPs back to original
	'''

	def __init__(self,snpGeneMap,genePathMap,params):
		'''
		@param snpGeneMap:
			snpGeneMap 	- snpGeneMap[snp_RS_id] = [gene1 gene2 ...] 
		@param genePathMap: 
			genePathMap - genePathMap[gene] = [path_1, path_2, ...]
		@params.MIN_SNPS_PER_PATHWAY: 
			minimum number of SNPs allowed per pathway			
		'''
		self.mapPathways(snpGeneMap,genePathMap) # perform pathway to SNP mapping and vice versa
		self.filterPathways(params) # apply filters to mapped pathways, e.g. filter out small ones
		self.printPathSNPmap(params)
		
		
	def mapPathways(self,snpGeneMap,genePathMap):
		'''
		map SNPs to pathways and pathways to SNPs
		'''
		
		print '\nMapping SNPs to pathways, and vice versa...'
		self.snpPathMap = {}
		self.NsnpsNotMapped = 0 # number of SNPs that don't map to a pathway
		self.pathSNPmap = {}

		for snp in snpGeneMap:
			if snp not in self.snpPathMap:
				self.snpPathMap[snp] = []
			else:
				print 'WARNING: DUPLICATE SNP', snp, 'IN SNP GENE MAP!'
			for gene in snpGeneMap[snp]: # genes mapped to this SNP
				if gene in genePathMap: # gene may not be in pathway database
					for pathway in genePathMap[gene]:
						if pathway not in self.snpPathMap[snp]: # pathway may already have been added to this SNP
							self.snpPathMap[snp].append(pathway)
						else: # create new key for this snp
							self.snpPathMap[snp] = [pathway]
						if pathway in self.pathSNPmap:
							# ****** ADDED 27/2/12 ************
							# avoid duplicated SNPs in pathway
							# arising from fact that SNPs may map to multiple genes in same pathway							
							if snp not in self.pathSNPmap[pathway]: 
								self.pathSNPmap[pathway].append(snp)
						else:
							self.pathSNPmap[pathway] = [snp]

			if self.snpPathMap[snp] == []: # no pathway mapped to this SNP
				del self.snpPathMap[snp]
				self.NsnpsNotMapped += 1
								
		self.P = len(self.snpPathMap)
		self.nMappedPathways = len(self.pathSNPmap)
		
		# print unfiltered mapping information
		print self.P, 'SNPs mapped to',self.nMappedPathways, 'pathways' 
		print self.NsnpsNotMapped, 'SNPs not mapped to any pathway'
		print 'SNPs map to max:', max(map(len,self.snpPathMap.values())), 'and min:', min(map(len,self.snpPathMap.values())), \
			'pathways'
		print 'pathways map to max:', max(map(len,self.pathSNPmap.values())), 'and min:', min(map(len,self.pathSNPmap.values())), \
			'SNPs'
		
	def filterPathways(self,params):
		'''
		APPLY FILTERING
		'''
		
		# A: pathway to SNP mapping
		print '\nApplying filtering...'
		# 1. remove pathways with < params.MIN_SNPS_PER_PATHWAY SNPs
		pathsToRemove = [path for path in self.pathSNPmap if len(self.pathSNPmap[path]) < params.MIN_SNPS_PER_PATHWAY]
		for delPath in pathsToRemove:
			del self.pathSNPmap[delPath]
		print len(pathsToRemove), 'pathways with less than', params.MIN_SNPS_PER_PATHWAY, 'SNPs removed'
		
		# 2. remove identical pathways, keeping list for future ref
		i = 0
		nIdPaths = 0
		self.duplicatePathways = {} # dict mapping retained pathways to deleted identical pathways
		while i < len(self.pathSNPmap):
			testPath = self.pathSNPmap.keys()[i] # path to test for matches
			for matchPath in self.pathSNPmap.keys()[i+1:]: # skip self
				if self.pathSNPmap[testPath] == self.pathSNPmap[matchPath]:
					if testPath in self.duplicatePathways: # duplicate for this pathway already found
						self.duplicatePathways[testPath].append(matchPath)
					else: # first duplicate found for this pathway
						self.duplicatePathways[testPath] = [matchPath]
					del self.pathSNPmap[matchPath]
					nIdPaths += 1
			i = i + 1

		print nIdPaths, 'duplicate pathways with identical SNPs removed; one representative pathway retained for all duplicate pathways'
		print 'list of representative and duplicate pathways will be saved to file'
		self.nMappedPathways = len(self.pathSNPmap)
		
		# B: SNP to pathway mapping
		self.snpPathMap = {} # revised SNP to pathway mapping after filtering
		for path in self.pathSNPmap:
			for snp in self.pathSNPmap[path]:
				if snp in self.snpPathMap:
					self.snpPathMap[snp].append(path)
				else:
					self.snpPathMap[snp] = [path]
					
		self.P = len(self.snpPathMap) # number of mapped SNPs after filtering
		
		print self.nMappedPathways, 'pathways mapped to', self.P, 'SNPs after filtering'
		print '\n* pathway to SNP mapping completed *'
		
		# extract list of mapped SNPs for genotype processing
		self.mappedSNPlist = self.snpPathMap.keys()
		
	
	def createUnexpandedGroupMap(self,snpRSidToIdxMap):
		'''
		change pathway names to integer indices
		keep track of these name -> idx mappings
		
		map pathways to SNPs using SNP indices corresponding to 
		columns in design matrix (snpRSidToIdxMap)
		'''
		
		print 'Creating pathway -> SNP maps with pathway and SNP names replaced by indices'
		print 'indices correspond to columns in design matrix...'

		# map pathway names to indices
		pathNames = self.pathSNPmap.keys()
		pathIndices = range(self.nMappedPathways)
		self.pathNameToPathIdxMap = dict(zip(pathNames,pathIndices))
		print '\npathway name -> pathway index map created'
				
		# convert pathSNPmap to indices
		self.groups = {} # dict containing index mappings
		for path in self.pathSNPmap:
			self.groups[self.pathNameToPathIdxMap[path]] = [snpRSidToIdxMap[snp] for snp in self.pathSNPmap[path]]		
		
		print 'unexpanded pathway -> SNP map with SNP indices corresponding to columns in design matrix created'

	
	def expand_Pathways(self):
		'''
		Expand groups to account for overlaps
		
		create pathway to SNP mapping without overlaps using the variable expansion
		process described in  
		   Jacob, Laurent, Guillaume Obozinski, and Jean-philippe Vert
		   "Group Lasso with Overlap and Graph Lasso
		   Proc 26th Intl Conf on Machine Learning (2009).
		
		also create mapping matrix to map back from latent to original variables
		'''

		self.eP = sum(map(len,self.groups.values())) # total SNPs in expanded var space
		print 'SNPs expanded from', self.P, 'to ', self.eP, 'SNPs'
		
		# eGroups is a list of pathway to SNP idx mapping in the expanded var space
		self.eGroups = []
		# orig2exp is a (p x eP) matrix that maps SNPs from the original var space
		# to the expanded space.  Each column of orig2exp represents a
		# (separately indexed) predictor, with a '1' in the row corresponding to
		# the 'original' predictor to which it refers.  It is thus row sparse
		self.orig2exp = sparse.lil_matrix((self.P,self.eP),dtype=np.uint8) 
	
		# perform expansion
		groupStartIdx = 0
		colIdx = 0
		for path in self.groups.keys():
			for SNP in self.groups[path]:
				self.orig2exp[SNP,colIdx] = 1
				colIdx += 1
			pathSize = len(self.groups[path])
			self.eGroups.append(range(groupStartIdx,groupStartIdx + pathSize,1))
			groupStartIdx = groupStartIdx + pathSize

		print 'expanded pathway -> SNP map created'
		print('%s matrix mapping SNP indices to indices in expanded variable space created'% str(self.orig2exp.shape))

	def printPathSNPmap(self, params):
		'''
		Print pathSNPmap to text file for ref
		'''
		pathSNPmapFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/' + params.PLINK_FILE_NAME.split('/')[-1] + "_pathways_SNP_mapping.txt"
		pathMapFh = open(pathSNPmapFn,'w')
		print 'saving pathway -> snp map to file:', pathSNPmapFn
		pathSNPmapKeys = self.pathSNPmap.keys()
		pathSNPmapKeys.sort()
		for path in pathSNPmapKeys:
			pathMapFh.write(path + '\t')
			for snp in self.pathSNPmap[path]:
				pathMapFh.write(' %s'%(snp))
			pathMapFh.write('\n')
		pathMapFh.close()