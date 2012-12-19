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
create pathway to gene and gene to pathway map from pathways db file
'''

class PathGenes(object):
	'''
	pathGeneMap - pathGeneMap[pathway] = [gene_1, gene_2, ...]
	genePathMap - genePathMap[gene] = [path_1, path_2, ...]
	nPaths		- number of pathways in database file
	nGenes		- number of unique genes in pathways database
	'''

	def __init__(self,pathDBfileName):
		'''
		@param pathDBfileName: 
			white space delimited text file with following format:
				pathwayName   gene_1   gene_2...
		'''

		print '\nMapping pathways to genes and genes to pathways...'
		print 'using', pathDBfileName
		self.pathDBfname = pathDBfileName
		self.generate_pathGeneMap()
		self.generate_genePathMap()
		print self.nGenes, 'unique genes mapped to ', self.nPaths, 'pathways'
		
		
	def generate_pathGeneMap(self):
		'''
		self.pathGeneMap: dict of genes mapped to each pathway
		'''
		pathDBfh = open(self.pathDBfname,'r')
		lines = pathDBfh.readlines()
		pathDBfh.close()

		self.pathGeneMap = {}
		for path in lines:
			self.pathGeneMap[path.split()[0]] = path.upper().split()[1:] # extract gene symbols for each pathway
		
		self.nPaths = len(self.pathGeneMap) # number of pathGeneMap
		
		
	def generate_genePathMap(self):
		'''
		self.gene: dict of pathways mapped to each gene
		'''
		self.genePathMap = {}
		for path in self.pathGeneMap:
			for gene in self.pathGeneMap[path]:
				if gene in self.genePathMap:
					if path not in self.genePathMap[gene]: # many genes duplicated in pathways files!
						self.genePathMap[gene].append(path)
				else:
					self.genePathMap[gene] = [path]
		self.nGenes = len(self.genePathMap) # number of pathGeneMap