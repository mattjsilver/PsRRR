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
create object containing all genotype and SNP to pathway mapping data
required for PsRRR analysis
'''

class Data(object):

	def __init__(self,params,genotypes,pathSNPmap):

		# various params used for mapping
		self.params = params
		# genotype data
		self.X = genotypes.X
		self.eX = genotypes.eX
		self.N = genotypes.N
		self.P = genotypes.P
		self.eP = genotypes.eP
		self.snpRStoSnpIdx = genotypes.snpRSidToSnpIdxMap
				
		# pathway to SNP mapping data
		self.groups = pathSNPmap.groups
		self.eGroups = pathSNPmap.eGroups
		self.orig2exp = pathSNPmap.orig2exp
		self.duplicatePathways = pathSNPmap.duplicatePathways
		self.pathNamesToPathIdx = pathSNPmap.pathNameToPathIdxMap