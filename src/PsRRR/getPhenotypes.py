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
create phenotype data instance
'''

import sys
import numpy as np
import cPickle

class Phenotypes(object):
	'''
	phenotypes first loaded from text file, then pickled
	pickled file used subsequently for speed
	'''
	
	def __init__(self,params):
				
		if params.PHENOTYPE_FILE.split('.')[-1] == 'pickle':
			# load from pickled file
			if params.verbose:
				print 'importing phenotype file:', params.PHENOTYPE_FILE, '...'
			phenotypesFh = open(params.PHENOTYPE_FILE)
			phenotypes = cPickle.load(phenotypesFh)
			phenotypesFh.close()
			self.Y = phenotypes.Y
			(self.N,self.Q) = phenotypes.Y.shape
			if params.verbose:
				print '(' + str(self.N) + ',' + str(self.Q) + ') phenotype matrix loaded'
			
		else: # phenotype file is not pickled - load from text file
		
			try:
				self.phenoFn = params.PHENOTYPE_FILE
				phenoFh = open(self.phenoFn)
				phenoFile_lines = phenoFh.readlines()
				phenoFh.close()
				
				# generate phenotype matrix from text file
				# omit column of subject IDs if necessary
				self.Y = np.genfromtxt(phenoFile_lines,dtype=np.float32)
				if params.PHENOTYPE_FIRST_COL_IDs:
					self.Y = self.Y[:,1:]
					
				if len(self.Y.shape) == 1: # univariate phenotype, reshape to (N,1) for compatibility
					self.Y.shape = (self.Y.shape[0],1)
				(self.N,self.Q) = self.Y.shape
				print '(' + str(self.N) + ',' + str(self.Q) + ') phenotype matrix loaded'
				
				# pickle for future speed
				fh = open(params.PHENOTYPE_FILE + '.pickle','w')
				cPickle.dump(self,fh,-1)
				fh.close()
				
			except:
				print 'failed to load phenotype file:', params.PHENOTYPE_FILE
				print 'exiting..'
				sys.exit()