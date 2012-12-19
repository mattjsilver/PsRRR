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
generate unexpanded design matrix, X, and equivalent expanded 
design matrix, eX, containing only SNPs mapped to pathways
'''

import os
import numpy as np
import sys

class DesignMatrix(object):
	'''
	P					- number of mapped SNPs
	N					- sample size
	X					- (N x P) matrix of SNP minor allele counts
	snpRSidToSnpIdxMap	- 1 to 1 mapping of SNP RSid to column in design matrix
	eP					- number of SNPs in expanded design matrix
	eX					- (N x eP) matrix of SNP minor allele counts in expanded var space
	
	'''
		
	def __init__(self,params,pathSNPmap):
		'''
		@param: params
			source data filenames etc			 
		@param: pathSNPmap
			pathway to SNP mapping object
		'''

		# create updated PLINK files using only mapped SNPs
		pedFn = self.processPLINKfiles(params,pathSNPmap)
		# generate design matrix from this file
		self.generateDesignMatrix(pedFn)
		
		print '\n* genotype file processing completed *'		
		
		
	def processPLINKfiles(self,params,pathSNPmap):
		
		print '\n\n-----------------------------------------'
		print 'Updating PLINK files and creating design matrix...'	
		
		# create new bed/bim/fam files containing SNPs mapped to pathways only
		# First write list of mapped SNPs to text file
		newPLINKfn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/' + params.PLINK_FILE_NAME.split('/')[-1]
		mappedSNPs_fn = newPLINKfn + '_mappedSNPs.txt'
		mappedSNPs_fh = open(mappedSNPs_fn,"w")
		
		
		for snp in pathSNPmap.mappedSNPlist:
			mappedSNPs_fh.write(snp+"\n")
		mappedSNPs_fh.close()
	
		newPLINKfn = newPLINKfn + '_mappedSNPsOnly'
		print '\nUpdating PLINK files to include only mapped SNPs'
		print '(see ' + newPLINKfn + '.log)'
	
		cmd = params.PLINK_PROG_FOLDER + 'plink --bfile ' + params.PLINK_FILE_NAME + ' --extract ' + mappedSNPs_fn + ' --make-bed --out ' + newPLINKfn + ' --noweb'
		sysError = os.system(cmd + ' > /dev/null') # '> /dev/null' suppresses PLINK output
		if sysError:
			print '!!PLINK command failed!!'
			sys.exit()
		
		# recode bed file to allele dosages, and convert to .ped file
		# note --recodeA recodes to dosage for *MINOR* allele
		print 'recoding PLINK file to minor allele dosage'
		cmd = params.PLINK_PROG_FOLDER + 'plink --bfile ' + newPLINKfn + ' --recodeA --out ' + newPLINKfn + '_recoded'  + ' --noweb'
		sysError = os.system(cmd + ' > /dev/null') # '> /dev/null' suppresses PLINK output
		if sysError:
			print '!!PLINK command failed!!'
			sys.exit()
		
		newPedfn = newPLINKfn + '_recoded.raw'
		print 'genotype file recoded and saved as ', newPedfn 	

		return newPedfn


	def generateDesignMatrix(self,pedFn):
		print 'Generating genotype design matrix for mapped SNPs only...'

		self.pedFn = pedFn
		pedFh = open(self.pedFn)
		pedFile_lines = pedFh.readlines()
		pedFh.close()
		
		# list of SNPs corresponding to columns in the design matrix (taken from ped file header)
		plinkOrderedSNPs = pedFile_lines[0].split()[6:] # header, excluding fam_id etc
		plinkOrderedSNPs = [snp.split('_')[0] for snp in plinkOrderedSNPs] # rs_id only (PLINK adds minor allele to this)
		
		# map SNP rsIDs to integer indices
		snpIndices = range(len(plinkOrderedSNPs))
		self.snpRSidToSnpIdxMap = dict(zip(plinkOrderedSNPs,snpIndices))

		
		# create design matrix - skip header and 1st 6 cols) - use int8 to save alot of space
		self.X = np.genfromtxt(pedFile_lines,skip_header=1,missing_values='NA', filling_values='0', \
										usecols=range(6,len(plinkOrderedSNPs)+6),dtype=np.uint8)
		(self.N,self.P) = self.X.shape

		print self.X.shape, 'design matrix generated for', self.N, 'subjects and', self.P, 'SNPs'
				
		
	def expandDesignMatrix(self,orig2exp):
		# eX is (n x eP) expanded design matrix with columns re-ordered according to 
		# their ordering in groups, AND with duplicated columns for overlapping predictors 
		# i.e. those that occur in multiple groups
		
		print '\n\n-----------------------------------------'
		print 'Creating expanded design matrix...\n'

		self.eX = self.X*orig2exp
		self.eP = self.eX.shape[1]


		print np.shape(self.eX), 'Expanded design matrix created'		