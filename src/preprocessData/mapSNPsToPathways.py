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
Map SNPs to pathways and create PLINK genotype files for mapped SNPs only
Expand genotype design matrix to account for overlaps
Record all SNP->gene->pathway mappings for future results processing
Save all data required for analysis in single data file

Launch from directory containing this file, and with the following command:

python mapSNPsToPathways.py paramFile_name

where paramFile_name.py is the parameter file located in ../params/
'''

import cPickle, os, sys

import pathwaysGenes, geneSNPs, setPathAndParams
import pathwaysSNPs, getDesignMatrix
import initialiseData, myLogger


def mapSNPsToPathways(paramFile=None):
	
	print '\n\n-----------------------------------------'
#	print 'Mapping SNPs to Pathways'
	print '\n'
	
	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile)
	
	# mirror terminal output to log file
	logFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/' + params.PLINK_FILE_NAME.split('/')[-1] + "_pathwayMapping.log"
	logger = myLogger.Logger(file=logFn)
	stdout = sys.stdout; stderr = sys.stderr
	sys.stdout = logger; sys.stderr = logger
	

	# START PATHWAY MAPPING
	
	# generate pathway to gene map
	pathGeneMap = pathwaysGenes.PathGenes(params.PATHWAYS_DB_FILE)

	# generate gene to SNP map
	geneSNPmap = geneSNPs.GeneSNPs(params.GENE_LOCATIONS_FILE,params.PLINK_FILE_NAME,params.SNP_TO_GENE_MAPPING_RANGE)

	# map SNPs to pathways, apply filtering and and obtain list of mapped SNPs for genotype extraction
	pathSNPmap = pathwaysSNPs.PathSNPs(geneSNPmap.snpGeneMap,pathGeneMap.genePathMap,params) 

	# created filtered snp to gene and gene to SNP mappings 
	# that contains only snps and genes that are mapped to pathways
	geneSNPmap.createFilteredGeneSNPmap(geneSNPmap.snpGeneMap,pathSNPmap,pathGeneMap.genePathMap,params)

	# create PLINK files in required format using only mapped SNPs
	genotypes = getDesignMatrix.DesignMatrix(params,pathSNPmap)

	# create unexpanded and expanded pathway -> SNP maps using indices instead of names
	# pathway/SNP name to index mappings are created for future ref 
	pathSNPmap.createUnexpandedGroupMap(genotypes.snpRSidToSnpIdxMap)
	pathSNPmap.expand_Pathways()
	
	# create design matrix
	genotypes.expandDesignMatrix(pathSNPmap.orig2exp)
	

	# PICKLE EVERYTHING
	# 1. path-gene mapping
	pickleDumpFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/' + params.PLINK_FILE_NAME.split('/')[-1]
	picklePathwayGenesFh = open(pickleDumpFn + '_unfilteredPathwayGenes.pickle','w')
	cPickle.dump(pathGeneMap, picklePathwayGenesFh,-1)
	picklePathwayGenesFh.close()
	# 2. gene-snp mapping
	pickleGeneSNPFh = open(pickleDumpFn + '_geneSNPs.pickle','w')
	cPickle.dump(geneSNPmap, pickleGeneSNPFh,-1)
	pickleGeneSNPFh.close()
	# 3. pathway-snp mapping
	picklePathwayMapFh = open(pickleDumpFn + '_pathwayMapping.pickle','w')
	cPickle.dump(pathSNPmap, picklePathwayMapFh,-1)
	picklePathwayMapFh.close()
	# 4. genotype info
	pickleGenotypesFh = open(pickleDumpFn + '_genotypes.pickle','w')
	cPickle.dump(genotypes, pickleGenotypesFh,-1)
	pickleGenotypesFh.close()
	

	
	print 'all phenotype, genotype and pathway mapping data saved to file'	
	
	
	
	# create container for all data required in analysis and pickle
# 	data = initialiseData.Data(params,genotypes,phenotype,pathSNPmap)
	data = initialiseData.Data(params,genotypes,pathSNPmap)
	dataFn = pickleDumpFn + '_dataFile.pickle'
	dataFh = open(dataFn,'w')
	cPickle.dump(data, dataFh,-1) # '-1' flag most efficient protocol
	dataFh.close()
	print '\nfinal data structure required for analysis saved to:\n', dataFn
	
	print 'all info saved to log file: ', logFn
	sys.stdout = stdout # revert to unlogged output
	sys.stderr = stderr # revert to unlogged output
	

	


if __name__ == "__main__":
	
	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()
	
	mapSNPsToPathways(paramFile=paramFile)