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
Run PsRRR.  Launch from directory containing this file, and with the following command:

python runPsRRR.py paramFile_name

where paramFile_name.py is the parameter file located in ../params/
'''

import sys, setPathAndParams, PsRRR


def runPsRRR(paramFile=None):
	
	print
	print '###########################################################'
	print '####### PATHWAYS SPARSE REDUCED RANK REGRESSION	##########'
	print '###########################################################'
	
	
	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile)
	
	
	print '\n------------------------------------------'
	print params.bPenalty,'PENALTY'
	print '------------------------------------------\n\n'
	
	params.adaptWeights = False

	if params.generateNull: # do pathway selection under the null, with correspondence between X and Y subjects broken
		print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
		print '!!!!!!!		WARNING		params.generateNull = True....			!!!!!!!!!!'
		print '!!!!!!!    	PERMUTING PHENOTYPE SUBJECTS TO GENERATE NULL DATA      !!!!!!'
		print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	
			
	results = PsRRR.setup_PsRRR(params=params, data=None, bPenalty=params.bPenalty, aPenalty=params.aPenalty, subsample=True)

	# save results
	fn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + params.bPenalty + '_' + params.time + '_pathway_results.txt'
	fh = open(fn,'w')

	if params.mode == 'pp':
		
		# parse threads and subsamps within each thread
		for thread in results:
			for subsampRes in thread()['selGroups']:
				for path in subsampRes:
					fh.write(str(path)+' ')
				fh.write('\n')
		fh.close()
		print 'pathway results saved to', fn


		if params.bPenalty == 'SGL':
			fn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + params.bPenalty + '_' + params.time + '_snp_results.txt'
			fh = open(fn,'w')
			for thread in results:
				for subsampRes in thread()['selSNPs']:
					for snp in subsampRes:
						fh.write(str(snp)+' ')
					fh.write('\n')
			fh.close()
			print 'SNP selection results saved to', fn
			
			
	else:
		for subsampRes in results:
			for path in subsampRes['selGroups'][0]:
				fh.write(str(path)+' ')
			fh.write('\n')
		fh.close()
		print 'pathway results saved to', fn

		if params.bPenalty == 'SGL':
			fn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + params.bPenalty + '_' + params.time + '_snp_results.txt'
			fh = open(fn,'w')
	
			for subsampRes in results:
				for snp in subsampRes['selSNPs'][0]:
					fh.write(str(snp)+' ')
				fh.write('\n')
			fh.close()
			print 'SNP selection results saved to', fn
	
	
	
	
	
	


if __name__ == "__main__":

	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()

	runPsRRR(paramFile)
