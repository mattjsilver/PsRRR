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
#**************************************************************************************************
'''
set up python path and get program parameters
using specified parameter file
return parameters
'''

import sys, os

def setup(paramFile=None):
	
	sys.path.append('../params/');
	# can be called src/preprocessData/ or src/PsRRR/ 
	sys.path.append('../preprocessData/')
	sys.path.append('../PsRRR/')
	try:
		exec('import ' + paramFile)
	except:
		print 'unable to import param file:',paramFile,'...'
		print 'check param file name and location!\n\n'
		sys.exit()

	params = eval(paramFile).Params(printParams=True)
	
	# create required output directories if don't exist already
	OUTPUT_DIRS = [params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/',
						params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/',
						params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'pathwayMapping/']
	for directory in OUTPUT_DIRS:
		if not os.path.exists(directory):
			os.makedirs(directory)


	
	
	
	return params
