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
create w ~ group weights vector
'''

import numpy as np
import cPickle as pickle
from math import sqrt

def createWeights(params, stdWeights=False, groups=[], weightsFn=''):
	'''
	stdWeights = False - use weights in weightsFn
	stdWeights = True - use standard size weighting (sqrt(S_l))
	'''
	if stdWeights:

		# import group size data for weight calculation
		dataFh = open(params.DATA_FILE)
		data = pickle.load(dataFh)
		dataFh.close()
		
		if params.verbose:
			print 'computing standard size weighting'
		
		gWeights = map(sqrt,(map(len,data.groups.values())))
		gWeights = np.array(gWeights)
		gWeights.shape = (len(gWeights),1)
		stdWeightsFn = params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'adaptWeights/' + params.WEIGHTS_FILE + '_iteration0.pickle'
		stdWeightsFh = open(stdWeightsFn,'w')
		pickle.dump(gWeights, stdWeightsFh, -1)
		stdWeightsFh.close()
		
		if params.verbose:
			print 'std weights saved to:', stdWeightsFn
	
	else:
	
		weightsFn = weightsFn + '.pickle'
		print 'loading weights file:', weightsFn
		weightsFh = open(weightsFn,'r')
		gWeights = pickle.load(weightsFh)
		weightsFh.close()

	return gWeights	
