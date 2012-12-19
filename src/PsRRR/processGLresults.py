'''
Created on Jun 24, 2012
process pathway results with new text file format
@author: mattsilver
'''

import sys, operator, cPickle
import setPathAndParams

def processResults(paramFile=None):
	
	pathsToView = 30 # number of top ranks, by selection frequency, to view
	resultsFn = 'GL_22_Aug_2012_16_36_pathway_results.txt' # name of results file in results/ directory 
	
	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile,printParams=False)

	resFh = open(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + resultsFn)
	empRes = resFh.readlines()
	resFh.close()
	
	# load pathway index to pathway name mappings
	dataFh = open(params.DATA_FILE)
	data = cPickle.load(dataFh)
	dataFh.close()
	
	namesToPathIdx = data.pathNamesToPathIdx
	pathIdxToNames = dict (zip(namesToPathIdx.values(),namesToPathIdx.keys()))
	nPaths = len(pathIdxToNames)

	empFreqs = getSelFreqs(empRes,nPaths)
	
	print 'pathway selection frequencies over', len(empRes), 'subsamples recorded\n'
		
	# print some results
	print 'TOP PATHWAY RANKS (top', pathsToView, 'only)\n' 
	print 'rank\tpathway\t\tselFreq'
	for idx in range(min(pathsToView,len(empFreqs))): # in case fewer pathways are ranked than required pathsToView 
		path = empFreqs[idx]
		print str(idx+1) + '\t' + pathIdxToNames[int(path[0])] + '\t' + str(path[1])
	
	

	
def getSelFreqs(results,nPaths):
	
	pathFreqs = dict.fromkeys(range(nPaths),0) # tally of pathway selection frequencies # tally of pathway selection frequencies
	nSubsamps = len(results)
	
	for subsamp in results:
		res = subsamp.split() # list of pathways selected at this subsamp
		for path in res:
			if path in pathFreqs:
				pathFreqs[path] += 1./nSubsamps
			else:
				pathFreqs[path] = 1./nSubsamps
	
	sortedPathFreqs = sorted(pathFreqs.iteritems(), key=operator.itemgetter(1))
	sortedPathFreqs.reverse()
	
	
	return (sortedPathFreqs)

if __name__ == '__main__':
	
	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()

	
	processResults(paramFile=paramFile)