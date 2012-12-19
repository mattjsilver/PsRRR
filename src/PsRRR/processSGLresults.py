'''
Created on Jun 24, 2012
process pathway results with new text file format
@author: mattsilver
'''

import operator,sys,cPickle
import numpy as np
import setPathAndParams


def processResults(paramFile=None):

	ranksToView = 30 # number of top ranks, by selection frequency, to view
	pathResultsFn = 'SGL_22_Aug_2012_16_36_pathway_results.txt' # name of pathwa results file in results/ directory
	snpResultsFn = 'SGL_22_Aug_2012_16_36_snp_results.txt' # name of SNP results file in results/ directory 
	
	# set python path, get params, and create output directories if necessary
	params = setPathAndParams.setup(paramFile=paramFile,printParams=False)

#	# load pathway results
#	pathResFh = open(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + pathResultsFn)
#	pathRes = resFh.readlines()
#	pathResFh.close()
#	
#	# load SNP results
#	snpResFh = open(params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + snpResultsFn)
#	snpRes = snpResFh.readlines()
#	snpResFh.close()


	printPathwayResults(pathResFile=params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + pathResultsFn, 
								ranksToView=ranksToView, params=params)

 	printSNPresults(snpResFile=params.PsRRR_ROOT_DIR + params.PsRRR_OUTPUT_DIR + 'results/' + snpResultsFn, 
											ranksToView=ranksToView,params=params)
	

	
	

def printSNPresults(snpResFile=None,ranksToView=None,params=None):

	print '\nProcessing SNP results...\n'

	# load results
	resFh = open(snpResFile)
	subsamps = resFh.readlines()
	resFh.close()
	
	print 'snp selection frequencies over', len(subsamps), 'subsamples recorded\n'

	# load snp idx to snp RS id map
	fh = open(params.DATA_FILE)
	data = cPickle.load(fh)
	fh.close()
	rsToIdx = data.snpRStoSnpIdx
	snpIdxToRS = dict (zip(rsToIdx.values(),rsToIdx.keys()))
	
	# load snp RS to gene map
	fh = open(params.DATA_FILE.split('dataFile.pickle')[0] + 'geneSNPs.pickle')
	geneSNPs = cPickle.load(fh)
	fh.close()
	snpGeneMap = geneSNPs.filteredSnpGeneMap


	# unlike pathways, we only record selected SNPs
	snpFreqs = {} # tally of SNP selection frequencies
	geneFreqs = {} # tally of gene selection frequencies
	nSelSNPs = 0 # cumulative tally of number of selected SNPs	
	nSelGenes = 0 # cumulative tally of number of selected genes
	nSubSamps = len(subsamps) # number of subsamps run
	
	for subsamp in subsamps:
#		res = eval(subsamp[:-2]) # list of SNPs selected at this subsamp
		res = subsamp.split() # list of SNPs selected at this subsamp
		nSelSNPs += len(res) # keep tally of total number of selected SNPs

		mappedGenesThisSubsamp = []
		
		# update SNP selfreqs
		for snp in res:
			if snp in snpFreqs:
				snpFreqs[snp] += 1/float(nSubSamps)
			else:
				snpFreqs[snp] = 1/float(nSubSamps)

			# update genes selected at this subsamp
			snpRS = snpIdxToRS[int(snp)]
			mappedGenesThisSubsamp.append(snpGeneMap[snpRS][0])
			
		# remove dupicate genes (same gene maps to multiple selected SNPs)
		mappedGenesThisSubsamp=list(set(mappedGenesThisSubsamp))
		nSelGenes += len(mappedGenesThisSubsamp)

		# updated gene selFreqs
		for gene in mappedGenesThisSubsamp:
			if gene in geneFreqs:
				geneFreqs[gene] += 1/float(nSubSamps)
			else:
				geneFreqs[gene] = 1/float(nSubSamps)

	
	sortedSNPfreqs = sorted(snpFreqs.iteritems(), key=operator.itemgetter(1))
	sortedSNPfreqs.reverse()
	sortedGeneFreqs = sorted(geneFreqs.iteritems(), key=operator.itemgetter(1))
	sortedGeneFreqs.reverse()
	
	print 'mean SNPs selected per subsample:', nSelSNPs/float(nSubSamps)
	print 'mean genes selected per subsample:', nSelGenes/float(nSubSamps)	
		
	
	# print SNP results
	print '============'
	print 'SNP RANKINGS\n'
	
	rank = 1
	print 'rank\tpathway\t\tselFreq\tmapped gene(s)'
	for snp in sortedSNPfreqs[0:ranksToView]:
		snpRS = snpIdxToRS[int(snp[0])]
		print rank, '\t', snpRS, '\t', snp[1], '\t',
		if len(snpGeneMap[snpRS])>1: # more than one gene mapped to this SNP
			for gene in snpGeneMap[snpRS]:
				print gene,'\t',
			print
		else:
			print snpGeneMap[snpRS][0]
		rank += 1


	# print gene results
	print '\n\n============'
	print 'GENE RANKINGS\n'
	rank = 1
	print 'rank\tgene\tselfreq'
	for gene in sortedGeneFreqs[0:ranksToView]:
		print rank, '\t', gene[0], '\t', gene[1]
		rank += 1
		
	return (snpFreqs,geneFreqs)



def printPathwayResults(pathResFile=None,ranksToView=20,params=None):

	print '\nProcessing pathway results...\n'
	
	# load pathway index info
 	fh = open(params.DATA_FILE)
 	data = cPickle.load(fh)
 	fh.close()
 	
 	pathNamesToIdx = data.pathNamesToPathIdx
	pathIdx = dict (zip(pathNamesToIdx.values(),pathNamesToIdx.keys()))
	nPaths = len(pathIdx) # number of pathways

	# load results
	resFh = open(pathResFile)
	subsamps = resFh.readlines()
	resFh.close()
	nSubsamps = len(subsamps)

	print 'pathway selection frequencies over', len(subsamps), 'subsamples recorded\n'
	
	pathFreqs = dict.fromkeys(range(nPaths),0) # tally of pathway selection frequencies
	
	nSelPaths = 0 
	for subsamp in subsamps:
		res = subsamp.split() # list of pathways selected at this subsamp
		nSelPaths += len(res)
		for path in res:
			pathFreqs[int(path)] += 1/float(nSubsamps)

	# calculate mean n selected pathways across all subsamps
	print 'mean sel paths =', nSelPaths/float(nSubsamps), '\n'
	
	sortedPathFreqs = sorted(pathFreqs.iteritems(), key=operator.itemgetter(1))
	sortedPathFreqs.reverse()

	# print some results
	print 'rank\tpathway\t\tselFreq'
	rank = 1
	for path in sortedPathFreqs[0:ranksToView]:
		print rank, '\t', pathIdx[int(path[0])], '\t', path[1]
		rank += 1
		
	return pathFreqs



if __name__ == '__main__':
	
	try:
		paramFile = sys.argv[1] # parameter file should be used as arg
	except:
		print 'please supply parameter file as command line argument'
		sys.exit()
	
	processResults(paramFile=paramFile)