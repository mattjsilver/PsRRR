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
Create gene to SNP and SNP to gene maps using info on gene and SNP locations
'''

import re

class GeneSNPs(object):
	'''
	geneLocs	- geneLocs[gene HGNC symbol] = [gene_start(bp)  gene_end(bp)  chr]
	nGeneLocs	- number of gene locations extracted
	SNPlocs	 	- SNPlocs[SNP_RS_id] = [SNP_location(bp)  chr]
	nSNPlocs	- number of SNP locations extracted
	geneSNPmap 	- geneSNPmap[gene] = [SNP1_RS_id SNP2_RS_id ...]
	snpGeneMap 	- snpGeneMap[snp_RS_id] = [gene1 gene2 ...]
	SNPsMappedToGenes 		- list of SNPs that map to genes within mappingRange
	nSNPsMappedToGenes 		- number of SNPs that map to genes within mappingRange
	nSNPsNotMappedToGenes 	- number of SNPs that don't map to a gene within mappingRange
	nGenesMappedToSNPs		- number of genes that map to SNPs within mappingRange
	'''
	
	def __init__(self,geneLocsFileName,SNPlocsFileName,mappingRange):
		'''
		@param geneLocsfileName: 
			fn of comma separated text file with following format:
			geneStart, geneEnd, chromsome, HGNC symbol
		@param SNPlocsfileName: 
			fn of PLINK bim/fam/bed file. *.bim is tab delimited file with following format:
			chromosome	SNP_RSid	distance (ignored)	location(bp)	major/minor alleles
		@param mappingRange: 
			maximum SNP to gene mapping range (bp) within which SNP can be mapped to 
			corresponding gene
		'''
		
		self.geneLocsFname = geneLocsFileName
		self.SNPLocsFname = SNPlocsFileName+'.bim' # PLINK .bim file contains SNP locations
		self.SNPtoGeneMappingRange = mappingRange
		print '\nMapping genotyped SNPs to genes in pathways database, and vice versa...'
		self.generate_gene_location_map()
		print self.nGeneLocs, 'gene locations extracted from', geneLocsFileName
		self.generate_SNP_location_map()
		print self.nSNPlocs, 'SNP locations extracted from', self.SNPLocsFname
		print 'mapping SNPs to genes (this may take some time)...'
		self.map_snps_to_genes()
		print self.nSNPsMappedToGenes, 'SNPs mapped to', self.nGenesMappedToSNPs, 'genes in pathways database'
		print self.nSNPsNotMappedToGenes, 'not mapped to any gene'

		
#		
	def generate_gene_location_map(self):
		geneLocsFh = open(self.geneLocsFname,'rU')
		self.lines = geneLocsFh.readlines()[1:] # skip header
		geneLocsFh.close()

		self.geneLocs = {}
		for line in self.lines:
			line = line.split(',') # comma separated values
			chrRegex = re.compile('^[1-9]$|^1[0-9]$|^2[0-2]$|^X$') # matches any chromosome id incl 'X'
			hgnc = line[3].rstrip() # hgnc gene symbol
			if hgnc != '': # only interested in lines with gene symbol - many lines don't
				if chrRegex.match(line[2]): # check valid chromosome id
					if line[2] == 'X': # recode X to '23' (PLINK coding)
						line[2] = '23'
					self.geneLocs[hgnc] = map(int,line[0:3])	# create dict entry for this gene

		self.nGeneLocs = len(self.geneLocs) # number of gene locations extracted		

	def generate_SNP_location_map(self):
		SNPlocsFh = open(self.SNPLocsFname,'rU')
		lines = SNPlocsFh.readlines()
		SNPlocsFh.close()

		self.SNPlocs = {}
		for line in lines:
			line = line.split()
			self.SNPlocs[line[1]] = [int(line[3]), int(line[0])]

		self.nSNPlocs = len(self.SNPlocs) # number of snp locations extracted		

	def map_snps_to_genes(self):
		self.snpGeneMap = {}
		self.geneSNPmap = {}
		self.SNPsMappedToGenes = [] # list of SNPs that do map to a gene
		self.nSNPsNotMappedToGenes = 0
		self.nSNPsMappedToGenes = 0 # number of unmapped SNPs
		self.nGenesMappedToSNPs = 0;
		justPrinted = 0 # avoid double printing screen
		# loop through SNPs
		
		for snp in self.SNPlocs:
			snpMapped = False;
			snpChr = int(self.SNPlocs[snp][1])
			snpLoc = int(self.SNPlocs[snp][0])
			# loop through genes
			for gene in self.geneLocs:
				geneChr = self.geneLocs[gene][2]
				if snpChr == geneChr:
					if (self.geneLocs[gene][0] <= snpLoc + self.SNPtoGeneMappingRange) and (self.geneLocs[gene][1] >= snpLoc - self.SNPtoGeneMappingRange):
# 						print 'snp mapped'
						snpMapped = True;
						if snp in self.snpGeneMap: # map genes to SNPs
							self.snpGeneMap[snp].append(gene) 
						else:
							self.snpGeneMap[snp] = [gene]
							self.SNPsMappedToGenes.append(snp)
							self.nSNPsMappedToGenes += 1
						if gene in self.geneSNPmap: # map SNPs to genes
							self.geneSNPmap[gene].append(snp) 
						else:
							self.geneSNPmap[gene] = [snp]
							self.nGenesMappedToSNPs += 1	
						
			if not(snpMapped):
				self.nSNPsNotMappedToGenes += 1;
			# print some output to screen as this is slow!
			nSNPsProgress = len(self.SNPsMappedToGenes)
			if not(nSNPsProgress%5000) and nSNPsProgress != justPrinted:
				print len(self.SNPsMappedToGenes),'SNPs mapped'
				justPrinted = nSNPsProgress


	def createFilteredGeneSNPmap(self,snpGeneMap,pathSNPmap,genePathMap,params):
		'''original snpGeneMap is unfiltered - that is it 
		includes snps that don't map to genes, and also
		genes that don't map to pathways
		Method produces a filtered snpGeneMap and geneSNPmap
		that contains only snps and genes that are mapped to pathways
		'''
	
		filteredSnpGeneMap = {}
		filteredGeneSnpMap = {}
		
		# include only SNPs mapped to pathways
		print '\nCreating filtered snp to gene mappings...'
		nMappedSNPs = 0
		nMappedGenes = 0
		print 'removing SNPs not mapped to pathways'
		for snp in snpGeneMap:
			if snp in pathSNPmap.snpPathMap:
				filteredSnpGeneMap[snp] = snpGeneMap[snp]
				nMappedSNPs += 1
		print 'filtered snpGeneMap with', nMappedSNPs, 'SNPs created (', len(snpGeneMap) - nMappedSNPs, 'SNPs removed)'
	

		# filter out genes not mapped to pathways
		print 'removing genes not mapped to pathways'
		nRemovedGenes = 0
		for snp in filteredSnpGeneMap:
			
			genesToRemove = []
			for gene in filteredSnpGeneMap[snp]:
			
				if gene not in genePathMap: # unmapped gene
					genesToRemove.append(gene)
					nRemovedGenes += 1
				else: # gene mapped to pathway
					# create gene to SNP map
					if gene in filteredGeneSnpMap:
						filteredGeneSnpMap[gene].append(snp)
					else:
						filteredGeneSnpMap[gene] = [snp]
						nMappedGenes += 1
			for delGene in genesToRemove:
				filteredSnpGeneMap[snp].remove(delGene)
		
		print nRemovedGenes, 'genes removed from snpGeneMap'
		print nMappedGenes, 'genes that map to genotyped SNPs and pathways remain'
		print 'filtered GeneSnpMap with', nMappedGenes,'genes created'
	
		self.filteredSnpGeneMap = filteredSnpGeneMap
		self.filteredGeneSnpMap = filteredGeneSnpMap