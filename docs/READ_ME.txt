PATHWAYS SPARSE REDUCED-RANK REGRESSION

Distributed under GNU general public licence (see COPYRIGHT.txt) 
Copyright (C) 2012 Matt Silver

WARRANTY DISCLAIMER
Please note there is no warranty for the program, to the extent permitted by applicable law.
Except when otherwise stated in writing the copyright holders and/or other parties provide the program
“as is” without warranty of any kind, either expressed or implied, including, but not limited to, 
the implied warranties of merchantability and fitness for a particular purpose. The entire risk
as to the quality and performance of the program is with you. Should the program
prove defective, you assume the cost of all necessary servicing, repair or correction.

Please send comments or bugs to matt.silver@lshtm.ac.uk

If you use PsRRR in a published paper, please cite at least one of the following 2 papers:

PATHWAYS ANALYSIS WITH UNIVARIATE Y
Silver, Matt, and Montana, G. “Fast Identification of Biological Pathways Associated with
a Quantitative Trait Using Group Lasso with Overlaps.” 
Statistical Applications in Genetics and Molecular Biology 11, no. 1 (2012).

PATHWAYS ANALYSIS WITH MULTIVARIATE Y
Silver, Matt, Eva Janousova, Xue Hua, Paul M. Thompson, and Giovanni Montana. 
“Identification of Gene Pathways Implicated in Alzheimer’s Disease Using Longitudinal 
Imaging Phenotypes with Sparse Regression.” 
NeuroImage 63, no. 3 (2012)
******************************************************************************************

1. OVERVIEW

PsRRR consists of 3 separate stages, each with its own code:

1. mapping SNPs to pathways 
2. adapting weights to remove bias, e.g. due to variation in number of genes per pathway or LD between SNPs 
3. final sparse regression analysis


2. SYSTEM REQUIREMENTS

The code is written in Python, and tested using the Enthought Python Distribution (EPD)
7.2-2. However, it should run under Python 2.7.  The program requires that the Numpy and
SciPy modules are installed.

Processor and memory requirements will depend on the size of the datasets to be analysed.


3. INSTALLATION 

Copy the code and all subdirectories to an appropriate directory on your
system.  Provided that your system meets the above requirements, it should then be ready to
run.  


4. DIRECTORY STRUCTURE

Code and associated data is arranged in the following way:

~/								'home' directory for program code and all associated files

~/docs/							installation info, program documentation etc

~/data/							source data files - genotypes, phenotypes, pathway mapping files, gene locations

~/src/preprocessData/			python source code for mapping SNPs to pathways

~/src/PsRRR/  					python source code for PsRRR analysis

~/src/params/						location of parameter file(s), detailing programme parameters, 
								file names and locations etc.


In addition, the following output directories are created at runtime, depending on the location
specified by the user in the params file.  So, in the example included with this code, 
where the output directory is 'sample_output', the following directory structure is created:

~/sample_output/					contains output file subdirectories
~/sample_output/pathwayMapping/		files created during SNP to pathway mapping process 
~/sample_output/adaptWeights/		files created while adapting weights
~/sample_output/adaptWeights/plots/	plots to visualise weight adaptation process
~/sample_output/results/			results files


5. RUNNING THE PIPELINE WITH SAMPLE DATA PROVIDED

Sample genotype, phenotype and pathway data is included with the code.  See docs/pathway_mapping.txt,
docs/adapting_weights.txt and docs/PsRRR.txt for details on how to run the code with this data.
