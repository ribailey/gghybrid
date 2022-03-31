# gghybrid
R package for analysis of hybrids and hybrid zones. Currently includes hybrid index and genomic cline estimation for bi-allelic genomic data.

Note: New version 2.0.0 20 March 2022.

#############################################################################################################
###Important updates for version 2.0.0 causing compatibility issues with code written for earlier versions###
#############################################################################################################
1. In 'read.data' the number of individuals must now be specified using the 'NUMINDS' option. The number 
   of loci does not need to be specified, although a warning will be produced if it is left blank.
2. Column names have changed for the genomic cline results in the ggcline$gc output object: 
 - 'v_mean' (old - the best posterior estimate for cline v) is now 'exp_mean_log_v', to reflect the fact that ggcline estimates log(v), then calculates best posterior v and its credible intervals at the end.
 - 'centre_mean' (old - the best posterior estimate for cline centre) is now 'invlogit_mean_logit_centre', to reflect the fact that ggcline estimates logit(centre), then calculates best posterior centre and its credible intervals at the end.
3. The new 'plot_clinecurve' function no longer plots data, only the cline curve. However, individual genotypes can be added to the plot. Furthermore, samples from the posterior for genomic cline estimates can be taken using the new function 'rtmvnormDT3', a cline calculated for each, and these can be added to the plot to indicate uncertainty. Please see examples in the new help file using '?plot_clinecurve'.
#############################################################################################################
#############################################################################################################
#############################################################################################################

To cite package ‘gghybrid’ in publications use:

Richard Ian Bailey. (2022). ribailey/gghybrid: gghybrid R package for Bayesian hybrid index and genomic cline estimation (v2.0.0). Zenodo. DOI: https://doi.org/10.5281/zenodo.3676498

Latest example code is in the file 'Example_code_RIBailey.R ', and the accompanying data set was downloaded and prepared from here: https://www.datadryad.org/resource/doi:10.5061/dryad.v6f4d

Examples of plots produced by the package are 'Figure_gcexample.pdf ' and 'Figure_hiexample.pdf '. All software comparisons and data simulation and subsequent analysis are included as R code files.

Basic functionality of the package is to read in SNP data in the form of structure files or similar, prepare data for analysis, carry out Bayesian MCMC hybrid index and genomic cline estimation, compare models (for either hybrid index or genomic clines) run on the same data set using the widely applicable information criterion (waic) or AIC, and make plots of hybrid indices or cline curves.

gghybrid can be downloaded from within R using the following two lines of code:

install.packages("devtools"); devtools::install_github("ribailey/gghybrid")


Functions should be run in the following order:
1. read.data #Read in a data file in structure format or similar#
2. data.prep #Prepare the data for hybrid index and genomic cline estimation#
3. split_data_prep #split the prepared data file into multiple subfiles, by (any number of) individual, locus or any other chosen character or factor column (optional)#
4. esth #Estimate hybrid indices#
5. plot_h #Plot hybrid indices (optional)#
6. ggcline #Estimate genomic clines#
7. plot_clinecurve #Plot fitted clines (optional)#
8. rtmvnormDT3 #Sample from the posterior distribution of fitted genomic clines for one or more test subjects (optional)#
9. compare.models #Compare two models (either from esth or ggcline) run on the same data set using the widely applicable information criterion (optional)#
10. calc_AIC #Calculate AIC for one esth or ggcline model#

For usage please see help files for individual functions by typing e.g. '?read.data'.

Synopsis:
To understand mechanisms of speciation and the evolutionary impacts of admixture it is vital to identify loci showing restricted or biased introgression among hybridizing taxa. Genomic cline analysis provides a means to do this, by examining patterns of introgression of loci into foreign genomic backgrounds. Here I present the R package gghybrid which allows hypothesis-testing on bi-allelic genomic data through Bayesian hybrid-index (proportion of allele copies coming from one of two parental reference sets) and logit-logistic genomic-cline (Fitzpatrick 2013) estimation. The package takes structure files or similar data tables as input, allows filtering of loci based on parental allele frequencies, and pooling and fixing of parameters followed by model comparison for both hybrid index and genomic clines with the Bayesian widely applicable information criterion (waic) or AIC. It therefore provides great flexibility in comparing, for example, populations, transects, genomic regions or gene networks for differing patterns of admixture and introgression. It also allows rapid creation of a genotype table, with genotypes scored according to the parent-of-origin of each allele, and contains plot functions for hybrid index and genomic cline estimates. I use an adaptive algorithm during burnin to optimize multivariate parameter proposal distributions, utilizing both the acceptance rate and the estimated parameter covariance matrix. Furthermore, given the intention for the package to be used on large whole-genome data sets, I employ recursive estimation of posterior distributions to avoid storage of the full set of posterior values and hence improve memory efficiency, and also provide a function to split the data analysis file unto multiple sub-files.

Hybrid index estimation uses the formulae from Buerkle (2005), plus a prior:

Buerkle, C. A. (2005). Maximum likelihood estimation of a hybrid index based on molecular markers. Molecular Ecology Notes, 5(3), 684-687.

Genomic cline analysis uses the logit-logistic cline function of Fitzpatrick (2013), with parameter 'centre' instead of 'u' (see '?ggcline'):

Fitzpatrick, B. M. (2013). Alternative forms for genomic clines. Ecology and evolution, 3(7), 1951-1966.

References for the data set:

Hermansen JS, Haas F, Trier CN, Bailey RI, Nederbragt AJ, Marzal A, Sætre G (2014) Hybrid speciation by sorting of parental incompatibilities in Italian sparrows. Molecular Ecology 23(23): 5831-5842. https://doi.org/10.1111/mec.12910.

Hermansen JS, Haas F, Bailey RI, Nederbragt AJ, Trier CN, Marzal A, Sætre G (2014) Data from: Hybrid speciation by sorting of parental incompatibilities in Italian sparrows. Dryad Digital Repository. https://doi.org/10.5061/dryad.v6f4d.
