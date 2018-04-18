# gghybrid
R package for analysis of hybrids and hybrid zones. Currently includes hybrid index and genomic cline estimation for bi-allelic genomic data.

To cite package ‘gghybrid’ in publications use:

  Richard Ian Bailey (2018). gghybrid: Evolutionary Analysis of Hybrids and
  Hybrid Zones. R package version 0.0.0.9000.

Example code is in the file 'Example_code_RIBailey_17April2018.R ', and the data set used by this code is 'RB_Italy_ggcline_precol_headers_haploidND2.csv ', which was downloaded and prepared from here: https://www.datadryad.org/resource/doi:10.5061/dryad.v6f4d

Examples of plots produced by the package are 'cline_curve_plot.pdf ' and 'hybrid_index_plot.pdf '.

Basic functionality of the package is to read in SNP data in the form of structure files or similar, prepare data for analysis, carry out Bayesian MCMC hybrid index and genomic cline estimation, compare models (for either hybrid index or genomic clines) run on the same data set using the widely applicable information criterion (waic), and make plots of hybrid indices or cline curves.

gghybrid can be downloaded from within R using the following two lines of code:

install.packages("devtools"); devtools::install_github("ribailey/gghybrid")


Functions should be run in the following order:
1. read.data #Read in a data file in structure format or similar#
2. data.prep #Prepare the data for hybrid index and genomic cline estimation#
3. esth #Estimate hybrid indices#
4. plot_h #Plot hybrid indices (optional)#
5. ggcline #Estimate genomic clines#
6. plot_clinecurve #Plot fitted clines (optional)#
7. compare.models #Compare two models (either from esth or ggcline) run on the same data set using the widely applicable information criterion (optional)#

For usage please see help files for individual functions by typing e.g. '?read.data'.

Synopsis:
To understand mechanisms of speciation and the evolutionary impacts of admixture it is vital to identify loci showing restricted or biased introgression among hybridizing taxa. Genomic cline analysis provides a means to do this, by examining patterns of introgression of loci into foreign genomic backgrounds. Here I present a new R package, gghybrid, which allows hypothesis-testing on bi-allelic genomic data through Bayesian hybrid-index (proportion of allele copies coming from one of two parental reference sets) and logit-logistic genomic-cline (Fitzpatrick 2013) estimation. The package takes structure files or similar data tables as input, allows filtering of loci based on parental allele frequencies, and pooling and fixing of parameters followed by model comparison for both hybrid index and genomic clines with the Bayesian widely applicable information criterion (waic). It therefore provides great flexibility in comparing, for example, populations, transects, genomic regions or gene networks for differing patterns of admixture and introgression. It also allows rapid creation of a genotype table, with genotypes scored according to the parent-of-origin of each allele, and contains plot functions for hybrid index and genomic cline estimates. I use an adaptive algorithm during burnin to optimize multivariate parameter proposal distributions, utilizing both the acceptance rate and the estimated parameter covariance matrix. Furthermore, given the intention for the package to be used on large whole-genome data sets, I employ recursive estimation of posterior distributions to avoid storage of the full set of posterior values and hence improve memory efficiency.

Hybrid index estimation uses the formulae from Buerkle (2005), plus a prior:

Buerkle, C. A. (2005). Maximum likelihood estimation of a hybrid index based on molecular markers. Molecular Ecology Notes, 5(3), 684-687.

Genomic cline analysis uses the logit-logistic cline function of Fitzpatrick (2013), with parameter 'centre' instead of 'u' (see '?ggcline'):

Fitzpatrick, B. M. (2013). Alternative forms for genomic clines. Ecology and evolution, 3(7), 1951-1966.

References for the data set:

Hermansen JS, Haas F, Trier CN, Bailey RI, Nederbragt AJ, Marzal A, Sætre G (2014) Hybrid speciation by sorting of parental incompatibilities in Italian sparrows. Molecular Ecology 23(23): 5831-5842. https://doi.org/10.1111/mec.12910.

Hermansen JS, Haas F, Bailey RI, Nederbragt AJ, Trier CN, Marzal A, Sætre G (2014) Data from: Hybrid speciation by sorting of parental incompatibilities in Italian sparrows. Dryad Digital Repository. https://doi.org/10.5061/dryad.v6f4d.
