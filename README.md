# gghybrid
R package for analysis of hybrids and hybrid zones.

Functions should be run in the following order:
1. read.data  #Read in a data file in structure format or similar#
2. data.prep #Prepare the data for hybrid index and genomic cline estimation#
3. esth #Estimate hybrid indices#
4. plot_h  #Plot hybrid indices (optional)#
5. ggcline  #Estimate genomic clines#
6. plot_clinecurve  #Plot fitted clines (optional)#
7. compare.models #Compare two models (either from esth or ggcline) run on the same data set using the widely applicable information criterion (optional)#
