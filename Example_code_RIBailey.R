#############################################################################
#############################################################################
####gghybrid example code####################################################
#############################################################################
#############################################################################

#############################################################################
#Richard Ian Bailey 11 March 2022############################################
#############################################################################

#The processes for hybrid index and genomic cline estimation are split into 
#multiple stages. 

#(1) Data are initially loaded from a file in rectangular table or structure file format 
#using the read.data function, and converted into a table with full header row and
#one column per locus, one row per test subject per allele copy.

#(2) The data.prep function then produces an object ready for analysis, with 
#one row per non-missing allele copy. This allows for flexible filtering, pooling and 
#splitting (for example for parallelization) of the data in preparation for analysis.

#(3) For both hybrid index and genomic cline estimation, the parameters for each 
#test subject (typically 'INDLABEL' and 'locus' respectively) are estimated independently, 
#hence allowing parallelization down to individual test subjects. The function split_data_prep 
#is included to automate the splitting and storage of data files ready for parallelized 
#analysis.

#(4) Hybrid index estimation is carried out independently of genomic cline 
#estimation. This allows different sets of loci to be used for each analysis if desired. 
#For example, the user may wish to use only LD-pruned loci for hybrid index 
#estimation (this is typically recommended for admixture estimation methods that 
#assume statistical independence of loci), but may wish to include linked sets of loci 
#in genomic cline estimation to increase resolution along chromosomes.

#(5) The data.prep object, along with the resulting object from hybrid 
#index estimation, are used as inputs for genomic cline estimation. As with hybrid 
#index estimation, the data.prep object can be filtered and split as required 
#prior to analysis.

#(6) Both hybrid index and genomic cline estimation can be set up to carry out 
#statistical model comparison using waic or AIC.

#(7) The function rtmvnormDT3 is provided for sampling from the joint posterior of 
#the two genomic cline parameters across many test subjects simultaneously, 
#using the output object from ggcline.

#############################################################################
#############################################################################
#############################################################################

#The typical procedure for hybrid index estimation is to upload a data object 
#that includes a set of 'TEST' samples that are putatively admixed, and 
#two sets of 'parental reference' samples, each of which is intended to 
#represent one of the two non-admixed populations/species that hybridized to 
#create the test populations. There is no lower limit to how differentiated the 
#parental reference populations have to be (as long as there is SOME differentiation), 
#so 'hybridization' can be among populations within species. The two parental reference 
#or 'source' samples (referred to as 'S0' and 'S1' in gghybrid) are used to estimate parental 
#allele frequencies, and these frequencies are then used to estimate the 
#hybrid index (proportion of alleles inherited from S1) and its uncertainty 
#in each test subject (also in the S0 and S1 individuals themselves, if 
#desired).

#Alternatively, the user has the choice to import parental 
#reference (S0 and S1) allele frequencies, without the requirement for 
#importing any S0 or S1 samples. This may be useful if: (a) the user wishes 
#to save memory, can calculate/estimate parental allele frequencies in another 
#software, and doesn't intend to include parental reference samples in genomic 
#cline estimation; (b) the user wishes to test custom parental allele frequencies based 
#on a hypothesis, rather than based on real population samples.


#############################################################################
#############################################################################
#Loading data using the 'read.data' function#################################
#############################################################################
#############################################################################

#Existing data files in any of the structure file formats can be loaded, by 
#specifying 'mainparams' information in the same way as for structure 
#software (see documentation). A rectangular data table with complete header row 
#or no header row can also be loaded.

#WARNING! gghybrid currently can't accept zeroes as missing data, 
#due to a bug (which has now been actioned and the fix will hopefully be implemented soon) 
#in the underlying data.table package. Structure files created in PLINK using 
#'--recode structure' have missing data coded as zero, so for now these zeroes need to 
#be converted to 'NA' before reading in PLINK structure files with 'read.data'. 

#It may often be useful to upload a marker info file alongside the SNP data, 
#with one row per locus and any number of columns with potentially useful 
#information on the markers. ***This is essential if no parental reference 
#samples are included in the dataset (see below)***. The marker info will be 
#joined to the prepared data analysis object, if desired, when running the second 
#function, 'data.prep', and can then be used e.g. for splitting or pooling loci 
#in different ways to test hypotheses.

##############################################################
##############################################################
#Example code#################################################
##############################################################
##############################################################


##############################################################
#Loading a data file that includes parental reference samples#
##############################################################

#Load the package#
library(gghybrid)
library(coda)#Used for calculating the Gelman diagnostic to compare multiple posterior samples for convergence#

#Documentation.
?read.data


#First create example data files and save to the working directory#

#1. Genomic (SNP) data as a regular table with full header row#

ex <- "INDLABEL,POPID,chr1:001,chr1:002\nind1,pop1,A,A\nind1,pop1,A,B\nind2,pop1,B,B\nind2,pop1,B,B\nind3,pop2,NA,A\nind3,pop2,B,A\n"
ex

#Save 'ex' to the working directory as 'ex.data', to be read in by read.data. Open in e.g. notepad or notepad++ to view in table format#

cat("INDLABEL POPID chr1:001 chr1:002","ind1 pop1 A A","ind1 pop1 A B","ind2 pop1 B B","ind2 pop1 B B","ind3 pop2 NA A","ind3 pop2 B A", file = "ex.data", sep = "\n")

#2. Now create and save a file, ex2.data, in PLINK structure format (but with missing data recoded from 0 to NA), with a MAPDISTANCE row present below the marker names#

cat("chr1:001 chr1:002","-3 3 2 4","ind1 pop1 A A A B","ind2 pop1 B B B B","ind3 pop2 NA B A A", file = "ex2.data", sep = "\n")


#Load ex.data, in this case without specifying the number of loci (useful if this is not known in advance). These are the minimum arguments required for a file in this format#

dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,MISSINGVAL=NA)#NUMLOCI is not specified so will be calculated, with a warning#

#Take a look. The object contains the data plus several other useful outputs - see the documentation for explanations#

dat

#The same but specifying the number of loci. The resulting object is identical#

dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA)#No warning this time as NUMLOCI is specified#
dat


#The minimum arguments required to read in the PLINK structure format file with missing data recoded as 'NA' (gives a warning because NUMLOCI is not specified)#

dat2 <- read.data(file="ex2.data",MISSINGVAL=NA,NUMINDS=3,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
dat2


#Same, but with all arguments including those where the default works for this file format#

dat2 <- read.data(
  file="ex2.data",
  mainparams = NULL,              #Default#
  extracol.names = NULL,          #Default#
  precol.headers = 0,
  nprecol=2,
  markername.dup = 0,             #Default#
  NUMLOCI.autoAccept = TRUE,      #Default, see below#
  EXTRACOL = 0,                   #Default#
  INDLABEL = 1,                   #Default#
  LOCDATA = 0,                    #Default#
  MAPDISTANCE = 1,
  MARKERNAME = 1,                 #Default#
  MISSINGVAL = NA,                #Cannot be zero for the time being#
  NUMINDS = 3,
  NUMLOCI = 2,                    #Only mandatory when PHASED=1; otherwise will be calculated internally if not known#
  ONEROW = 1,
  PHASED = 0,                     #Default#
  PHENOTYPE = 0,                  #Default#
  PLOIDY = 2,                     #Default#
  POPFLAG = 0,                    #Default#
  POPID = 1,                      #Default#
  RECESSIVEALLELE = 0,            #Default#
  marker.info.file = NULL,        #Default#
  sourceAbsent = FALSE            #Default#
)

dat2


#If NUMLOCI is not specified, the number of loci will be calculated with a warning, which will 
#not interfere with downstream processes. However, if 'NUMLOCI.autoAccept = FALSE' is set, the 
#user is required to manually accept the calculated number of loci. This option is included in 
#case the user wishes to verify that the calculated number of loci is accurate. Stops with an 
#error if the estimated NUMLOCI is manually rejected. Example:

dat2 <- read.data(file="ex2.data",NUMLOCI.autoAccept = FALSE,MISSINGVAL=NA,NUMINDS=3,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)



#############
#Marker info#
#############

#Example: Create a marker info file indicating whether the locus is intronic or exonic, and save to the working directory#

cat("locus type","chr1:001 intronic","chr1:002 exonic", file = "ex_marker_info.data", sep = "\n")

#Read it in alongside the SNP data#

dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info.data")

dat#The marker info is loaded into the object, in preparation for joining to the SNP data when running 'data.prep'.


######################################################################
#Loading a data file that does not include parental reference samples#
######################################################################

#For the situation when parental (S0 and S1) allele frequencies are to be loaded, but not 
#parental reference samples#

#NOTE: It's usually beneficial to include parental reference samples, especially because steep 
#clines that occur outside the range of the test samples, such as between the hybrid index of 
#test samples and parental reference samples, can only be accurately estimated if parental 
#reference samples are not only included in the dataset, but their hybrid index estimated and 
#they are subsequently included in genomic cline analysis.

#A marker.info.file can be loaded regardless of whether S0 and S1 samples are present in the 
#dataset (see above), but it is obligatory if they are not present#

#When 'sourceAbsent = TRUE', as a minimum the following columns with the exact headers in the 
#first set of speech marks below are required in the marker.info.file (more columns are allowed).

cat("locus refAllele alternateAllele S0.prop_r S1.prop_r type","chr1:001 A B 0.1 0.9 intronic","chr1:002 B A 0.8 0.2 exonic", file = "ex_marker_info2.data", sep = "\n")

#"...prop_r" means the allele frequency of the reference allele. Choice of reference and 
#alternate alleles is arbitrary.

dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info2.data", sourceAbsent = TRUE)

dat#The data.prep function described below will then determine which allele has higher frequency in S1, as it would if parental reference samples were included#



#############################################################################
#############################################################################
#Preparing data for analysis using the 'data.prep' function##################
#############################################################################
#############################################################################

?data.prep

###############################################################
#Using a dataset with no parental reference (S0 or S1) samples#
###############################################################

#Read in a file with no parental reference samples.

dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info2.data", sourceAbsent = TRUE)

dat

#Run data.prep to produce the data analysis table and some other useful outputs.


prepdata=data.prep(
 data=dat$data,               #part of the read.data output object#
 loci=dat$loci,               #part of the read.data output object#
 sourceAbsent = TRUE,         #default is FALSE#
 marker.info=dat$marker.info, #must be specified if sourceAbsent = TRUE; make sure it contains the headers and their contents as specified above when 'sourceAbsent=TRUE'#
 alleles=dat$alleles,         #part of the read.data output object#
 #S0,                         #leave blank in the absence of parental reference samples#
 #S1,                         #leave blank in the absence of parental reference samples#
 precols=dat$precols,         #part of the read.data output object#
 return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1 than S0#
 return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values calculated within this function#
)

#Take a look.

prepdata

names(prepdata$locus.data)

###########################################################################
#Make a new data input file with more entries, including S0 and S1 samples#
###########################################################################

#Create the input file, read it in with read.data, then run data.prep#

#Create ex3.data in PLINK structure format with missing data recoded from 0 to NA and with a MAPDISTANCE row present below the marker names#

cat("chr1:001 chr1:002",
 "-3 3 2 4",
 "ind1 pop1 A A A B","ind2 pop1 B B B B","ind3 pop2 NA B A A", 
 "ind4 pop3 A A A A","ind5 pop3 A B A A","ind6 pop3 A A A A","ind7 pop3 NA A A A","ind8 pop4 A NA B B","ind9 pop4 B B B B","ind10 pop4 B B B B","ind11 pop5 A A B B",
 file = "ex3.data", sep = "\n")

#pop3 represents S0; pop4 and pop5 represent S1#

#Read in the data#

dat3 <- read.data(file="ex3.data",MISSINGVAL=NA,NUMINDS=11,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
dat3

prepdata= data.prep(
 data=dat3$data,
 loci=dat3$loci,
 sourceAbsent = FALSE,        #default is FALSE, for the existence of parental reference samples in the dataset#
 #marker.info=dat$marker.info,#this time we're not uploading any extra marker info#
 alleles=dat3$alleles,
 S0="pop3",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
 S1=c("pop4","pop5"),                          #second parental reference set; must be specified when sourceAbsent = FALSE#
 precols=dat3$precols,
 return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1#
 return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including those already uploaded plus values calculated within this function#
)

#prepdata$locus.data has more columns in the presence of parental reference samples, because uncertainty in parental allele frequencies can be calculated 
#by assuming a Bayesian beta-distributed posterior. '...shape1' and '...shape2' are the two parameters of the beta distribution for the posterior.

prepdata


#########################################
#Now completely exclude the POPID column#
#########################################

#Just to show that the functions work in the absence of this column#

#First with no parental reference samples#

#modify the plink format 'ex3.data' above to remove the POPID column, and use the existing marker info file, 'ex_marker_info2.data'#

cat("chr1:001 chr1:002",
 "-3 3 2 4",
 "ind1 A A A B","ind2 B B B B","ind3 NA B A A", 
 "ind4 A A A A","ind5 A B A A","ind6 A A A A","ind7 NA A A A","ind8 A NA B B","ind9 B B B B","ind10 B B B B","ind11 A A B B",
 file = "ex4.data", sep = "\n")

#you must specify 'POPID=0' (meaning the column is absent)#

dat4 <- read.data(file="ex4.data",POPID=0,nprecol=1,NUMINDS=11,NUMLOCI=2,MISSINGVAL=NA,marker.info.file="ex_marker_info2.data",sourceAbsent=TRUE,ONEROW=1,MAPDISTANCE=1,precol.headers=0)

dat4

#Now prepare the data for analysis#

prepdata= data.prep(
 data=dat4$data,
 loci=dat4$loci,
 sourceAbsent = TRUE,         #default is FALSE#
 marker.info=dat4$marker.info,#must be specified if sourceAbsent = TRUE; make sure it contains the headers and their contents as specified above for 'sourceAbsent=TRUE'#
 alleles=dat4$alleles,
 #S0,                         #leave blank in the absence of parental reference samples#
 #S1,                         #leave blank in the absence of parental reference samples#
 POPID.name = "INDLABEL",     #***must be specified if there's no POPID column, to avoid an error***#
 precols=dat4$precols,
 return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1#
 return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including those already uploaded plus values calculated within this function#
)

prepdata

#Now including parental reference samples#

#Create and save a marker info file.

cat("locus type","chr1:001 intronic","chr1:002 exonic", file = "ex_marker_info.data", sep = "\n")

#Read data in again first, but with a marker info file ('ex_marker_info.data') 
#not including the obligatory columns for sourceAbsent=TRUE. Plus set sourceAbsent=FALSE so that the code doesn't search for the obligatory columns and throw up an error.

dat5 <- read.data(file="ex4.data",POPID=0,nprecol=1,NUMINDS=11,NUMLOCI=2,MISSINGVAL=NA,marker.info.file="ex_marker_info.data",sourceAbsent=FALSE,ONEROW=1,MAPDISTANCE=1,precol.headers=0)

dat5

#Now prepare the analysis object#

prepdata= data.prep(
 data=dat5$data,
 loci=dat5$loci,
 sourceAbsent = FALSE,         #default is FALSE#
 marker.info=dat5$marker.info, #this now contains only the locus (obligatory) and type columns#
 alleles=dat5$alleles,
 S0=c("ind4","ind5","ind6","ind7"),                  #must be specified when sourceAbsent = FALSE#
 S1=c("ind8","ind9","ind10","ind11"),                #must be specified when sourceAbsent = FALSE#
 POPID.name = "INDLABEL",                            #***must be specified if there's no POPID column; the S0 and S1 sample references must be present in the column specified here***#
 precols=dat5$precols,
 return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1#
 return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including those already uploaded plus values calculated within this function#
)

prepdata

#################################################################################
#Locus filtering based on parental reference allele frequencies and sample sizes#
#################################################################################

#Take a look at the filtering options using '?data.prep':

#The defaults (no filtering, but still removes loci with exactly zero allele frq difference
#between S0 and S1):
 # max.S.MAF = 0.5,
 # min.diff = 0,
 # min.allele.copies.S0 = 0,
 # min.allele.copies.S1 = 0,
 # AF.CIoverlap = TRUE,

#Now we switch to using a real dataset with more loci to demonstrate these options#

#The file 'RB_Italy_ggcline_precol_headers_haploidND2.csv' is included in the gghybrid download from github. Make sure it's in the working directory, or change the file path below#

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569)

dat

#The file contains data for 77 diploid markers for 569 male sparrow individuals#

#Run data.prep including the different filtering options. Try hashing out or changing these 
#options to see how it affects the number of loci in the resulting analysis table (see 
#prepdata$Nloci.postfilter in the data.prep object).

prepdata=data.prep(data=dat$data,
 loci=dat$loci,
 alleles=dat$alleles,
 S0=c("Kralove","Oslo"), #POPID names for the first parental reference set#
 S1="LesinaSPANISH",     #POPID names for the second parental reference set#
 precols=dat$precols,
###Filtering below###
 max.S.MAF = 0.2,               #minor allele frq must be below this value in at least one of S0 or S1; loci with smallest MAF among parental reference sets above this will be removed#
 min.diff = 0.2,                #loci with parental allele frequency difference < this value will be removed#
 min.allele.copies.S0 = 30,     #loci with fewer non-missing allele copies in the S0 parental reference set will be removed#
 min.allele.copies.S1 = 12,     #in this dataset the S1 sample size is much smaller than S0, so I'm being less strict with filtering for sample size#
 AF.CIoverlap = FALSE,          #***RECOMMENDED*** filtering option IF parental reference samples are included - removes all loci for which there is overlap between S0 and S1 in the Bayesian posterior 95% credible intervals of allele frequency# 
###Filtering above###
 return.genotype.table=T,
 return.locus.table=T)

prepdata

#I recommend using 'AF.CIoverlap = FALSE' because parental allele frequencies are estimated 
#with error, and this option will retain only those loci for which we are at least certain the 
#allele frequency difference is real and in the expected direction. For small allele frequency 
#differences with strongly ovelapping posteriors for example, the true difference could be in 
#the opposite direction, which affects assignment of alleles to S1 and S0 and subsequent hybrid 
#index and genomic cline estimation.

#Try using 'AF.CIoverlap = FALSE' in the absence of parental reference samples, and read the 
#resulting error message. The user must supply the credible/confidence intervals in the 
#marker.info object in this case.

#We use one of the existing fake-data objects#

prepdata= data.prep(
 data=dat4$data,
 loci=dat4$loci,
 sourceAbsent = TRUE,
 marker.info=dat4$marker.info,
 alleles=dat4$alleles,
 POPID.name = "INDLABEL", 
 precols=dat4$precols,
 AF.CIoverlap = FALSE,         #Without the appropriate entries in the marker.info object this attempt at filtering will cause an error#
 return.genotype.table=TRUE,
 return.locus.table=TRUE 
)

######################################################################################################################
######################################################################################################################
#Dividing up the data.prep object in preparation for parallel runs of either hybrid index or genomic cline estimation#
######################################################################################################################
######################################################################################################################

?split_data_prep

#Especially useful for users with large datasets and access to a server allowing 
#parallelization. Also useful for running individuals or loci serially on a standalone 
#computer if memory is not sufficient to run all at once.

#This function splits the data.prep analysis object into multiple files, each with equal 
#numbers of test subjects, except for the final file, the size of which depends on how many 
#test subjects remain (a 'test subject' is e.g. 'INDLABEL' for hybrid index estimation on each 
#individual, or 'locus' for genomic cline estimation on individual loci).

#The resulting objects are written to the working directory as csv files and removed from the 
#workspace. These can then be loaded into R and used as the data.prep.object in downstream 
#functions. csv files are quite memory-heavy but load quickly.

#Make a set of data.prep files, one for each individual.

split_data_prep(
 data.prep.object=prepdata$data.prep,   #The data analysis table#
 splitBy="INDLABEL",                    #The planned test.subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
 keepN=1                                #The number of test.subjects you want in each resulting file#
)

#Make another set of files, one for each set of 10 loci.

split_data_prep(
 data.prep.object=prepdata$data.prep,   #The data analysis table#
 splitBy="locus",                       #The planned test.subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
 keepN=10                               #The number of test.subjects you want in each resulting file (the final file will contain fewer if the total is not a multiple of this number)#
)


#The entry for the 'splitBy' option will be included in the filename for each resulting file. 
#Therefore, to create a list of files for analysis e.g. of the 'locus' files above:

files <- list.files(pattern="_locus_");

files

#############################################################################
#############################################################################
#Running hybrid index estimation#############################################
#############################################################################
#############################################################################

?esth

#Use the sparrow data with loci filtered by parental allele frequency credible interval overlap.

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569)
prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE)

###############################################
#In the presence of parental reference samples#
###############################################

hindlabel=esth(
 data.prep.object=prepdata$data.prep,
 read.data.precols=dat$precols,
 include.Source=TRUE,	                 #Leave at default TRUE if you want hybrid indices for the parental reference individuals, which is often useful#
plot.ind = c("PD07-254","PD07-160","PD07-159","PI07-110","PI08-498","PH08-285"),  #Optionally plot some individuals in real time. Merely shows how well the adaptive burnin is working#
plot.col = c("blue","green","cyan","purple","magenta","red"),
 nitt=3000,                              #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)


#There will be a warning if any NAs are produced for the best posterior estimate 
#(h_posterior_mode).

setkey(hindlabel$hi,beta_mean)#beta_mean is the mean of the posterior and is always >0 and <1, which is sometimes useful, but the mode is the best posterior estimate in the presence of a prior#
hindlabel#Run this line twice if hindlable$hi doesn't show up the first time#

#If you want to see all the results#

View(hindlabel$hi)

#It is possible to get a result of NA or NaN for the posterior hybrid index estimate 
#(h_posterior_mode). If this happens it probably means the variance of the proposal 
#distribution during burnin was too high. The default is to use the variance 
#(1/(n_allele_copies_per_test_subject/2)/10 after the first 100 iterations. If this fails 
#for any individuals, add another zero and set the option 'init.var2' to this value. And/or 
#use the MCMC plots within esth to try and identify the problem. I often find that if I run an 
#individual a second time even without changing the settings, it works fine.

#nitt=3000, burnin=1000 should be sufficient for convergence of the posterior estimates across 
#multiple runs according to the Gelman-Rubin diagnostic, which should have a value < 1.2 
#to indicate convergence. Both nitt and burnin can be increased if the G-R diagnostic indicates 
#some individuals with poor convergence across runs. It's unlikely more than nitt=5000, 
#burnin=2000 would be necessary, and even nitt=2000 is very likely to be sufficient.

#Run esth twice more on the same data and carry out the Gelman-Rubin diagnostic test#

hindlabel2=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,include.Source=TRUE,nitt=3000,burnin=1000)
hindlabel3=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,include.Source=TRUE,nitt=3000,burnin=1000)

#Join the three results objects#

setkey(hindlabel$hi,INDLABEL);setkey(hindlabel2$hi,INDLABEL);setkey(hindlabel3$hi,INDLABEL)
hindall=hindlabel$hi[,.(Source,INDLABEL,POPID,beta_shape1,beta_shape2)][hindlabel2$hi[,.(INDLABEL,beta_shape1,beta_shape2)]][hindlabel3$hi[,.(INDLABEL,beta_shape1,beta_shape2)]]

hindall[,rn:=seq(1,.N)];

#Take a sample of N=10000 from the posterior of each run for each individual, logit transform 
#(the G-R diagnostic assumes normally distributed posteriors), and calculate the Gelman-Rubin 
#diagnostic per individual.

hindall[,h_gelman:=coda::gelman.diag(
 mcmc.list(
  mcmc(qlogis(rbeta(10000,beta_shape1,beta_shape2))), 
  mcmc(qlogis(rbeta(10000,i.beta_shape1,i.beta_shape2))),
  mcmc(qlogis(rbeta(10000,i.beta_shape1.1,i.beta_shape2.1)))
          )
                              )$psrf[1],
 by=rn];

#Plot. Ideally all values should be < 1.2#

hist(hindall$h_gelman)

#Take a look at which are most extreme. With this dataset it's always those with hi estimates 
#very close to 0 or 1. This may be an issue with the logit-transform of these beta-distributed 
#posterior samples being asymmetrical and therefore breaking the assumptions of the G-R 
#diagnostic, rather than a problem with convergence.

hindall[order(h_gelman)]

#For example, plot a histogram of the logit-transformed posterior samples from the individual 
#with the most extreme G-R value#

setkey(hindall,h_gelman);hindall
hist(hindall[569,qlogis(rbeta(10000,beta_shape1,beta_shape2))])


###################################################
#esth in the absence of parental reference samples#
###################################################

#We will first create a marker info file in the correct format to be loaded alongside the data.
#The marker info file should include parental allele frequencies when there are no parental 
#reference individuals in the dataset. These parental allele frequencies are needed for both 
#hybrid index and genomic cline estimation.

#For this example, first load and prepare the data, and include producing a 'locus.data' table 
#with the data.prep function. This locus data table will then be modified, saved, and used as 
#the marker info file in the absence of declared parental reference samples.

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",
 precols=dat$precols,AF.CIoverlap = FALSE,return.locus.table = TRUE);

#Create a marker info file for this dataset in the correct format for sourceAbsent = TRUE and 
#save to the working directory. Column names for the parental allele frequency (either allele 
#can be chosen) must match those specified second within 'setnames' below.

sparrowMarkers=prepdata$locus.data[,.(locus,S0_allele,S1_allele,S0.prop_0,S1.prop_0)];
setnames(sparrowMarkers,c("S0_allele","S1_allele","S0.prop_0","S1.prop_0"),c("refAllele","alternateAllele","S0.prop_r","S1.prop_r"));

fwrite(sparrowMarkers,"sparrowMarkers.csv")#By default 'fwrite' saves as a comma-separated file, whatever the file extension#

#Now run the workflow without specifying parental reference samples#

#Read in the sparrow data again, this time alongside the new marker info file#

dat <- read.data(file="RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,NUMINDS=569,MISSINGVAL=NA,marker.info.file = "sparrowMarkers.csv", sourceAbsent = TRUE);

#Run data.prep for the absence of parental reference samples and the same filtering as before. 
#We don't have uncertainty measures for parental allele frequencies, so this time we filter by 
#difference in parental allele frequency using 'min.diff' (<=0.2).

prepdata= data.prep(
 data=dat$data,
 loci=dat$loci,
 sourceAbsent = TRUE,        #Required when no parental reference samples are declared#
 marker.info=dat$marker.info,
 alleles=dat$alleles,
 precols=dat$precols,
 min.diff = 0.2
)


#esth can now be run as above.

hindlabelx=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,nitt=3000,burnin=1000)

hindlabelx$hi[order(beta_mean)]#beta_mean is the posterior mean (the mode is the best posterior estimate), but allows for more refined sorting of estimates#

#################################################################
#Estimating a hybrid index per locus across all TEST individuals#
#################################################################

#The cline centre (see ggcline below) is an estimate of the genome-wide hybrid index at which 
#the focal locus allele frequency is halfway between the parental allele frequencies. It 
#therefore does not clearly indicate the frequency of the S1 allele across all TEST 
#individuals, which also depends on cline steepness (v). Furthermore, for very extreme cline 
#centres (e.g. centre = 1e-6), the cline shape becomes distorted and and a cline with v < 1 
#(typically thought of as a shallow cline) can appear steep (see ggcline below). Hence, the 
#clearest way to estimate the S1 allele frequency across a set of individuals is either to 
#run genomic cline analysis with v fixed at 1 (see ggcline function below) or to estimate a 
#per-locus hybrid index for the pooled individuals, or both.

#Use a dataset including parental reference samples, and then run per-locus hybrid index 
#estimation on the 'TEST' individuals only (i.e. excluding S0 and S1), to get the proportion 
#of S1 allele copies across TEST individuals only.

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);

hindlabel=esth(
 data.prep.object=prepdata$data.prep,
 read.data.precols=dat$precols,
 test.subject="locus",                                       #This switches the analysis from estimating h-index per individual (INDLABEL, the default) to per locus#
 include.Source=FALSE,	                                     #Exclude the parental reference individuals (which are included by default) to get the h-index per locus for TEST individuals only#
plot.ind = c("A2M","ACO1","CHD1Z","NFIL3","SECISBP2","HECTD1"),  #Optional#
plot.col = c("blue","green","cyan","purple","magenta","red"),    #Optional#
 nitt=3000,                                                  #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)

#The TEST sample includes some Spanish and house sparrows (the parental species) so we don't 
#expect any really extreme values in this case.
View(hindlabel$hi[order(beta_mean)])

##############################################################################
#Running esth on one of the files produced by the 'split_data_prep' function##
##############################################################################

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);

#Produce one data analysis file per individual and save to WD.

split_data_prep(
 data.prep.object=prepdata$data.prep,   #The data analysis table#
 splitBy="INDLABEL",                    #The planned test.subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
 keepN=1                                #The number of test.subjects you want in each resulting file#
)

#Take the first one and estimate hybrid index.

prep=fread("prepdata_INDLABEL_1.csv")

hindlabel=esth(
 data.prep.object=prep,           #The analysis object#
 read.data.precols=dat$precols,
plot.ind = prep[1,1],
plot.col = c("purple"),
 nitt=3000,                       #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)



#################################################################
#################################################################
#Model comparison using waic and AIC for hybrid index estimation#
#################################################################
#################################################################

?compare.models

######################################################################################
#Compare the fit between running esth at the individual level versus population level#
######################################################################################

#This will produce a separate comparison for each different POPID (the level with lower 
#resolution).

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);

#First run with test.subject="INDLABEL". For both runs we need to store the likelihood 
#calculations for each data point using 'return.likmeans=TRUE'.

hindlabel=esth(
 data.prep.object=prepdata$data.prep,     #The analysis object#
 read.data.precols=dat$precols,
 test.subject = "INDLABEL",               #default#
 test.subject.compare = "POPID",          #***the test.subject you want to compare with the current test.subject, which only requires an entry when running the finer resolution test.subject***#
 include.Source=TRUE,	                  #Set to TRUE if you want hybrid indices for the parental reference individuals#
 return.likmeans = TRUE,                  #***likelihood calculations required for model comparison***#
 nitt=3000,                               #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)

hindlabel

#Then test.subject="POPID"#

hindlabel2=esth(
 data.prep.object=prepdata$data.prep,      #The analysis object#
 read.data.precols=dat$precols,
 test.subject = "POPID",                   #***Changed to "POPID"***#
 test.subject.compare = "INDLABEL",        #***leave this blank when running esth for the lower resolution model (POPID in this case), otherwise you'll get a warning message***#
 include.Source=TRUE,	                   #Set to TRUE if you want hybrid indices for the parental reference individuals#
 return.likmeans = TRUE,                   #***likelihood calculations required for model comparison***#
 nitt=3000,                                #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)

hindlabel2


####################
#Now compare models#
####################

#####
#AIC#
#####

#AIC is calculated per declared test subject, so for direct comparison, use the lower-resolution 
#test subject to calculate AIC for both models.

aic_INDLABEL=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel,
 test.subject="POPID")

aic_POPID=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel2,
 test.subject="POPID")

#Now join the tables and calculate the difference in AIC.

setkey(aic_INDLABEL,POPID);setkey(aic_POPID,POPID)
aic_both=aic_INDLABEL[aic_POPID]#The columns named 'i...' refer to the second table, aic_POPID

aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (i.e. they mean that AIC was higher across all individuals when hybrid index was estimated at the individual level)

aic_both[order(aicdiff)]

#########################
#Now the waic comparison#
#########################

#In this case the model comparison is done within the function.

comp1=compare.models(esth.object1=hindlabel,esth.object2=hindlabel2);

comp1[order(waicdiff_npar_AICscale)]

#The values and their differences are scaled to be directly comparable to the AIC results.

#Positive waicdiff_npar_AICscale values indicate stronger support for the simpler model.

#Negative values support the model with more parameters, and values close to zero indicate 
#no clear support for either model, as with AIC.

#By far the strongest evidence for differentiation among individuals within a population is 
#for the Lago di Burano population, which appears to contain an individual recently admixed 
#with Spanish sparrows:

hindlabel$hi[POPID=="LagodiBurano"]

#A quick comparison of AIC versus waic model comparison results.

setkey(aic_both,POPID);setkey(comp1,POPID)
comp2=aic_both[comp1]

plot(comp2[,waicdiff_npar_AICscale]~comp2[,aicdiff])
abline(a=0,b=1,col="red",lty=2)
abline(h=-4,v=-4,col="grey",lty=2)

#In this example, while the points with difference < -4, i.e. strongly in favour of the more 
#complex model (grey lines) mostly correspond, waic seems to otherwise more strongly favour 
#the simpler model.


##########################
##########################
#Genomic cline estimation#
##########################
##########################

#Genomic cline analysis estimates 2 parameters, v and centre:

#centre = the genome-wide hybrid index at which the focal test.subject allele frequency is 
#halfway between those of the parental reference sets.

#v = cline steepness. It would be 1/cline width IF this were geographic cline analysis 
#(which has an unconstrained x axis of distance along a transect). With the x axis constrained 
#to [0,1], i.e. the set of genome-wide hybrid index estimates, v is not exactly 1/width but 
#remains a measure of cline steepness.

?ggcline

#See above for preparing all the objects needed for genomic cline analysis.

#The full set of arguments for ggcline including all defaults (comments indicate which entries 
#are not defaults).

gc=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = TRUE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  fix.subject.v = FALSE,
  fix.value.v,
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("A2ML1_SNP1","GTF2H2"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre = c(0, sqrt(50)),     #Default#
  prior.logv = c(0, sqrt(10)),             #Default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)

#The main output, gc$gc, includes rounded parameter estimates on the original scale 
#(exp_mean_log_v and invlogit_mean_logit_centre) and their 95% credible intervals, and the 
#(not rounded) posterior parameter means (mean_log_v and mean_logit_centre) and multivariate 
#normally distributed covariance matrix on the latent scale (var_log_v, var_logit_centre, 
#cov_log_v_logit_centre).

gc

#ggcline includes Bayesian p values for the cline parameters but does not include false 
#discovery rate or other adjustments, which are left to the user. There are several R packages 
#for FDR calculations.

#Be aware that occasionally p is estimated as exactly 0 due to very narrow credible intervals, 
#and this may cause FDR calculations to fail. 

#Using 'beta_mean' as the h-index estimate instead of h_posterior_mode. 
#This would be appropriate if you had set 'prior=c(0,0)' (i.e. no prior) when running esth.


gc1=ggcline(
  data.prep.object=prepdata$data.prep,  #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  esth.colname="beta_mean",             #***If the column containing the hybrid index values is not called 'h_posterior_mode', it must be entered here***#
  include.Source = TRUE,                #Default is FALSE#
  read.data.precols=dat$precols,        #Needs an entry#
  plot.test.subject = c(1,2),#c("CHD1Z","HSDL2"),
  plot.col = c("blue","red"),
  plot.ylim = c(-3, 5),
  plot.pch.v.centre = c(1, 3),
  nitt=5000,                            #Needs an entry#
  burnin=2000                           #Needs an entry#
)


##########################
#Try another test.subject#
##########################

#Let's try estimating a single parameter for cline steepness, v, across a whole chromosome, 
#while still allowing cline centre to be estimated individually per locus on that chromosome.

#Due to issues with rapid vectorized sampling of the proposal distribution when pooling occurs 
#by group rather than across the whole dataset, if you want to for example estimate a single v 
#across all loci per chromosome for multiple chromosomes, the chromosomes must be run 
#separately (sequentially or in parallel).

#First, add a column to prepdata$data.prep indicating the chromosome for each marker.

chrom=fread("markerchr.csv")#This contains only the loci retained in prepdata$data.prep after filtering#

setkey(prepdata$data.prep,locus);setkey(chrom,locus)
prepdata$data.prep=prepdata$data.prep[chrom]

#Make a data.prep object for the Z chromosome markers only#

prepZ=prepdata$data.prep[chrom=="Z"]

unique(prepZ$locus)

#First run ggcline as before, but Z chromosome only. It's necessary to match this dataset to the 
#run with v pooled across loci (see below) in order to correctly calculate the total number 
#of parameters for model comparison.

gc_unpooled_z=ggcline(
  data.prep.object=prepZ,               #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = TRUE,
  return.likmeans = TRUE,               #***Important to set this to TRUE for model comparison***#
  read.data.precols=dat$precols,        #Needs an entry#
  fix.subject.v = FALSE,
  fix.value.v,
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("ADFP","ZCCHC6"),
  plot.col = c("orange","cyan"),
  plot.ylim = c(-1, 2),
  plot.pch.v.centre = c(1, 3),
  prior.logitcentre = c(0, sqrt(20)),
  prior.logv = c(0, sqrt(5)),
  nitt=5000,                           #Needs an entry#
  burnin=2000,                         #Needs an entry#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)

gc_unpooled_z


#Now run ggcline for this Z chromosome-only dataset with the v parameter pooled across all loci. 
#This can then be compared for model fit with the Z-linked loci run individually above.

#The MCMC is less efficient with pooling, so using longer nitt and burnin.

#Defaults hashed out#

gc_poolv_Z=ggcline(
  data.prep.object=prepZ,               #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  #esth.colname="h_posterior_mode",
  #test.subject = "locus",                 
  poolv = TRUE,                         #***Set to TRUE***#
  #poolcentre = FALSE,                  #Keep the default, to estimate a different cline centre per locus#
  include.Source = TRUE,
  return.likmeans = TRUE,               #***Important to set this to TRUE for model comparison***#
  read.data.precols=dat$precols,        #Needs an entry#
  #fix.subject.v = FALSE,
  #fix.value.v,
  #fix.subject.centre = FALSE,
  #fix.value.centre,
  plot.test.subject = c("ADFP","ZCCHC6"),
  plot.col = c("orange","cyan"),
  plot.ylim = c(-2, 2),
  plot.pch.v.centre = c(1, 3),
  #prior.logitcentre = c(0, sqrt(20)),
  #prior.logv = c(0, sqrt(5)),
  nitt=12000,                          #Use a longer post-burnin (5000 in this case)#
  burnin=7000,                         #Use a longer burnin#
  #start.v = NULL,
  #start.centre = NULL,
  #init.var.v = NULL,
  #init.var.centre = NULL,
  #init.cov.vcentre = NULL,
  #print.k = 50
);

gc_poolv_Z

#Use compare.models.

?compare.models

gclinecomp=compare.models(
 ggcline.object1=gc_unpooled_z,
 ggcline.object2=gc_poolv_Z,
 ggcline.pooled=TRUE              #Default is FALSE. This must be set to TRUE when one of the ggcline runs includes pooling#
)

#Always gives a warning.

gclinecomp

#The result, a large negative value of waicdiff_npar_AICscale, indicates strong support for 
#the more complex model with independent estimates of v per locus.

#calc_AIC can also be used as for esth.

aic_unpooled_z=calc_AIC(
 data.prep.object=prepZ,
 esth.object=hindlabel,
 ggcline.object=gc_unpooled_z,
 test.subject="locus",
 ggcline.pooled=TRUE
)

aic_poolv_Z=calc_AIC(
 data.prep.object=prepZ,
 esth.object=hindlabel,
 ggcline.object=gc_poolv_Z,
 test.subject="locus",
 ggcline.pooled=TRUE
)

aic_unpooled_z$AIC - aic_poolv_Z$AIC

#Similar difference to waic.

################################################
#A model comparison for ggcline with no pooling#
################################################

#First run ggcline with both parameters fixed to their null values. Only 2 iterations necessary.

gc_v_centre_fixed=ggcline(
  data.prep.object=prepdata$data.prep,  #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  include.Source = TRUE,
  return.likmeans = TRUE,
  read.data.precols=dat$precols,        #Needs an entry#
  fix.subject.v = gc$gc$locus,          #***Vector of locus names***#
  fix.value.v = 1,                      #***The null parameter value on the data scale***#
  fix.subject.centre = gc$gc$locus,     #***Vector of locus names***#
  fix.value.centre = 0.5,               #***The null parameter value on the data scale***#
  nitt=2,                               #MCMC estimation not needed when both parameters are fixed for all test subjects, only nitt=2 is necessary to avoid errors#
  burnin=0                              #MCMC estimation not needed when both parameters are fixed for all test subjects#
)

gc_v_centre_fixed

#Now run model comparison against 'gc' above, which has independent v and centre estimates 
#for each locus.

gclinecomp_fixed=compare.models(
 ggcline.object1=gc_v_centre_fixed,
 ggcline.object2=gc,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

gclinecomp_fixed[order(waicdiff_npar_AICscale)]

#Few of the loci show strong support for the null model (positive 'waicdiff_npar_AICscale' 
#values above about 2), and many show strong support for the 2-parameter model.

#How does this compare with p values? Remember that the above comparison involved fixing both 
#parameters.

setkey(gclinecomp_fixed,locus);setkey(gc$gc,locus)
gcomp=gclinecomp_fixed[gc$gc]

#Red points are centre_pvalue on the y axis; blue crosses are v_pvalue. So each locus has a 
#red point and blue cross with identical values on the x axis.

plot(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],type="n",ylab="-log10(p value) centre and v")
points(gcomp[,-log10(v_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="blue",pch=4)
points(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="red")
abline(h=-log10(0.05),col="grey",lty=2)                   #Grey line indicates p = 0.05#


#A closer look at those borderline significant, indicated by the grey lines 
#(x axis > -4 and y axis < 1.3 is relatively weak evidence for deviation from null).

plot(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],type="n",xlim=c(-50,5),ylim=c(0,40),ylab="-log10(p value) centre and v")
points(gcomp[,-log10(v_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="blue",pch=4)
points(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="red")
abline(h=-log10(0.05),v=-4,col="grey",lty=2)                   #Grey line indicates p = 0.05#


#The comparison would be cleaner if we only fixed one parameter:

gc_v_fixed=ggcline(
  data.prep.object=prepdata$data.prep,  #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  include.Source = TRUE,
  return.likmeans = TRUE,
  read.data.precols=dat$precols,        #Needs an entry#
  plot.test.subject = c("ADFP","ZCCHC6"),
  plot.col = c("orange","cyan"),
  plot.ylim = c(-1, 2),
  plot.pch.v.centre = c(1, 3),
  fix.subject.v = gc$gc$locus,          #***Vector of locus names***#
  fix.value.v = 1,                      #***The null parameter value on the data scale (could also be fixed to some other value)***#
  #fix.subject.centre = gc$gc$locus,    #Not fixing centre this time#
  #fix.value.centre = 0.5,              #Not fixing centre this time#
  nitt=5000,
  burnin=2000
)


gclinecomp_vfixed=compare.models(
 ggcline.object1=gc_v_fixed,
 ggcline.object2=gc,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

gclinecomp_vfixed[order(waicdiff_npar_AICscale)]


setkey(gclinecomp_vfixed,locus);setkey(gc$gc,locus)
gcomp=gclinecomp_vfixed[gc$gc]

#This time, the blue crosses (v parameter) should more closely match the waic results.

plot(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],type="n",ylab="-log10(p value) centre and v")
points(gcomp[,-log10(v_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="blue",pch=4)
points(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="red")
abline(h=-log10(0.05),col="grey",lty=2)                   #Grey line indicates p = 0.05#

#A closer look at those borderline significant, indicated by the grey lines 
#(x axis > -4 and y axis < 1.3 is relatively weak evidence for deviation from null).

plot(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],type="n",xlim=c(-50,5),ylim=c(0,40),ylab="logit(p value) centre and v")
points(gcomp[,-log10(v_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="blue",pch=4)
points(gcomp[,-log10(centre_pvalue)]~gcomp[,waicdiff_npar_AICscale],col="red")
abline(h=-log10(0.05),v=-4,col="grey",lty=2)                   #Grey line indicates p = 0.05#


#####
#AIC#
#####

#Try v fixed against nothing fixed.

aic_v_fixed=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel,
 ggcline.object=gc_v_fixed,
 test.subject="locus",
ggcline.pooled=FALSE
)

aic_no_fixed=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel,
 ggcline.object=gc,
 test.subject="locus",
ggcline.pooled=FALSE
)

#Now join the tables and calculate the difference in AIC for each locus.
setkey(aic_no_fixed,locus);setkey(aic_v_fixed,locus);
aic_both=aic_no_fixed[aic_v_fixed];#The columns named 'i...' refer to the second table, aic_v_fixed#
aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (v fixed at 1)#

#Show the results in order of aicdiff for each POPID. Those at the top of the ordered table 
#more strongly favour the individual-level model.
aic_both[order(aicdiff)]

############################
#The chromosome-level model#
############################

#Run ggcline using chromosome as test.subject rather than locus.

gc_by_chr=ggcline(
  data.prep.object=prepdata$data.prep,  #Needs an entry#
  esth.object=hindlabel,                #Needs an entry#
  test.subject = "chrom",               #***Changed from default***#
  include.Source = TRUE,
  return.likmeans = TRUE,
  read.data.precols=dat$precols,        #Needs an entry#
  plot.test.subject = c("1","Z"),
  plot.col = c("orange","cyan"),
  plot.ylim = c(-1, 2),
  plot.pch.v.centre = c(1, 3),
  nitt=5000,                            #Needs an entry#
  burnin=2000                           #Needs an entry#
)

gc_by_chr

#The column "chrom" is absent from the 'gc' results object and must first be added to 
#both $gc and $likmeans.

chrom=fread("markerchr.csv")#This contains only the loci retained in prepdata$data.prep#

setkey(gc$gc,locus);setkey(chrom,locus);setkey(gc$likmeans,locus)
gc$gc=gc$gc[chrom]
gc$likmeans=gc$likmeans[chrom]

#Now compare models at the chromosome level.

gclinecomp_chr=compare.models(
 ggcline.object1=gc,
 ggcline.object2=gc_by_chr,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

gclinecomp_chr[order(waicdiff_npar_AICscale)]

#Again, negative values for waicdiff_npar_AICscale favour the more complex model 
#(independent cline parameters per locus within each chromosome).

#####
#AIC#
#####

#The column "chrom" is absent from the 'gc' results object and must first be added to gc$gc.
chrom=fread("markerchr.csv")#This contains only the loci retained in prepdata$data.prep#

setkey(gc$gc,locus);setkey(chrom,locus);setkey(gc$likmeans,locus)
gc$gc=gc$gc[chrom]
gc$likmeans=gc$likmeans[chrom]

#Check if there's a column called 'chrom' in prepdata$data.prep (there will be if you've run
#all the previous code). If it's absent, add it here.

setkey(prepdata$data.prep,locus);setkey(chrom,locus);
prepdata$data.prep=prepdata$data.prep[chrom]


aic_gc=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel,
 ggcline.object=gc,
 test.subject="chrom",
ggcline.pooled=FALSE
)

aic_gc_by_chr=calc_AIC(
 data.prep.object=prepdata$data.prep,
 esth.object=hindlabel,
 ggcline.object=gc_by_chr,
 test.subject="chrom",
ggcline.pooled=FALSE
)

#Now join the tables and calculate the difference in AIC for each locus.

setkey(aic_gc,locus);setkey(aic_gc_by_chr,locus);
aic_both=aic_gc[aic_gc_by_chr];#The columns named 'i...' refer to the second table, aic_v_fixed#
aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (v fixed at 1)#

#Show the results in order of aicdiff for each POPID. Those at the top of the ordered table 
#more strongly favour the individual-level model.
aic_both[order(aicdiff)]


########################################################
########################################################
#Sampling from the genomic cline posterior distribution#
########################################################
########################################################

#Take 500 samples from the posterior parameter distribution of each of 5 loci#

#The function to do this is called 'rtmvnormDT3'. It is designed to simultaneously sample 
#multiple two-dimensional multivariate (truncated) normal distributions from a data.table. 
#The number 3 is because there are two internal functions in gghybrid called 'rtmvnormDT' and 
#'rtmvnormDT2', both designed for sampling from the parameter proposal distribution during MCMC 
#within ggcline.

#rtmvnormDT3 can be used more generally to sample from a set of bivariate normal distributions, 
#including truncated. It produces a long-form data.table with one row per sample per group.

#The function must be run from within a data.table, in the square brackets after the comma.

gcsamp=gc$gc[locus%in%c("ACO1","CLTA","LNPEP","HSDL2","MCCC2"),  #in this case I'm subsetting by row to only include these 5 loci#
 rtmvnormDT3(
  nsamp=500,                                 #default is 1 sample per locus#
  meanlogv=mean_log_v,                       #no default, column must be specified#
  meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
  varlogv=var_log_v,                         #no default, column must be specified#
  varlogitcentre=var_logit_centre,           #no default, column must be specified#
  covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
  lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
  upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
            ),                               #end of rtmvnormDT3 function#
 by="locus"];                                #indicate a grouping variable for the sampling#

gcsamp

#Plot the posterior samples for one of these loci on the latent scale#

plot(gcsamp[locus=="ACO1",log_v]~gcsamp[locus=="ACO1",logit_centre])


#Same plot this time on the original scale (i.e plotting the actual values of v and centre)#

plot(gcsamp[locus=="ACO1",exp(log_v)]~gcsamp[locus=="ACO1",plogis(logit_centre)])



################################################################
################################################################
#Running the Gelman diagnostic on multiple runs of ggcline######
################################################################
################################################################

#To save memory, gghybrid does not store individual posterior samples, but instead stores 
#their means, variances and covariances. 

#If multiple independent runs are carried out with the same data set, these can be tested for 
#convergence by sampling post hoc from the posterior distribution.

#Run ggcline two more times, no pooling or fixing.

gcrun2=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  include.Source = TRUE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  plot.test.subject = c("A2ML1_SNP1","GTF2H2"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  nitt=5000,
  burnin=2000
)

gcrun3=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  include.Source = TRUE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  plot.test.subject = c("A2ML1_SNP1","GTF2H2"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  nitt=5000,
  burnin=2000
)

#Take 1000 posterior samples across all loci for all three runs (including the earlier 'gc').

gc_samp=gc$gc[,
 rtmvnormDT3(
  nsamp=1000,                                 #default is 1 sample per locus#
  meanlogv=mean_log_v,                       #no default, column must be specified#
  meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
  varlogv=var_log_v,                         #no default, column must be specified#
  varlogitcentre=var_logit_centre,           #no default, column must be specified#
  covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
  lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
  upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
            ),                               #end of rtmvnormDT3 function#
 by="locus"];                                #indicate a grouping variable for the sampling#
#
gcrun2_samp=gcrun2$gc[,
 rtmvnormDT3(
  nsamp=1000,                                 #default is 1 sample per locus#
  meanlogv=mean_log_v,                       #no default, column must be specified#
  meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
  varlogv=var_log_v,                         #no default, column must be specified#
  varlogitcentre=var_logit_centre,           #no default, column must be specified#
  covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
  lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
  upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
            ),                               #end of rtmvnormDT3 function#
 by="locus"];                                #indicate a grouping variable for the sampling#
#
gcrun3_samp=gcrun3$gc[,
 rtmvnormDT3(
  nsamp=1000,                                 #default is 1 sample per locus#
  meanlogv=mean_log_v,                       #no default, column must be specified#
  meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
  varlogv=var_log_v,                         #no default, column must be specified#
  varlogitcentre=var_logit_centre,           #no default, column must be specified#
  covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
  lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
  upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
            ),                               #end of rtmvnormDT3 function#
 by="locus"];                                #indicate a grouping variable for the sampling#
#

#Calculate the gelman diagnostic for log v for each locus.

gelman_res=data.table(locus=gc$gc$locus)

for(i in 1:length(gelman_res$locus)){
loci=gelman_res[i,locus];
gelman_res[i,gelman:=coda::gelman.diag(mcmc.list(mcmc(gc_samp[locus==loci]$log_v),mcmc(gcrun2_samp[locus==loci]$log_v),mcmc(gcrun3_samp[locus==loci]$log_v)))$psrf[1]]

}#endfor#

#Plot the results.

hist(gelman_res$gelman)#All values should be < 1.2#



#######################
#######################
#Plotting hybrid index#
#######################
#######################

?plot_h

abc = plot_h(data=hindlabel$hi,
 test.subject=hindlabel$test.subject,
 mean.h.by="POPID",			             #Calculate the mean hybrid index for each value of the "POPID" column#
 sort.by=c("mean_h","POPID","h_posterior_mode"),  #Order test subjects along the x axis by the mean hybrid index calculated above and also by individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
 col.group="POPID",
 group.sep="POPID",
 fill.source=TRUE,
 basic.lines=FALSE,
 source.col=c("blue","red"),
 source.limits=c("blue","red"),
 custom.abline=abline(h=0.76,col="grey",lty=2),  #The lowest hybrid index estimate for (Sardinian) Spanish sparrows#
 cex=1,pch=16,
 cex.lab=1.5,cex.main=1.5,ylim=c(0,1))


#Stretch and shape the plot as you want, then add a legend using the plot_h object 'abc' 
#created above.

 setkey(abc,rn);		      #Order data by row number#

 legend("topleft",	      #Place the legend in the top left of the figure#
 abc[,POPID], 		      #Name of the field by which data point colours are grouped#
 bg="white",			#Background colour#
 text.col=c("black"),         #Text colour#
 pch=22, 				#Text size#
 col=abc[,col.Dark2],         #Name of the field containing colour information#
 pt.bg=abc[,col.Dark2],	      #Name of the field containing colour information#
 ncol=5,				#Number of columns for splitting the group names#
 cex=0.7, pt.cex=0.7);



#########################
#########################
#Plotting genomic clines#
#########################
#########################

?plot_clinecurve


#The code below plots the cline for the locus CLTA, then adds the data points to the plot in 
#the form of diploid genotypes scaled to 0-1, 1 being homozygote for the allele with higher 
#frequency in parental reference S1. It also adds clines for 1000 random samples from the parameter 
#posteriors, as a form of shading.


pdf(file="mycline.pdf")#This results in a nice smooth line with the plot saved in the working directory#
########################################################################
#Plot a genomic cline for one locus and add uncertainty and data points#
########################################################################

#Plot the curve for locus CLTA, using a previously created ggcline object (see examples for ggcline).
plot_clinecurve(
ggcline.object=gc$gc,
cline.locus="CLTA",
locus.column="locus",
cline.col="#E495A5",
cline.centre.line="CLTA",
cline.centre.col="black"
)

#Add a title and axis labels.
title(main = "CLTA",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)

###################################################
#Add sampled cline curves to represent uncertainty#
###################################################

#There is no option for shading credible intervals within plot_clinecurve, so here is a suggestion for how to add some: 
#sample parameter values from the joint posterior, and add a grey curve to the plot for each posterior sample. 
#This will overlay the original cline curve plotted above, so it needs to be plotted again.

#First generate random samples from the posterior using gghybrid's rtmvnormDT3 function.

gcsamp=gc$gc[locus%in%c("CLTA"),
 rtmvnormDT3(
  nsamp=1000,                                #default is 1 sample per locus#
  meanlogv=mean_log_v,                       #no default, column must be specified#
  meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
  varlogv=var_log_v,                         #no default, column must be specified#
  varlogitcentre=var_logit_centre,           #no default, column must be specified#
  covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
  lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
  upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
            ),                               #end of rtmvnormDT3 function#
 by="locus"];                                #indicate a grouping variable for the sampling#

#Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
locipost=gc$gc[locus=="CLTA",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
gcsamp=rbind(gcsamp,locipost)

#Add the individual curves to the plot in a 'for' loop.
for(i in 1:nrow(gcsamp)){
        v = gcsamp[i,exp(log_v)]
        u = gcsamp[i,logit_centre*exp(log_v)]
        S1.prop_1 = gc$gc[locus=="CLTA",S1.prop_1];
        S0.prop_1 = gc$gc[locus=="CLTA",S0.prop_1];

        par(new=T);

if(i < nrow(gcsamp)){
        curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
	        from=0,to=1,axes=F,xlab="",ylab="", col="grey",lwd=0.5,ylim=c(0,1));
}else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
        curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
	        from=0,to=1,axes=F,xlab="",ylab="", col="#E495A5",lwd=3,ylim=c(0,1));
};
                         };#end of for loop#

######################################################
#Add data to the plot as genotypes scaled from 0 to 1#
######################################################

ploidy=2

#Create a data.table containing individual reference, hybrid index, locus name and genotype. See the ggcline examples for creating prepdata and hindlabel.
setkey(prepdata$data.prep,INDLABEL);setkey(hindlabel$hi,INDLABEL);
genodat=hindlabel$hi[,.(INDLABEL,h_posterior_mode)][prepdata$data.prep[locus=="CLTA",sum(Source_allele)/ploidy,by=c("INDLABEL","locus")]]

#Add the resulting points to the cline plot.
points(genodat[,V1]~genodat[,h_posterior_mode],pch=16,cex=0.7,col="#E495A5")
###
dev.off()#END OF PLOT CODE FOR LOCUS CLTA************************************************


######################################################
#Plot cline curves for multiple loci on the same plot#
######################################################

#No data points or uncertainty included this time.

plot_clinecurve(
ggcline.object=gc$gc,
cline.locus=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
locus.column="locus",
cline.col=c("orange","blue","green","red","magenta"),
null.line.locus=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
null.line.col=c("orange","blue","green","red","magenta"),
cline.centre.line=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
cline.centre.col=c("orange","blue","green","red","magenta")
)

#Add a title and axis labels:
title(main = "ACO1,EGR1,A2M,HSDL2,RPS4",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)
###


#It looks much nicer saved as a PDF!

