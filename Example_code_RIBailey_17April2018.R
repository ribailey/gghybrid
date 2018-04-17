##############################################################################
###########gghybrid example code##############################################
##############################################################################

############Richard Ian Bailey 17 April 2018##################################

#Install the package from GitHub#
install.packages("devtools")
devtools::install_github("ribailey/gghybrid")

#Attach it#
library(gghybrid)

#Take a look at the function help files#
?read.data	#Read in a data file in structure format or similar#
?data.prep	#Prepare the data for analysis#
?esth		#Run hybrid index estimation#
?plot_h	#Plot estimated hybrid indices and credible intervals#
?ggcline	#Run genomic cline estimation#
?plot_clinecurve	#Plot one or more fitted cline curves, with or without individual data points#
?compare.models	#compare two models run on the same data set using the widely applicable information criterion#

#(Note: gghybrid relies on the data.table package for data manipulation)#


#Example: Using a data file in structure format but with a complete header row.
#The file contains 1 column per marker (equivalent to ONEROW=0 in structure) for a set
#of diploid markers, plus haploid (mitochondrial) ND2. For ND2, the second (non-existent)
#allele copy is coded as missing data. Other file formats would require more information
#when running 'read.data' (see ?read.data).

#The data file contains two mandatory columns to the left of the marker columns, named
#INDLABEL (the individual references) and POPID (the sample population). These columns
#don't need to be named in the header row (or could have different names), but it makes it easier if they are.

#I downloaded the data (and then created a single input file in the right format)
#from Dryad here: https://www.datadryad.org/resource/doi:10.5061/dryad.v6f4d

#Read in the data file (This is the simplest file format for reading in)#

dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",
 nprecol=2,MISSINGVAL=NA)

#Take a look at the help file to get a handle on what the object 'dat' now contains.

###

#Data preparation and filtering. Here I'm filtering out loci that have a minor allele
#frequency greater than 0.1 in both parental reference sets. There are also options for
#filtering by difference in parental allele frequencies, and for number of allele copies
#in each parental reference set (variable among loci due to missing data).

#The function uses objects produced by 'read.data'#

prepdata=data.prep(data=dat$data,
    loci=dat$loci,
    alleles=dat$alleles,
    S0=c("Kralove","Oslo"), #POPID names for the first parental reference set#
    S1="LesinaSPANISH", #POPID names for the second parental reference set#
    precols=dat$precols,
    max.S.MAF = 0.1,	#Filtering by parental minor allele frequency#
    return.genotype.table=T,
    return.locus.table=T)

#'return.genotype.table=T' makes an optional table of genotypes, where for each locus
#an individual's genotype (assuming diploidy) will be 0 (two copies of the allele with
#relatively higher frequency in the 'S0' parental set), 1 (heterozygote), 2 (two copies
#of the designated 'S1' allele). This table isn't needed in downstream functions, but
#could be useful e.g. for estimating parental linkage disequilibria (associations of
#alleles from the same parent species).

#'return.locus.table=T' is also optional and not needed downstream. It's just a table
#with one row per marker, giving some information on parental allele frequencies, sample
#sizes etc.

###

#Next, run hybrid index estimation#

#This function uses objects used by both the previous functions#

hindlabel= esth(data.prep.object = prepdata$data.prep,
read.data.precols = dat$precols,
include.Source = TRUE,	#Set to true if you want hybrid indices for the parental reference individuals#
plot.ind = c("P?08-141","PD11-255","PH08-442","PI07-243",
        "PI08-501","PX08-520"),
plot.col = c("blue","green","cyan","purple","magenta","red"),
nitt=10000,burnin=5000)

#The plots ('plot.ind' and 'plot.col' options) are optional. They plot accepted posterior hybrid 
#index values in real time, in this case for 5 randomly chosen individuals.

#'esth' has more functionality - this above just shows the basics#

#Take a look at the results#

hindlabel

#data.tables sometimes have a strange habit of not showing up the first
#time - if that happens just run the above line again.

###

#Plot a subset of the estimated hybrid indices (the resulting object 'abc' is useful for making a legend)#

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

#
abc = plot_h(data=hindlabel$hi[c("Kralove","Susa","Figline","Accettura","Crotone","Pula","LesinaSPANISH")],#Subset of POPIDs#
test.subject=hindlabel$test.subject,
mean.h.by="POPID",			#Calculate the mean hybrid index for each value of the "POPID" column#
sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid 
								#index calculated above and also by individual hybrid index#
col.group="POPID",
group.sep="POPID",
fill.source=TRUE,
basic.lines=FALSE,
source.col=c("blue","red"),
source.limits=c("blue","red"),
cex=1,pch=16,
cex.lab=1.5,cex.main=1.5,ylim=c(0,1))
#

#Reshape the plot window as you want#

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)		#Order data by row number#
legend("topleft",	#Place the legend in the top left of the figure#
abc[,POPID], 		#Name of the field by which data point colours are grouped#
bg="white",			#Background colour#
text.col=c("black"), #Text colour#
pch=22, 				#Text size#
col=abc[,col.Dark2], #Name of the field containing colour information#
pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
ncol=2,				#Number of columns for splitting the group names#
cex=1, pt.cex=1)

###

#Run genomic cline analysis#

gc1=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=dat$precols,
    plot.test.subject=c("A2ML1_SNP2"),
    plot.col=c("orange"),
    plot.ylim=c(-2,3),
    nitt = 10000, burnin = 5000, print.k = 50)

#Again the real-time plots are optional. Open circles are parameter v (width), 
#'plus' symbols are logit(centre). centre ranges from 0 to 1 (the value of the hybrid index
#at which allele frequencies for the locus are halfway between the parental values).
#The null value for centre is 0.5, but logit(centre) is plotted here to 
#try and improve clarity. logit(centre) ranges from -Inf to Inf and the null value is 0.

#Take a look#

gc1

###

#Plot the cline curve for one locus, and add a title and axis labels#

plot_clinecurve(
    ggcline.object=gc1$gc,
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel$hi,
    cline.locus="CHD1Z",
    cline.col="green",
    null.line.locus="CHD1Z",	#A diagonal dashed line representing the null expectation#
    null.line.col="black",
    plot.data="CHD1Z",		#Plot the data too#
    data.col="black",
    plot.genotype=TRUE,		#Plot genotypes, scaled from 0 to 1, rather than alleles (which would all be 0 or 1)#
    PLOIDY=2,			#Needed when plotting genotypes#
    cline.centre.line="CHD1Z",#Plot the cline centre#
    cline.centre.col="green")

title(main = "CHD1Z",xlab="Hybrid index",
      ylab="Genotype frequency",cex.main=1.5,cex.lab=1.5)


#These functions have more options than shown here, and there are also possibilities for 
#pooling or fixing parameters in a variety of different ways,
#followed by model comparison of pairs of models using the 'compare.models' function.