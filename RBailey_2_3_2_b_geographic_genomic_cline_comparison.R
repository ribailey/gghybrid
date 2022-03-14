#############################################################################################################
#2.3.2(b) Code for comparison between gghybrid (genomic cline software) and hzar (geographic cline software)#
#############################################################################################################

####################################################################################################
#Richard Ian Bailey 14 February 2022################################################################
####################################################################################################
library(gghybrid)
library(hzar)
library(ggplot2)
library(cowplot)
library(colorspace)
####################################################################################################

#############################################
#Direct comparison of the two functions######
#############################################

#The cline centre parameter is equivalent between geographic and logit-logistic genomic cline 
#analysis, but because it is estimated on the logit scale in gghybrid, it is constrained to 
#between 0 and 1 on the data scale.

#Cline width in geographic cline analysis is defined as 1/maximum slope. The parameter v in 
#logit-logistic genomic cline analysis is equivalent to the maximum slope of a geographic cline, 
#after rescaling the x axis from [0,1] to [-Inf,Inf].

#Therefore, due to the x axis [0,1] constraint in genomic cline analysis, the true steepness of v 
#on the [0,1] scale changes somewhat with changing cline centre (it gets steeper as centre moves 
#away from 0.5).

#However, with cline centre = 0.5 and a steep cline fully encompassed within the [0,1] range 
#(e.g. v = 10; width = 1/10), the two clines are indistinguishable.

#Hence:

#The two curves are only (roughly) identical when centre = 0.5#

#The genomic cline gets steeper for the same v the closer the centre is to 0 or 1#

#########################################################################
#########################################################################
#Simulating haplotypes from the geographic and genomic cline ############
#functions followed by parameter estimation using both gghybrid and hzar#
#########################################################################
#########################################################################

##########################################################################
#Genomic data are simulated independently from each cline function #######
#(genomic and geographic) and for each sample size (10, 100, 1000, 10000)#
##########################################################################

#In all cases, width = 0.1, v = 10 (1/width), centre = 0.5, parentals fixed for alternate alleles, 
#individuals sampled from 21 populations evenly distributed between 0-1 (which represents either the 
#geographic location along a 1D transect, or hybrid index for geographic and genomic cline analysis 
#respectively).

############################################################################################
#############################sampsize 10####################################################
############################################################################################
pMin=0                                 #Parental reference assumed fixed#
pMax=1                                 #Parental reference assumed fixed#
centre=0.5
x=seq(0,1,length=21)                   #length = number of sampling sites#
sampsize=10                            #Total samples from all sites combined#

set.seed(12345)
x2=sample(x,sampsize,replace=T)

x2=sort(x2)
rm(x)

clineres=data.table(x=x2,centre=centre,width=0.1);
clineres[,v:=1/width][,u:=qlogis(centre)*(1/width)]#Preparing for sampling from the genomic cline function, with its parameters u and v#

set.seed(23456)
clineres[,ygeogsamp:=rbinom(sampsize,1,prob=pMin + (pMax - pMin) * (1/(1 + exp(-((x - centre) * 4/width)))))]   #Sampling from a geographic cline, using the same formula as hzar#
set.seed(34567)
clineres[,ygenomsamp:=rbinom(sampsize,1,prob=pMin + (x^v/(x^v + (1 - x)^v*exp(u)))*(pMax - pMin))]              #Sampling from a genomic cline, using the same formula as gghybrid#

#These are the population mean allele frequencies, needed for hzar#

ygeogmeans=clineres[,mean(ygeogsamp),by=x];setnames(ygeogmeans,"V1","ygeogsampmeans");
ygenommeans=clineres[,mean(ygenomsamp),by=x];setnames(ygenommeans,"V1","ygenomsampmeans");

#Open circles are population allele frequencies for simulated haplotypes from the geographic cline function; 
#pluses from the genomic cline function#

plot(ygeogmeans$ygeogsampmeans~ygeogmeans$x,ylim=c(0,1),xlim=c(0,1))
points(ygenommeans$ygenomsampmeans~ygenommeans$x,col="blue",pch=3)

#########################################################
#Preparing for GC (genomic cline analysis with gghybrid)#
#########################################################

#For this we need a data.prep.object (the genotype data ready for analysis, in the format produced by the 
#data.prep function) and an esth.object (the individual hybrid indices, in the format produced by the esth 
#function).

clineres[,INDLABEL:=paste("a",seq(1,.N),sep="")]
clineres[,POPID:=paste("a",x,sep="")]

clineres_10=clineres

#Making a custom data.prep object.
prepdata=list()
prepdata$data.prep=clineres[,.(INDLABEL,POPID,ygenomsamp,ygeogsamp)]
prepdata$data.prep[,Source_allele:=ygenomsamp]       #The 'Source_allele' column indicates the individual haplotype. Below I estimate genomic clines for the genomic cline function sample first#
prepdata$data.prep[,Source:="TEST"][,S0.prop_1:=pMin][,S1.prop_1:=pMax][,locus:="L1"]

#Making a custom esth object.
hindlabel=list()
hindlabel$hi=clineres[,.(x)]

setnames(hindlabel$hi,"x","h_posterior_mode")

hindlabel$hi[,INDLABEL:=prepdata$data.prep[,INDLABEL]]
hindlabel$hi[,POPID:=prepdata$data.prep[,POPID]]
hindlabel$test.subject="INDLABEL"


#Take a look at the two objects. The hybrid index is determined by the population of origin of each individual.

prepdata  #Source_allele = ygenomsamp#
hindlabel

################################################################
#Run genomic clines first on the data sampled from the #########
#genomic cline, then the data sampled from the geographic cline#
################################################################

#The results may not be identical because the two datasets were sampled independently (but using the same 
#set of individual hybrid indices).

set.seed(12345)

gc_10_ygenomsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-2,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_10_ygenomsamp

###

#Now change the 'Source_allele' column, which indicates the individual haplotypes, to represent the data 
#sampled from the geographic cline function. Run genomic cline estimation again.

prepdata$data.prep[,Source_allele:=ygeogsamp]

set.seed(23456)

gc_10_ygeogsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-2,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_10_ygeogsamp


##############################################################################################
#Now run hzar, first for the genomic cline simulated data, then the geographic cline data#####
##############################################################################################

sampN=clineres[,.N,by=x]
setkey(sampN,x);setkey(ygeogmeans,x);setkey(ygenommeans,x);
sampmeans=ygeogmeans[ygenommeans][sampN]
sampmeans_10=sampmeans

#########
#genomic#
#########
mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygenomsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

#mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1); #Hashed out because I end up with infinite confidence intervals for centre#

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_10_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_10_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelData_10_ygenomsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_10_ygenomsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_10_ygenomsamp,c("center","width")));
###

############
#geographic#
############

mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygeogsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

#mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_10_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_10_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelData_10_ygeogsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_10_ygeogsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_10_ygeogsamp,c("center","width")));
###


###########################################################
#Everything is now repeated for all the other sample sizes#
###########################################################


#############################################################################################
#############################sampsize 100####################################################
#############################################################################################

pMin=0
pMax=1
centre=0.5
x=seq(0,1,length=21)               #length = number of sampling sites
sampsize=100                     #Total samples from each site

set.seed(15432)
x2=sample(x,sampsize,replace=T)

x2=sort(x2)
rm(x)

clineres=data.table(x=x2,centre=centre,width=0.1);
clineres[,v:=1/width][,u:=qlogis(centre)*(1/width)]

set.seed(54321)
clineres[,ygeogsamp:=rbinom(sampsize,1,prob=pMin + (pMax - pMin) * (1/(1 + exp(-((x - centre) * 4/width)))))]#Sampling from a geographic cline#

set.seed(43215)
clineres[,ygenomsamp:=rbinom(sampsize,1,prob=pMin + (x^v/(x^v + (1 - x)^v*exp(u)))*(pMax - pMin))]#Sampling from a genomic cline#


ygeogmeans=clineres[,mean(ygeogsamp),by=x];setnames(ygeogmeans,"V1","ygeogsampmeans");
ygenommeans=clineres[,mean(ygenomsamp),by=x];setnames(ygenommeans,"V1","ygenomsampmeans");

plot(ygeogmeans$ygeogsampmeans~ygeogmeans$x,ylim=c(0,1),xlim=c(0,1))
points(ygenommeans$ygenomsampmeans~ygenommeans$x,col="blue",pch=3)


##################
#Preparing for GC#
##################

clineres[,INDLABEL:=paste("a",seq(1,.N),sep="")]
clineres[,POPID:=paste("a",x,sep="")]

clineres_100=clineres

prepdata=list()
prepdata$data.prep=clineres[,.(INDLABEL,POPID,ygenomsamp,ygeogsamp)]
prepdata$data.prep[,Source_allele:=ygenomsamp]
hindlabel=list()
hindlabel$hi=clineres[,.(x)]

setnames(hindlabel$hi,"x","h_posterior_mode")

prepdata$data.prep[,Source:="TEST"][,S0.prop_1:=pMin][,S1.prop_1:=pMax][,locus:="L1"]

hindlabel$hi[,INDLABEL:=prepdata$data.prep[,INDLABEL]]
hindlabel$hi[,POPID:=prepdata$data.prep[,POPID]]
hindlabel$test.subject="INDLABEL"
################
#Try running gc#
################

set.seed(32154)

gc_100_ygenomsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_100_ygenomsamp



###
prepdata$data.prep[,Source_allele:=ygeogsamp]

set.seed(21543)

gc_100_ygeogsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-2,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_100_ygeogsamp


##################
#Now run hzar#####
##################
rm(sampmeans)

sampN=clineres[,.N,by=x]
setkey(sampN,x);setkey(ygeogmeans,x);setkey(ygenommeans,x);
sampmeans=ygeogmeans[ygenommeans][sampN]
sampmeans_100=sampmeans

###
###
mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygenomsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

#mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_100_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_100_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelData_100_ygenomsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_100_ygenomsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_100_ygenomsamp,c("center","width")));
###


mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygeogsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_100_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_100_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelData_100_ygeogsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_100_ygeogsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_100_ygeogsamp,c("center","width")));
###


##############################################################################################
#############################sampsize 1000####################################################
##############################################################################################

pMin=0
pMax=1
centre=0.5
x=seq(0,1,length=21)               #length = number of sampling sites
sampsize=1000                     #Total samples from each site

set.seed(1234)
x2=sample(x,sampsize,replace=T)

x2=sort(x2)
rm(x)


clineres=data.table(x=x2,centre=centre,width=0.1);
clineres[,v:=1/width][,u:=qlogis(centre)*(1/width)]

set.seed(4123))
clineres[,ygeogsamp:=rbinom(sampsize,1,prob=pMin + (pMax - pMin) * (1/(1 + exp(-((x - centre) * 4/width)))))]#Sampling from a geographic cline#

set.seed(3412)
clineres[,ygenomsamp:=rbinom(sampsize,1,prob=pMin + (x^v/(x^v + (1 - x)^v*exp(u)))*(pMax - pMin))]#Sampling from a genomic cline#


ygeogmeans=clineres[,mean(ygeogsamp),by=x];setnames(ygeogmeans,"V1","ygeogsampmeans");
ygenommeans=clineres[,mean(ygenomsamp),by=x];setnames(ygenommeans,"V1","ygenomsampmeans");

plot(ygeogmeans$ygeogsampmeans~ygeogmeans$x,ylim=c(0,1))
points(ygenommeans$ygenomsampmeans~ygenommeans$x,col="blue")


##################
#Preparing for GC#
##################

clineres[,INDLABEL:=paste("a",seq(1,.N),sep="")]
clineres[,POPID:=paste("a",x,sep="")]

clineres_1000=clineres

prepdata=list()
prepdata$data.prep=clineres[,.(INDLABEL,POPID,ygenomsamp,ygeogsamp)]
prepdata$data.prep[,Source_allele:=ygenomsamp]
hindlabel=list()
hindlabel$hi=clineres[,.(x)]

setnames(hindlabel$hi,"x","h_posterior_mode")

prepdata$data.prep[,Source:="TEST"][,S0.prop_1:=pMin][,S1.prop_1:=pMax][,locus:="L1"]

hindlabel$hi[,INDLABEL:=prepdata$data.prep[,INDLABEL]]
hindlabel$hi[,POPID:=prepdata$data.prep[,POPID]]
hindlabel$test.subject="INDLABEL"
################
#Try running gc#
################

set.seed(2341)

gc_1000_ygenomsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_1000_ygenomsamp



###
prepdata$data.prep[,Source_allele:=ygeogsamp]

set.seed(11111)

gc_1000_ygeogsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-2,15),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_1000_ygeogsamp


##################
#Now run hzar#####
##################
rm(sampmeans)

sampN=clineres[,.N,by=x]
setkey(sampN,x);setkey(ygeogmeans,x);setkey(ygenommeans,x);
sampmeans=ygeogmeans[ygenommeans][sampN];
sampmeans_1000=sampmeans

###
###
mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygenomsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_1000_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_1000_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelData_1000_ygenomsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_1000_ygenomsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_1000_ygenomsamp,c("center","width")));
###


mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygeogsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_1000_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_1000_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelData_1000_ygeogsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_1000_ygeogsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_1000_ygeogsamp,c("center","width")));
###




###############################################################################################
#############################sampsize 10000####################################################
###############################################################################################

pMin=0
pMax=1
centre=0.5
x=seq(0,1,length=21)               #length = number of sampling sites
sampsize=10000                     #Total samples from each site

set.seed(123456)
x2=sample(x,sampsize,replace=T)

x2=sort(x2)
rm(x)


clineres=data.table(x=x2,centre=centre,width=0.1);
clineres[,v:=1/width][,u:=qlogis(centre)*(1/width)]


set.seed(612345)
clineres[,ygeogsamp:=rbinom(sampsize,1,prob=pMin + (pMax - pMin) * (1/(1 + exp(-((x - centre) * 4/width)))))]#Sampling from a geographic cline#

set.seed(561234)
clineres[,ygenomsamp:=rbinom(sampsize,1,prob=pMin + (x^v/(x^v + (1 - x)^v*exp(u)))*(pMax - pMin))]#Sampling from a genomic cline#


ygeogmeans=clineres[,mean(ygeogsamp),by=x];setnames(ygeogmeans,"V1","ygeogsampmeans");
ygenommeans=clineres[,mean(ygenomsamp),by=x];setnames(ygenommeans,"V1","ygenomsampmeans");

plot(ygeogmeans$ygeogsampmeans~ygeogmeans$x,ylim=c(0,1))
points(ygenommeans$ygenomsampmeans~ygenommeans$x,col="blue")


##################
#Preparing for GC#
##################

clineres[,INDLABEL:=paste("a",seq(1,.N),sep="")]
clineres[,POPID:=paste("a",x,sep="")]

clineres_10000=clineres

prepdata=list()
prepdata$data.prep=clineres[,.(INDLABEL,POPID,ygenomsamp,ygeogsamp)]
prepdata$data.prep[,Source_allele:=ygenomsamp]
hindlabel=list()
hindlabel$hi=clineres[,.(x)]

setnames(hindlabel$hi,"x","h_posterior_mode")

prepdata$data.prep[,Source:="TEST"][,S0.prop_1:=pMin][,S1.prop_1:=pMax][,locus:="L1"]

hindlabel$hi[,INDLABEL:=prepdata$data.prep[,INDLABEL]]
hindlabel$hi[,POPID:=prepdata$data.prep[,POPID]]
hindlabel$test.subject="INDLABEL"
################
#Try running gc#
################

set.seed(456123)

gc_10000_ygenomsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,12),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_10000_ygenomsamp



###
prepdata$data.prep[,Source_allele:=ygeogsamp]

set.seed(234561)

gc_10000_ygeogsamp=ggcline(
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel,
    read.data.precols=c("INDLABEL","POPID"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,12),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50)

gc_10000_ygeogsamp


##################
#Now run hzar#####
##################
rm(sampmeans)

sampN=clineres[,.N,by=x]
setkey(sampN,x);setkey(ygeogmeans,x);setkey(ygenommeans,x);
sampmeans=ygeogmeans[ygenommeans][sampN]
sampmeans_10000=sampmeans

###
###
mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygenomsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_10000_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_10000_ygenomsamp <- hzar.dataGroup.add(mknAdaAmodelData_10000_ygenomsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_10000_ygenomsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_10000_ygenomsamp,c("center","width")));
###


mknAdaA <-
hzar.doMolecularData1DPops(sampmeans$x,
sampmeans$ygeogsampmeans,
sampmeans$N);

mknAdaAmodel <- hzar.makeCline1DFreq(mknAdaA, tails="none");#default scaling="none", so pMin and pMax are 0 and 1, as they were in the original sampling#

mknAdaAmodel <- hzar.model.addBoxReq(mknAdaAmodel,0,1);

mknAdaAmodel

mknAdaAmodelFitR <- hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,mknAdaA,verbose=FALSE);

mknAdaAmodelFitR

mknAdaAmodelFitR$mcmcParam$chainLength <- 2e5;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e3;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData_10000_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelFit);
## 
mknAdaAmodelData_10000_ygeogsamp <- hzar.dataGroup.add(mknAdaAmodelData_10000_ygeogsamp,hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
## 
mknAdaAmodelData_10000_ygeogsamp$ML.cline$param.free
print(hzar.getLLCutParam(mknAdaAmodelData_10000_ygeogsamp,c("center","width")));
###


############################################
############################################
#End of simulations and model fits##########
############################################
############################################




################################
#Now plot results comparison####
################################

#v, 1/width results#

restable_1000=data.table(group=c("gc_ygeogsamp","gc_ygenomsamp","hzar_ygeogsamp","hzar_ygenomsamp"))
restable_1000[,sampsize:=1000]

restable_1000[group=="gc_ygeogsamp",v:=exp(gc_1000_ygeogsamp$gc$mean_log_v)]
restable_1000[group=="gc_ygenomsamp",v:=exp(gc_1000_ygenomsamp$gc$mean_log_v)]
restable_1000[group=="hzar_ygeogsamp",v:=1/(mknAdaAmodelData_1000_ygeogsamp$ML.cline$param.free[[2]])]
restable_1000[group=="hzar_ygenomsamp",v:=1/(mknAdaAmodelData_1000_ygenomsamp$ML.cline$param.free[[2]])]

restable_1000[group=="gc_ygeogsamp",v_lower_95:=gc_1000_ygeogsamp$gc$v_lower_95]
restable_1000[group=="gc_ygenomsamp",v_lower_95:=gc_1000_ygenomsamp$gc$v_lower_95]
restable_1000[group=="hzar_ygeogsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_1000_ygeogsamp,c("center","width"))[[4]])]
restable_1000[group=="hzar_ygenomsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_1000_ygenomsamp,c("center","width"))[[4]])]

restable_1000[group=="gc_ygeogsamp",v_upper_95:=gc_1000_ygeogsamp$gc$v_upper_95]
restable_1000[group=="gc_ygenomsamp",v_upper_95:=gc_1000_ygenomsamp$gc$v_upper_95]
restable_1000[group=="hzar_ygeogsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_1000_ygeogsamp,c("center","width"))[[3]])]
restable_1000[group=="hzar_ygenomsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_1000_ygenomsamp,c("center","width"))[[3]])]

###
restable_10=data.table(group=c("gc_ygeogsamp","gc_ygenomsamp","hzar_ygeogsamp","hzar_ygenomsamp"))
restable_10[,sampsize:=10]

restable_10[group=="gc_ygeogsamp",v:=exp(gc_10_ygeogsamp$gc$mean_log_v)]
restable_10[group=="gc_ygenomsamp",v:=exp(gc_10_ygenomsamp$gc$mean_log_v)]
restable_10[group=="hzar_ygeogsamp",v:=1/(mknAdaAmodelData_10_ygeogsamp$ML.cline$param.free[[2]])]
restable_10[group=="hzar_ygenomsamp",v:=1/(mknAdaAmodelData_10_ygenomsamp$ML.cline$param.free[[2]])]

restable_10[group=="gc_ygeogsamp",v_lower_95:=gc_10_ygeogsamp$gc$v_lower_95]
restable_10[group=="gc_ygenomsamp",v_lower_95:=gc_10_ygenomsamp$gc$v_lower_95]
restable_10[group=="hzar_ygeogsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10_ygeogsamp,c("center","width"))[[4]])]
restable_10[group=="hzar_ygenomsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10_ygenomsamp,c("center","width"))[[4]])]

restable_10[group=="gc_ygeogsamp",v_upper_95:=gc_10_ygeogsamp$gc$v_upper_95]
restable_10[group=="gc_ygenomsamp",v_upper_95:=gc_10_ygenomsamp$gc$v_upper_95]
restable_10[group=="hzar_ygeogsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10_ygeogsamp,c("center","width"))[[3]])]
restable_10[group=="hzar_ygenomsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10_ygenomsamp,c("center","width"))[[3]])]

###
restable_100=data.table(group=c("gc_ygeogsamp","gc_ygenomsamp","hzar_ygeogsamp","hzar_ygenomsamp"))
restable_100[,sampsize:=100]

restable_100[group=="gc_ygeogsamp",v:=exp(gc_100_ygeogsamp$gc$mean_log_v)]
restable_100[group=="gc_ygenomsamp",v:=exp(gc_100_ygenomsamp$gc$mean_log_v)]
restable_100[group=="hzar_ygeogsamp",v:=1/(mknAdaAmodelData_100_ygeogsamp$ML.cline$param.free[[2]])]
restable_100[group=="hzar_ygenomsamp",v:=1/(mknAdaAmodelData_100_ygenomsamp$ML.cline$param.free[[2]])]

restable_100[group=="gc_ygeogsamp",v_lower_95:=gc_100_ygeogsamp$gc$v_lower_95]
restable_100[group=="gc_ygenomsamp",v_lower_95:=gc_100_ygenomsamp$gc$v_lower_95]
restable_100[group=="hzar_ygeogsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_100_ygeogsamp,c("center","width"))[[4]])]
restable_100[group=="hzar_ygenomsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_100_ygenomsamp,c("center","width"))[[4]])]

restable_100[group=="gc_ygeogsamp",v_upper_95:=gc_100_ygeogsamp$gc$v_upper_95]
restable_100[group=="gc_ygenomsamp",v_upper_95:=gc_100_ygenomsamp$gc$v_upper_95]
restable_100[group=="hzar_ygeogsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_100_ygeogsamp,c("center","width"))[[3]])]
restable_100[group=="hzar_ygenomsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_100_ygenomsamp,c("center","width"))[[3]])]


###
restable_10000=data.table(group=c("gc_ygeogsamp","gc_ygenomsamp","hzar_ygeogsamp","hzar_ygenomsamp"))
restable_10000[,sampsize:=10000]

restable_10000[group=="gc_ygeogsamp",v:=exp(gc_10000_ygeogsamp$gc$mean_log_v)]
restable_10000[group=="gc_ygenomsamp",v:=exp(gc_10000_ygenomsamp$gc$mean_log_v)]
restable_10000[group=="hzar_ygeogsamp",v:=1/(mknAdaAmodelData_10000_ygeogsamp$ML.cline$param.free[[2]])]
restable_10000[group=="hzar_ygenomsamp",v:=1/(mknAdaAmodelData_10000_ygenomsamp$ML.cline$param.free[[2]])]

restable_10000[group=="gc_ygeogsamp",v_lower_95:=gc_10000_ygeogsamp$gc$v_lower_95]
restable_10000[group=="gc_ygenomsamp",v_lower_95:=gc_10000_ygenomsamp$gc$v_lower_95]
restable_10000[group=="hzar_ygeogsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10000_ygeogsamp,c("center","width"))[[4]])]
restable_10000[group=="hzar_ygenomsamp",v_lower_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10000_ygenomsamp,c("center","width"))[[4]])]

restable_10000[group=="gc_ygeogsamp",v_upper_95:=gc_10000_ygeogsamp$gc$v_upper_95]
restable_10000[group=="gc_ygenomsamp",v_upper_95:=gc_10000_ygenomsamp$gc$v_upper_95]
restable_10000[group=="hzar_ygeogsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10000_ygeogsamp,c("center","width"))[[3]])]
restable_10000[group=="hzar_ygenomsamp",v_upper_95:=1/(hzar.getLLCutParam(mknAdaAmodelData_10000_ygenomsamp,c("center","width"))[[3]])]


#center/centre results#


restable_1000[group=="gc_ygeogsamp",centre:=plogis(gc_1000_ygeogsamp$gc$mean_logit_centre)]
restable_1000[group=="gc_ygenomsamp",centre:=plogis(gc_1000_ygenomsamp$gc$mean_logit_centre)]
restable_1000[group=="hzar_ygeogsamp",centre:=(mknAdaAmodelData_1000_ygeogsamp$ML.cline$param.free[[1]])]
restable_1000[group=="hzar_ygenomsamp",centre:=(mknAdaAmodelData_1000_ygenomsamp$ML.cline$param.free[[1]])]

restable_1000[group=="gc_ygeogsamp",centre_lower_95:=gc_1000_ygeogsamp$gc$centre_lower_95]
restable_1000[group=="gc_ygenomsamp",centre_lower_95:=gc_1000_ygenomsamp$gc$centre_lower_95]
restable_1000[group=="hzar_ygeogsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_1000_ygeogsamp,c("center","width"))[[1]])]
restable_1000[group=="hzar_ygenomsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_1000_ygenomsamp,c("center","width"))[[1]])]

restable_1000[group=="gc_ygeogsamp",centre_upper_95:=gc_1000_ygeogsamp$gc$centre_upper_95]
restable_1000[group=="gc_ygenomsamp",centre_upper_95:=gc_1000_ygenomsamp$gc$centre_upper_95]
restable_1000[group=="hzar_ygeogsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_1000_ygeogsamp,c("center","width"))[[2]])]
restable_1000[group=="hzar_ygenomsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_1000_ygenomsamp,c("center","width"))[[2]])]

###

restable_10[group=="gc_ygeogsamp",centre:=plogis(gc_10_ygeogsamp$gc$mean_logit_centre)]
restable_10[group=="gc_ygenomsamp",centre:=plogis(gc_10_ygenomsamp$gc$mean_logit_centre)]
restable_10[group=="hzar_ygeogsamp",centre:=(mknAdaAmodelData_10_ygeogsamp$ML.cline$param.free[[1]])]
restable_10[group=="hzar_ygenomsamp",centre:=(mknAdaAmodelData_10_ygenomsamp$ML.cline$param.free[[1]])]

restable_10[group=="gc_ygeogsamp",centre_lower_95:=gc_10_ygeogsamp$gc$centre_lower_95]
restable_10[group=="gc_ygenomsamp",centre_lower_95:=gc_10_ygenomsamp$gc$centre_lower_95]
restable_10[group=="hzar_ygeogsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_10_ygeogsamp,c("center","width"))[[1]])]
restable_10[group=="hzar_ygenomsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_10_ygenomsamp,c("center","width"))[[1]])]

restable_10[group=="gc_ygeogsamp",centre_upper_95:=gc_10_ygeogsamp$gc$centre_upper_95]
restable_10[group=="gc_ygenomsamp",centre_upper_95:=gc_10_ygenomsamp$gc$centre_upper_95]
restable_10[group=="hzar_ygeogsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_10_ygeogsamp,c("center","width"))[[2]])]
restable_10[group=="hzar_ygenomsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_10_ygenomsamp,c("center","width"))[[2]])]

###

restable_100[group=="gc_ygeogsamp",centre:=plogis(gc_100_ygeogsamp$gc$mean_logit_centre)]
restable_100[group=="gc_ygenomsamp",centre:=plogis(gc_100_ygenomsamp$gc$mean_logit_centre)]
restable_100[group=="hzar_ygeogsamp",centre:=(mknAdaAmodelData_100_ygeogsamp$ML.cline$param.free[[1]])]
restable_100[group=="hzar_ygenomsamp",centre:=(mknAdaAmodelData_100_ygenomsamp$ML.cline$param.free[[1]])]

restable_100[group=="gc_ygeogsamp",centre_lower_95:=gc_100_ygeogsamp$gc$centre_lower_95]
restable_100[group=="gc_ygenomsamp",centre_lower_95:=gc_100_ygenomsamp$gc$centre_lower_95]
restable_100[group=="hzar_ygeogsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_100_ygeogsamp,c("center","width"))[[1]])]
restable_100[group=="hzar_ygenomsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_100_ygenomsamp,c("center","width"))[[1]])]

restable_100[group=="gc_ygeogsamp",centre_upper_95:=gc_100_ygeogsamp$gc$centre_upper_95]
restable_100[group=="gc_ygenomsamp",centre_upper_95:=gc_100_ygenomsamp$gc$centre_upper_95]
restable_100[group=="hzar_ygeogsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_100_ygeogsamp,c("center","width"))[[2]])]
restable_100[group=="hzar_ygenomsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_100_ygenomsamp,c("center","width"))[[2]])]


###

restable_10000[group=="gc_ygeogsamp",centre:=plogis(gc_10000_ygeogsamp$gc$mean_logit_centre)]
restable_10000[group=="gc_ygenomsamp",centre:=plogis(gc_10000_ygenomsamp$gc$mean_logit_centre)]
restable_10000[group=="hzar_ygeogsamp",centre:=(mknAdaAmodelData_10000_ygeogsamp$ML.cline$param.free[[1]])]
restable_10000[group=="hzar_ygenomsamp",centre:=(mknAdaAmodelData_10000_ygenomsamp$ML.cline$param.free[[1]])]

restable_10000[group=="gc_ygeogsamp",centre_lower_95:=gc_10000_ygeogsamp$gc$centre_lower_95]
restable_10000[group=="gc_ygenomsamp",centre_lower_95:=gc_10000_ygenomsamp$gc$centre_lower_95]
restable_10000[group=="hzar_ygeogsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_10000_ygeogsamp,c("center","width"))[[1]])]
restable_10000[group=="hzar_ygenomsamp",centre_lower_95:=(hzar.getLLCutParam(mknAdaAmodelData_10000_ygenomsamp,c("center","width"))[[1]])]

restable_10000[group=="gc_ygeogsamp",centre_upper_95:=gc_10000_ygeogsamp$gc$centre_upper_95]
restable_10000[group=="gc_ygenomsamp",centre_upper_95:=gc_10000_ygenomsamp$gc$centre_upper_95]
restable_10000[group=="hzar_ygeogsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_10000_ygeogsamp,c("center","width"))[[2]])]
restable_10000[group=="hzar_ygenomsamp",centre_upper_95:=(hzar.getLLCutParam(mknAdaAmodelData_10000_ygenomsamp,c("center","width"))[[2]])]


##################################################

###

restable=rbind(restable_10,restable_100,restable_1000,restable_10000)

restable[,logsampsize:=log(sampsize)]

restable[group=="gc_ygenomsamp",group2:=as.factor(1)]
restable[group=="hzar_ygenomsamp",group2:=as.factor(2)]
restable[group=="gc_ygeogsamp",group2:=as.factor(3)]
restable[group=="hzar_ygeogsamp",group2:=as.factor(4)]




####################################################################################
#Closed circles are estimates from gghybrid; closed squares are estimates from hzar#
####################################################################################
setkey(restable,logsampsize,group2)
restable[,pch:=rep(c(16,15),8)]


########################################################################################################################################################
#The colour represents the cline function from which the haplotypes were sampled (pinkish = genomic cline function; bluish = geographic cline function)#
########################################################################################################################################################
restable[group%like%"genom",colfac:="genom"]
restable[group%like%"geog",colfac:="geog"]

restable[,color:=factor(colfac,levels=c("genom","geog"),labels=rainbow_hcl(2))]


 geom_line(data=clineres2[centre==0.5],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.1],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2)



##########################################################################################
##########################################################################################
#Combine with 3 plots of logit-logistic curves only#######################################
##########################################################################################
##########################################################################################


#The cline function, with null defaults for v and centre:

eq = function(x,v=1,centre=0.5){u=qlogis(centre)*v;(x^v/(x^v + (1 - x)^v*exp(u)))}


###############################
#v remaining at the null value#
###############################

#No s-shaped curve when the slope is not steeper than null.

centres=seq(0.000001,0.999999,length=9)




##################################################################################
##################################################################################
#All 6 plots combined for Figure 1################################################
##################################################################################
##################################################################################



base1 <- ggplot() + xlim(0, 1)
base1 = base1 + geom_function(fun = eq, args = list(v=1,centre=centres[5]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[1]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[2]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[3]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[4]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[6]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[7]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[8]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=1,centre=centres[9]),colour="#E495A5",lwd=1.2) +
geom_vline(xintercept=centres,col=c("grey","grey","grey","grey","black","grey","grey","grey","grey"),lty=2) +
geom_point(aes(x=centres,y=rep(0.5,9))) +
theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="genome-wide hybrid index", y = "locus-specific allele frequency", tag = "B")

###############################
#Steep shifted clines##########
###############################

base2 <- ggplot() + xlim(0, 1)
base2 = base2 + geom_function(fun = eq, args = list(v=30,centre=centres[5]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[1]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[2]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[3]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[4]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[6]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[7]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[8]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=30,centre=centres[9]),colour="#E495A5",lwd=1.2) +
geom_vline(xintercept=centres,col=c("grey","grey","grey","grey","black","grey","grey","grey","grey"),lty=2) +
geom_point(aes(x=centres,y=rep(0.5,9))) +
theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="genome-wide hybrid index", y = "locus-specific allele frequency", tag = "C")

###############################
#Shallow shifted clines########
###############################

#Same slope for all these lines.

base3 <- ggplot() + xlim(0, 1)
base3 = base3 + geom_function(fun = eq, args = list(v=0.5,centre=centres[5]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[1]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[2]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[3]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[4]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[6]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[7]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[8]),colour="#E495A5",lwd=1.2) +
geom_function(fun = eq, args = list(v=0.5,centre=centres[9]),colour="#E495A5",lwd=1.2) +
geom_vline(xintercept=centres,col=c("grey","grey","grey","grey","black","grey","grey","grey","grey"),lty=2) +
geom_point(aes(x=centres,y=rep(0.5,9))) +
theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="genome-wide hybrid index", y = "locus-specific allele frequency", tag = "A")

###
p1= ggplot() + xlim(0, 1) + ylim(0,1)
p1 = p1 + geom_line(data=clineres2[centre==0.5],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.5],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.1],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.1],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.2],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.2],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.3],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.3],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.4],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.4],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.5],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.5],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.6],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.6],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.7],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.7],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.8],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.8],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) + 
 geom_line(data=clineres2[centre==0.9],aes(x=x,y=ygenom),col="#E495A5",lwd=1.2) +
 geom_line(data=clineres2[centre==0.9],aes(x=x,y=ygeog),col="#39BEB1",lty=2,lwd=1.2) +
 geom_vline(xintercept=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),col=c("grey","grey","grey","grey","black","grey","grey","grey","grey"),lty=2) +
geom_point(aes(x=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),y=rep(0.5,9))) +
 theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="genome-wide hybrid index", y = "locus-specific allele frequency", tag = "D")

###
p2=ggplot(data=restable,aes(x=logsampsize,y=v, group=group2, color=colfac),position="jitter")
p2 = p2 + geom_point(aes(x=logsampsize,y=v), position=position_dodge(width=0.5), pch=restable$pch) +
 geom_errorbar(aes(ymin=v_lower_95, ymax=v_upper_95),width=0, position=position_dodge(width=0.5)) + 
 ylim(0,125) +
 geom_hline(yintercept=10) +
#scale_color_discrete_qualitative(palette = "Set 2") +
 #scale_color_manual(values=restable$color) +#c("blue","blue","orange","orange")) +
 theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="log(sample size)", y = "v (gghybrid) or 1/width (hzar)", tag = "E")

###
p3=ggplot(data=restable,aes(x=logsampsize,y=centre, group=as.factor(group2), color=colfac),position="jitter")
p3 = p3 + geom_point(aes(x=logsampsize,y=centre), position=position_dodge(width=0.5), pch=restable$pch) +
 geom_errorbar(aes(ymin=centre_lower_95, ymax=centre_upper_95),width=0, position=position_dodge(width=0.5)) + 
 ylim(0,1) +
 geom_hline(yintercept=0.5) +
scale_color_discrete_qualitative(palette = "Set 2") +
 #scale_color_manual(values=rainbow_hcl(2)) +#c("blue","blue","orange","orange")) +
 theme_classic() + 
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="log(sample size)", y = "cline centre", tag = "F")


##############################################################################

p=plot_grid(base3,base1,base2,p1, p2, p3,ncol=3)

p

ggsave("logitlogistic_gghybrid_hzar_comparison_plot.pdf",p)

#"Saving 12 x 7.36 in image"









