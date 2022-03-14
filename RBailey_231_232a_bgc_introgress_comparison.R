##############################################################
#Richard Ian Bailey 11 March 2022#############################
##############################################################

#Numbered the same as in the paper.

##############################################################################
##############################################################################
#2.3.1 Comparison of gghybrid hybrid index estimation with introgress and bgc#
##############################################################################
##############################################################################

###################################################################
#Comparisons using Hermansen et al Italy analysis dataset##########
###################################################################


#Run gghybrid and INTROGRESS for hybrid index estimation and compare, alongside existing bgc hybrid index estimates.

###################################################################
###################################################################
###################################################################
library(gghybrid)
library(introgress)#Now obsolete!
library(ggplot2)
library(cowplot)
###################################################################
#Make sure I'm using the same individuals and loci#################
###################################################################

#Table S1 Hermansen et al contains SNP names and their parental allele frequencies.

#Figure S3 Hermansen et al shows clines for individual loci - states there were 77 loci.

#Loading files downloaded from Dryad for Hermansen et al 2014.


#Set the working directory to the one containing the bgc files from Hermansen et al.


#alpha estimates per locus italy#
aout=fread("a.out_italy.txt",header=T)#No locus IDs#

#beta estimates per locus italy#
bout=fread("b.out_italy.txt",header=T)#No locus IDs#

#hybrid index estimates italy#
hiout=fread("hi.out_italy.txt",header=T)#433 individuals, no IDs#

#These below are cline parameter estimates on the data scale#
qaout=fread("qa.out_italy.txt",header=T)
qbout=fread("qb.out_italy.txt",header=T)

################################
#Cross reference the locus data#
################################
aout[,locusnum:=tstrsplit(param,"_",keep=3L)]
bout[,locusnum:=tstrsplit(param,"_",keep=3L)]
qaout[,locusnum:=tstrsplit(param,"_",keep=3L)]
qbout[,locusnum:=tstrsplit(param,"_",keep=3L)]

setnames(aout,c("mean","median","0.950_CI_LB","0.950_CI_UB"),c("a_mean","a_median","a_0.950_CI_LB","a_0.950_CI_UB"))
setnames(bout,c("mean","median","0.950_CI_LB","0.950_CI_UB"),c("b_mean","b_median","b_0.950_CI_LB","b_0.950_CI_UB"))
setnames(qaout,c("mean","median","0.950_CI_LB","0.950_CI_UB"),c("qa_mean","qa_median","qa_0.950_CI_LB","qa_0.950_CI_UB"))
setnames(qbout,c("mean","median","0.950_CI_LB","0.950_CI_UB"),c("qb_mean","qb_median","qb_0.950_CI_LB","qb_0.950_CI_UB"))

###
loci=fread("loci.csv")
loci[,locusnum:=tstrsplit(param,"_",keep=3L)]

#Joining all the loaded files.
setkey(loci,locusnum);setkey(aout,locusnum);setkey(bout,locusnum);setkey(qaout,locusnum);setkey(qbout,locusnum)
loci=loci[aout][bout][qaout][qbout]

loci[,locusnum:=as.numeric(locusnum)]

loci[order(locusnum)]

####################################################################
#Now cross reference the hybrid index estimates with individual ID##
####################################################################

spanishid=fread("id_spanish.csv",header=T)#52 individuals#
spanishid

houseid=fread("id_house.csv",header=T)#85 individuals#
houseid


#bgc only uses parental allele frequencies, not their individual IDs, so the IDs in the bgc hybrid index file should correspond to the italy file below.

italyid=fread("id_italy.csv",header=T)
italyid

str(italyid)
str(hiout)

hiout[,BGC:=tstrsplit(param,"_",keep=3L)]
hiout[,BGC:=as.numeric(BGC)]

setkey(italyid,BGC);setkey(hiout,BGC)
italyhi=italyid[hiout]

#italyhi now contains the hybrid index estimates for all test individuals (i.e. not parental reference) in the Italy dataset.

#############################################################################
#Estimate hybrid index using gghybrid and introgress for comparison with bgc#
#############################################################################

#Load the dataset I prepared elsewhere from the Hermansen et al dataset.

####################################################################################
#Using dataset with Mt locus ND2 coded as diploid homozygote, to match bgc analysis#
####################################################################################

################
#gghybrid first#
################
dat=read.data("RB_Italy_ggcline_precol_headers_diploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569)

dat#77 loci & 569 male individuals, including parental reference#


prepdata=data.prep(data=dat$data,
 loci=dat$loci,
 alleles=dat$alleles,
 S0=c("Kralove","Oslo"), #POPID names for the first parental reference set#
 S1=c("LesinaSPANISH","Lauro_Lesina","MasseriaFontanarosa","MasseriaMustazzo"),     #POPID names for the second parental reference set#
 precols=dat$precols,
 return.genotype.table=T,
 return.locus.table=T)

prepdata

#Some loci have fewer Spanish parental reference samples than others because they weren't gentoyped in some populations.

##################################################
#Check the number of test individuals is correct##
##################################################

prepdata$data.prep[Source=="TEST",uniqueN(INDLABEL)]#433, same as italyhi#

prepdata$data.prep[Source=="S0",uniqueN(INDLABEL)]#84 - Hermansen et al had 85. They have one duplicated individual!!#
prepdata$data.prep[Source=="S1",uniqueN(INDLABEL)]#52

####################################################
#Check N parental copies of each allele corresponds#
####################################################

prepdata$locus.data[,S0.N_0:=N_allele_copies.S0 - S0.N_1]
prepdata$locus.data[,S1.N_0:=N_allele_copies.S1 - S1.N_1]

loci[,.(chrom,locus,locusnum)]
prepdata$locus.data[,.(locus,S0.N_0,S0.N_1,S1.N_0,S1.N_1)]

setkey(loci,locus);setkey(prepdata$locus.data,locus)
allcomp=loci[,.(chrom,locus,locusnum)][prepdata$locus.data[,.(locus,S0.N_0,S0.N_1,S1.N_0,S1.N_1)]]

setkey(allcomp,locusnum)

#Now all the per locus info is in prepdata$locus.data

#############################################################
#Compare locus.data to the parental reference files from bgc#
#############################################################

#Spanish parental allele frequency estimated using information from 8 individuals for the markers listed below. 
#Spanish sparrow allele frequencies for rest of the markers estimated using information from 52 individuals.

#PIK3C3,TNPO1,RSF1,RABGAP1l,UnTrC8

#Again using files downloaded from Dryad for Hermansen et al. 2014.

spanish=fread("spanish.txt")
spanish=spanish[rep(c(FALSE,TRUE),length=.N)]
spanish[,locusnum:=seq(0,76)]

house=fread("house.txt")
house=house[rep(c(FALSE,TRUE),length=.N)]
house[,locusnum:=seq(0,76)]

#ND2 (locusnum=0 in the above files) appears to be coded as haploid in the house parental file (house.txt). 
#The locus is fixed so this won't affect the parental allele frequency#

#Otherwise everything looks fine.


##########################
#Now run esth#############
##########################

#########################################################################################
#For comparison with bgc results from Hermansen et al and new introgress analysis, below#
#########################################################################################

set.seed(98765)
hindlabel=esth(
 data.prep.object=prepdata$data.prep,
 read.data.precols=dat$precols,
 include.Source=TRUE,
 plot.ind = c("PD07-254","PD07-160","PD07-159","PI07-110","PI08-498","PH08-285"),
 plot.col = c("blue","green","cyan","purple","magenta","red"),
 nitt=3000,
 burnin=1000
)

View(hindlabel$hi[order(beta_mean)])#Ordering by beta_mean because it's never exactly 0 or 1#

###############################################
#Now run hybrid index estimation in introgress#
###############################################

#Preparing the input files. prepdata$geno.data was created above by gghybrid's data.prep function.

geno1.loci=data.table(locus=names(prepdata$geno.data[,!c("INDLABEL","POPID","Source")]))
geno1.loci[,type:="C"]

#Save to the working directory.
fwrite(geno1.loci,"geno1_loci.txt",na=NA)

#Converting the allele ids from numeric to character.
geno1=copy(prepdata$geno.data)
geno1[]=lapply(geno1,function(x) as.character(x))

str(geno1)

#Replacing entries with the format required by introgress.

locfix0=function(x){
 replace(x,x=="0","0/0")
}
locfix1=function(x){
 replace(x,x=="1","1/0")
}
locfix2=function(x){
 replace(x,x=="2","1/1")
}

geno1=geno1[,lapply(.SD,locfix0)]
geno1=geno1[,lapply(.SD,locfix1)]
geno1=geno1[,lapply(.SD,locfix2)]

#placing parental reference and test individuals in separate, correctly formatted files.
geno1.admix=geno1[Source=="TEST"]
geno1.S0=geno1[Source=="S0"]
geno1.S1=geno1[Source=="S1"]

geno1.S0[,INDLABEL:=NULL][,POPID:=NULL][,Source:=NULL]
geno1.S1[,INDLABEL:=NULL][,POPID:=NULL][,Source:=NULL]

geno1.admix[,Source:=NULL]
dim(geno1.admix)
setcolorder(geno1.admix,c(2,1,3:79))

geno1.admix=t(geno1.admix)
fwrite(geno1.admix,"geno1_admix.txt",na=NA)#NEED TO OPEN THE FILE AND REMOVE THE HEADER ROW BEFORE RUNNING THE CODE BELOW#

geno1.S0=t(geno1.S0)
fwrite(geno1.S0,"geno1_S0.txt",na=NA)#NEED TO OPEN THE FILE AND REMOVE THE HEADER ROW BEFORE RUNNING THE CODE BELOW#

geno1.S1=t(geno1.S1)
fwrite(geno1.S1,"geno1_S1.txt",na=NA)#NEED TO OPEN THE FILE AND REMOVE THE HEADER ROW BEFORE RUNNING THE CODE BELOW#

#######################################
#Run introgress########################
#######################################

## read in data for individuals from admixed population
gen1.data.sim <- read.csv(file="geno1_admix.txt", header=FALSE)

## read in marker information
loci1.data.sim <- read.csv(file="geno1_loci.txt", header=TRUE)

## read in parental data sets
p1.data1 <- read.csv(file="geno1_S0.txt", header=FALSE)
p2.data1 <- read.csv(file="geno1_S1.txt", header=FALSE)

## code to convert genotype data into a matrix of allele counts,
## the results are saved to the list data object count.matrix

count.matrix1 <- prepare.data(admix.gen=gen1.data.sim, loci.data=loci1.data.sim,
                             parental1=p1.data1, parental2=p2.data1)


## estimate hybrid index values and save the results to the
## data.frame data object hi.index.sim1

set.seed(76543)
hi1.index.sim <- est.h(introgress.data=count.matrix1, loci.data=loci1.data.sim)


######################################################
#Comparison between gghybrid and introgress###########
######################################################

inds=data.table(t(gen1.data.sim[1:2,]))
setnames(inds,c("1","2"),c("POPID","INDLABEL"))

hi1.index.sim=data.table(hi1.index.sim)

hi1.index.sim=cbind(inds,hi1.index.sim)#The introgress results#

hindlabel#The gghybrid results#

setkey(hindlabel$hi,INDLABEL);setkey(hi1.index.sim,INDLABEL)
introres=hindlabel$hi[hi1.index.sim]

par(mfrow=c(2,2))
plot(introres[,qlogis(beta_mean)]~introres[,qlogis(h)])#this loses zeroes#
abline(a=0,b=1,col="red",lty=2)
plot(introres[,qlogis(h_posterior_mode)]~introres[,qlogis(h)],col="blue",pch=3)#this loses zeroes#
abline(a=0,b=1,col="red",lty=2)
plot(introres[,beta_mean]~introres[,h])
abline(a=0,b=1,col="red",lty=2)
plot(introres[,h_posterior_mode]~introres[,h],col="blue",pch=3)
abline(a=0,b=1,col="red",lty=2)


###
cor.test(introres$h,introres$h_posterior_mode)

        Pearson's product-moment correlation

data:  introres$h and introres$h_posterior_mode
t = 978.32, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9997281 0.9998137
sample estimates:
      cor 
0.9997749 
###
cor.test(introres$h,introres$beta_mean)

        Pearson's product-moment correlation

data:  introres$h and introres$beta_mean
t = 1316.9, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9998499 0.9998972
sample estimates:
      cor 
0.9998758 
###




################################
#Comparison with bgc############
################################

italyhi#bgc results#

hindlabel$hi#gghybrid results#

setkey(hindlabel$hi,INDLABEL);setkey(italyhi,individual)
introres_bgc=hindlabel$hi[italyhi]

introres_bgc


par(mfrow=c(2,2))
plot(introres_bgc[,qlogis(beta_mean)]~introres_bgc[,qlogis(mean)])#this loses zeroes#
abline(a=0,b=1,col="red",lty=2)
plot(introres_bgc[,qlogis(h_posterior_mode)]~introres_bgc[,qlogis(mean)],col="blue",pch=3)#this loses zeroes#
abline(a=0,b=1,col="red",lty=2)
plot(introres_bgc[,beta_mean]~introres_bgc[,mean])
abline(a=0,b=1,col="red",lty=2)
plot(introres_bgc[,h_posterior_mode]~introres_bgc[,mean],col="blue",pch=3)
abline(a=0,b=1,col="red",lty=2)


###
cor.test(introres_bgc$mean,introres_bgc$h_posterior_mode)

        Pearson's product-moment correlation

data:  introres_bgc$mean and introres_bgc$h_posterior_mode
t = 156.33, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9894955 0.9927906
sample estimates:
      cor 
0.9912969 

###
cor.test(introres_bgc$mean,introres_bgc$beta_mean)

        Pearson's product-moment correlation

data:  introres_bgc$mean and introres_bgc$beta_mean
t = 159.15, df = 431, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9898592 0.9930406
sample estimates:
      cor 
0.9915986 
###


##########################################################
##########################################################
#2.3.2 (a) comparison of gghybrid and bgc cline estimates#
##########################################################
##########################################################

#I don't know whether the bgc run included parental reference individuals, so run both with and without#

#And see if it's in the methods of Hermansen et al!

#I think not, because only the parental allele frequencies are given in bgc, the parental samples aren't included.

set.seed(121212)
gc_bgc_withSource=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  include.Source = TRUE,                  #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  plot.test.subject = c("CHD1Z_1","HSDL2_1"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre = c(0,sqrt(50)),     #old Default#
  prior.logv = c(0,sqrt(10)),             #old Default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000
)

set.seed(3121212)
gc_bgc_withoutSource=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  include.Source = FALSE,                  #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  plot.test.subject = c("CHD1Z_1","HSDL2_1"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre = c(0,sqrt(50)),     #old Default#
  prior.logv = c(0,sqrt(10)),             #old Default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000
)


###

#####################
#noparent comparison#
#####################

setkey(loci,locus);setkey(gc_bgc_withoutSource$gc,locus)
loci_noparent=loci[gc_bgc_withoutSource$gc]

loci_noparent

#The parameters are not directly comparable so no surprise that they don't always correspond.

#
par(mfrow=c(1,2))
#
plot(loci_noparent[,b_median]~loci_noparent[,exp_mean_log_v])
#
plot(loci_noparent[,-a_median]~loci_noparent[,invlogit_mean_logit_centre])
#

#plot(loci_noparent[,qb_median]~loci_noparent[,mean_log_v])
#
#plot(loci_noparent[,-qa_median]~loci_noparent[,mean_logit_centre])
#

###
cor.test(loci_noparent$b_median,loci_noparent$exp_mean_log_v)#Not presented in the paper#

        Pearson's product-moment correlation

data:  loci_noparent$b_median and loci_noparent$exp_mean_log_v
t = 9.4276, df = 75, p-value = 2.298e-14
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6136888 0.8244323
sample estimates:
      cor 
0.7364413 
###

cor.test(loci_noparent$a_median,loci_noparent$invlogit_mean_logit_centre)

        Pearson's product-moment correlation

data:  loci_noparent$a_median and loci_noparent$invlogit_mean_logit_centre
t = -12.66, df = 75, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.8856164 -0.7377590
sample estimates:
       cor 
-0.8253543 
###


###################
#parent comparison#
###################

setkey(loci,locus);setkey(gc_bgc_withSource$gc,locus)
loci_parent=loci[gc_bgc_withSource$gc]

loci_parent

#
par(mfrow=c(2,4))
#
plot(loci_parent[,b_median]~loci_parent[,mean_log_v])
#
plot(loci_parent[,b_median]~loci_parent[,exp_mean_log_v])
#
plot(loci_parent[,-a_median]~loci_parent[,mean_logit_centre])
#
plot(loci_parent[,-a_median]~loci_parent[,invlogit_mean_logit_centre])
#

plot(loci_parent[,qb_median]~loci_parent[,mean_log_v])
#
plot(loci_parent[,qb_median]~loci_parent[,exp_mean_log_v])
#
plot(loci_parent[,-qa_median]~loci_parent[,mean_logit_centre])
#
plot(loci_parent[,-qa_median]~loci_parent[,invlogit_mean_logit_centre])
#


###
cor.test(loci_parent$b_median,loci_parent$exp_mean_log_v)

        Pearson's product-moment correlation

data:  loci_parent$b_median and loci_parent$exp_mean_log_v
t = 11.512, df = 75, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.7005489 0.8677901
sample estimates:
      cor 
0.7991367 
###

cor.test(loci_parent$a_median,loci_parent$invlogit_mean_logit_centre)

        Pearson's product-moment correlation

data:  loci_parent$a_median and loci_parent$invlogit_mean_logit_centre
t = -12.911, df = 75, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.8890814 -0.7450976
sample estimates:
       cor 
-0.8304808 
###


################################################
#Now compare identification of significant loci#
################################################

#Significantly steep#

loci_noparent[v_pvalue<0.05&exp_mean_log_v>1]#22 loci#
loci_noparent[v_lower_95>1]#22 loci#

loci_noparent[b_0.950_CI_LB>0]#13 loci#

##
loci_parent[v_pvalue<0.05&exp_mean_log_v>1]#22 loci#
loci_parent[v_lower_95>1]#22 loci#

loci_parent[b_0.950_CI_LB>0]#13 loci#


#Significant centre#

loci_noparent[centre_pvalue<0.05]#52 loci#

loci_noparent[a_0.950_CI_LB>0&a_0.950_CI_UB>0|a_0.950_CI_LB<0&a_0.950_CI_UB<0]#55 loci#

##
loci_parent[centre_pvalue<0.05]#52 loci#

loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0|a_0.950_CI_LB<0&a_0.950_CI_UB<0]#55 loci#


#############
#Calculate u#
#############

loci_parent[,u_mean:=mean_logit_centre*exp_mean_log_v]



############################################################
##################Figure S1 below###########################
############################################################


###
p3 = ggplot()
p3 = p3 + geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95<=1],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),pch=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95>1],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="red",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95>1],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="blue",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95<=1],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="brown",pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB>0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB>0],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="red",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB<=0],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="brown",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB>0],aes(x=invlogit_mean_logit_centre,y=exp_mean_log_v),col="blue",pch=3) +
geom_vline(xintercept=0.5,lty=2,col="grey") +
geom_hline(yintercept=1,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="gghybrid cline centre", y = "gghybrid cline steepness v", tag = "A")

p3
###


#Same but plotting u instead of centre


###
p33 = ggplot()
p33 = p33 + geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95<=1],aes(x=u_mean,y=exp_mean_log_v),pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0],aes(x=u_mean,y=exp_mean_log_v),pch=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95>1],aes(x=u_mean,y=exp_mean_log_v),col="red",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95>1],aes(x=u_mean,y=exp_mean_log_v),col="blue",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95<=1],aes(x=u_mean,y=exp_mean_log_v),col="brown",pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB>0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB>0],aes(x=u_mean,y=exp_mean_log_v),col="red",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB<=0],aes(x=u_mean,y=exp_mean_log_v),col="brown",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB>0],aes(x=u_mean,y=exp_mean_log_v),col="blue",pch=3) +
geom_vline(xintercept=0,lty=2,col="grey") +
geom_hline(yintercept=1,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="gghybrid u (logit(centre)*v)", y = "gghybrid cline steepness v", tag = "B")


p33
###


#Same plot but plotting the bgc parameter estimates

###
p32 = ggplot()
p32 = p32 + geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95<=1],aes(x=-a_median,y=b_median),pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0],aes(x=-a_median,y=b_median),pch=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95>1],aes(x=-a_median,y=b_median),col="red",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue>=0.05&v_lower_95>1],aes(x=-a_median,y=b_median),col="blue",pch=1,cex=3) +
          geom_point(data=loci_parent[centre_pvalue<0.05&v_lower_95<=1],aes(x=-a_median,y=b_median),col="brown",pch=1,cex=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB>0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB>0],aes(x=-a_median,y=b_median),col="red",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB>0&a_0.950_CI_UB>0&b_0.950_CI_LB<=0|a_0.950_CI_LB<0&a_0.950_CI_UB<0&b_0.950_CI_LB<=0],aes(x=-a_median,y=b_median),col="brown",pch=3) +
          geom_point(data=loci_parent[a_0.950_CI_LB<0&a_0.950_CI_UB>0&b_0.950_CI_LB>0],aes(x=-a_median,y=b_median),col="blue",pch=3) +
geom_vline(xintercept=0,lty=2,col="grey") +
geom_hline(yintercept=0,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="bgc cline shift alpha", y = "bgc cline steepness beta", tag = "C")


p32
###
###


plot_grid(p3,p33,p32,nrow=1)













