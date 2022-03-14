################################################################################
#Statistical properties of ggcline tests########################################
################################################################################

################################################################################
#Richard Bailey 17 Feb 2022#####################################################
################################################################################
library(gghybrid)
library(ggplot2)
library(cowplot)
################################################################################


##################################################################################################################
##################################################################################################################
#2.3.3a High-resolution, evenly spaced hybrid indices. N loci all with centre=0.5 & v=1. Test false positive rate#
##################################################################################################################
##################################################################################################################

#1000 individuals with evenly distributed hybrid index values#

#1000 diploid loci fixed for alternate alleles in the parents#


#Hybrid index##############################################################################################

hybind=data.table(h_posterior_mode=seq(0,1,length=1000))#Assume this is a set of test individuals
hybind[,INDLABEL:=paste("A",seq(.N),sep="")]

hindex=list()
hindex$hi=hybind
hindex$hi[,Source:="TEST"]
hindex$test.subject="INDLABEL"

hindex

#Simulated locus data######################################################################################

simparam=data.table(v=rep(1,1000),centre=rep(0.5,1000))
simparam[,logv:=log(v)]

simparam[,index:=seq(1,.N)]
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]


#1 row per allele copy (diploid)#

indloc=expand.grid(hybind$INDLABEL,simparam$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(simparam,locus);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]

setkey(prepdata,locus)
prepdata=simparam[prepdata]

prepdata[,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,u:=qlogis(centre)*v]

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]

prepdata[is.na(Source_allele)]#No missing#

prepdata[,v:=NULL][,centre:=NULL][,logv:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#

prepdata[,Source:="TEST"]

####################
#Run genomic clines#
####################


###
gc_null=ggcline(
  data.prep.object=prepdata,    #Needs an entry#
  esth.object=hindex,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #All individuals are TEST#
  return.likmeans = TRUE,                 #For waic#
  read.data.precols="INDLABEL",          #Needs an entry#
  fix.subject.v = FALSE,
  fix.value.v,
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("L1","L1000"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-1, 2),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
###
gc_null
##########################################################################





###########################################
#Other models, for model comparison########
###########################################

locuslist=gc_null$gc$locus


###
gc_null_v_fixed=ggcline(
  data.prep.object=prepdata,    #Needs an entry#
  esth.object=hindex,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols="INDLABEL",          #Needs an entry#
  fix.subject.v = locuslist,
  fix.value.v=1,
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-1, 2),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)



###
gc_null_centre_fixed=ggcline(
  data.prep.object=prepdata,    #Needs an entry#
  esth.object=hindex,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols="INDLABEL",          #Needs an entry#
  fix.subject.v = FALSE,
  fix.value.v,
  fix.subject.centre = locuslist,
  fix.value.centre=0.5,
  plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-1, 2),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)


###
gc_null_both_fixed=ggcline(
  data.prep.object=prepdata,    #Needs an entry#
  esth.object=hindex,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols="INDLABEL",          #Needs an entry#
  fix.subject.v = locuslist,
  fix.value.v=1,
  fix.subject.centre = locuslist,
  fix.value.centre=0.5,
  plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-1, 2),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=2,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=0,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)

###################
###################
###################

##############
#AIC and waic#
##############

aic_gc_null=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_null,
 test.subject="locus")

aic_gc_null_centre_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_null_centre_fixed,
 test.subject="locus")

aic_gc_null_v_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_null_v_fixed,
 test.subject="locus")

aic_gc_null_both_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_null_both_fixed,
 test.subject="locus")



comp_gc_null_gc_null_centre_fixed=compare.models(
 ggcline.object1=gc_null,
 ggcline.object2=gc_null_centre_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_null_gc_null_v_fixed=compare.models(
 ggcline.object1=gc_null,
 ggcline.object2=gc_null_v_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_null_gc_null_both_fixed=compare.models(
 ggcline.object1=gc_null,
 ggcline.object2=gc_null_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_null_centre_fixed_gc_null_both_fixed=compare.models(
 ggcline.object1=gc_null_centre_fixed,
 ggcline.object2=gc_null_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_null_v_fixed_gc_null_both_fixed=compare.models(
 ggcline.object1=gc_null_v_fixed,
 ggcline.object2=gc_null_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)



#Combine the AIC results#
setnames(aic_gc_null,"AIC","AIC_gc_null")
setnames(aic_gc_null_centre_fixed,"AIC","AIC_gc_null_centre_fixed")
setnames(aic_gc_null_v_fixed,"AIC","AIC_gc_null_v_fixed")
setnames(aic_gc_null_both_fixed,"AIC","AIC_gc_null_both_fixed")
###
setkey(aic_gc_null,locus);setkey(aic_gc_null_centre_fixed,locus);setkey(aic_gc_null_v_fixed,locus);setkey(aic_gc_null_both_fixed,locus)
gc_nullwAIC=aic_gc_null[,.(locus,AIC_gc_null)][aic_gc_null_centre_fixed[,.(locus,AIC_gc_null_centre_fixed)]][aic_gc_null_v_fixed[,.(locus,AIC_gc_null_v_fixed)]][aic_gc_null_both_fixed[,.(locus,AIC_gc_null_both_fixed)]]

#Add p values from the full model#
setkey(gc_null$gc,locus)
gc_nullwAIC=gc_null$gc[,.(locus,v_pvalue,centre_pvalue)][gc_nullwAIC]

#Add the true v, logv, centre, logitcentre#ALL THE SAME#
#setkey(gcsim,locus);setkey(gcwAIC,locus)
#gcwAIC=gcsim[,.(locus,v,logv,centre,logitcentre)][gcwAIC]


#Now add the waic values#
setnames(comp_gc_null_gc_null_centre_fixed,"waic1_npar1_AICscale","waic_gc_null")
setnames(comp_gc_null_gc_null_centre_fixed,"waic2_npar2_AICscale","waic_gc_null_centre_fixed")
setnames(comp_gc_null_gc_null_v_fixed,"waic2_npar2_AICscale","waic_gc_null_v_fixed")
setnames(comp_gc_null_gc_null_both_fixed,"waic2_npar2_AICscale","waic_gc_null_both_fixed")

setkey(comp_gc_null_gc_null_centre_fixed,locus);setkey(comp_gc_null_gc_null_v_fixed,locus);setkey(comp_gc_null_gc_null_both_fixed,locus);setkey(gc_nullwAIC,locus)
gc_nullwAIC=gc_nullwAIC[comp_gc_null_gc_null_centre_fixed[,.(locus,waic_gc_null,waic_gc_null_centre_fixed)]][comp_gc_null_gc_null_v_fixed[,.(locus,waic_gc_null_v_fixed)]][comp_gc_null_gc_null_both_fixed[,.(locus,waic_gc_null_both_fixed)]]



#AIC comparisons#
gc_nullwAIC[AIC_gc_null_both_fixed<AIC_gc_null&AIC_gc_null_both_fixed<AIC_gc_null_centre_fixed&AIC_gc_null_both_fixed<AIC_gc_null_v_fixed]#676 both fixed is best#
gc_nullwAIC[AIC_gc_null_centre_fixed<AIC_gc_null&AIC_gc_null_centre_fixed<AIC_gc_null_both_fixed&AIC_gc_null_centre_fixed<AIC_gc_null_v_fixed]#160 centre fixed only is best#
gc_nullwAIC[AIC_gc_null_v_fixed<AIC_gc_null&AIC_gc_null_v_fixed<AIC_gc_null_both_fixed&AIC_gc_null_v_fixed<AIC_gc_null_centre_fixed]#142 v fixed only is best#
gc_nullwAIC[AIC_gc_null<AIC_gc_null_v_fixed&AIC_gc_null<AIC_gc_null_centre_fixed&AIC_gc_null<AIC_gc_null_both_fixed]#22 neither fixed is best#



#waic comparisons#
gc_nullwAIC[waic_gc_null_both_fixed<waic_gc_null&waic_gc_null_both_fixed<waic_gc_null_centre_fixed&waic_gc_null_both_fixed<waic_gc_null_v_fixed]#685 both fixed is best#
gc_nullwAIC[waic_gc_null_centre_fixed<waic_gc_null&waic_gc_null_centre_fixed<waic_gc_null_both_fixed&waic_gc_null_centre_fixed<waic_gc_null_v_fixed]#157 centre fixed only is best#
gc_nullwAIC[waic_gc_null_v_fixed<waic_gc_null&waic_gc_null_v_fixed<waic_gc_null_both_fixed&waic_gc_null_v_fixed<waic_gc_null_centre_fixed]#137 v fixed only is best#
gc_nullwAIC[waic_gc_null<waic_gc_null_v_fixed&waic_gc_null<waic_gc_null_centre_fixed&waic_gc_null<waic_gc_null_both_fixed]#21 neither fixed is best#

#############################################################################################################
#Plots#######################################################################################################
#############################################################################################################


###################################################
#Actual percentages favouring each of the 4 models#
###################################################

#####
#AIC#
#####

(676/1000)*100#67.6%

(160/1000)*100#16%

(142/1000)*100#14.2%

(22/1000)*100#2.2%



######
#waic#
######

(685/1000)*100#68.5%

(157/1000)*100#15.7%

(137/1000)*100#13.7%

(21/1000)*100#2.1%

###
mean(c(67.6,68.5))#68.05%
mean(c(16,14.2,15.7,13.7))#14.9%
mean(c(2.2,2.1))#2.15%


0.6805
0.298
0.0215


?rbinom

dbinom(c(0,1,2),size=2,prob=0.1)

pbinom(c(0.6805,0.298,0.0215),size=2,prob=0.1)

###
prs=seq(0,0.25,length=10001)

LL=matrix(0,length(prs),3)

for(i in 1:length(prs)){
LL[i,]=dbinom(c(0,1,2),size=2,prob=prs[i])
}
###

prLL=data.table(cbind(prs,LL))

prLL[,p0:=0.6805]
prLL[,p1:=0.298]
prLL[,p2:=0.0215]

prLL[,diff0:= (p0 - V2)^2]
prLL[,diff1:= (p1 - V3)^2]
prLL[,diff2:= (p2 - V4)^2]

prLL[,ssq:=diff0 + diff1 + diff2]


plot(ssq~prs,data=prLL)

prLL[ssq==min(ssq)]#So the overall probability of 1 success is 0.177#

##############################
#p values from the full model#
##############################

#
par(mfrow=c(1,2))
hist(gc_null$gc[,centre_pvalue],breaks=21)
abline(v=0.05,col="red",lty=2,lwd=2)
hist(gc_null$gc[,v_pvalue],breaks=21)
abline(v=0.05,col="red",lty=2,lwd=2)
#

gc_null$gc[centre_pvalue<=0.05]#46 loci
gc_null$gc[v_pvalue<=0.05]#63 loci


(46/1000)*100

(63/1000)*100


###################
###################
###################

#positive values mean the first model is favoured#

gc_nullwAIC[,diff_AIC_gc_null_gc_null_centre_fixed:= AIC_gc_null_centre_fixed - AIC_gc_null]
gc_nullwAIC[,diff_AIC_gc_null_gc_null_v_fixed:= AIC_gc_null_v_fixed - AIC_gc_null]
gc_nullwAIC[,diff_AIC_gc_null_gc_null_both_fixed:= AIC_gc_null_both_fixed - AIC_gc_null]

gc_nullwAIC[,diff_AIC_gc_null_v_fixed_gc_null_centre_fixed:= AIC_gc_null_centre_fixed - AIC_gc_null_v_fixed]
gc_nullwAIC[,diff_AIC_gc_null_centre_gc_null_v_fixed:= AIC_gc_null_v_fixed - AIC_gc_null_centre_fixed]

gc_nullwAIC[,diff_AIC_gc_null_v_fixed_gc_null_both_fixed:= AIC_gc_null_both_fixed - AIC_gc_null_v_fixed]
gc_nullwAIC[,diff_AIC_gc_null_centre_gc_null_both_fixed:= AIC_gc_null_both_fixed - AIC_gc_null_centre_fixed]

#

gc_nullwAIC[,diff_waic_gc_null_gc_null_centre_fixed:= waic_gc_null_centre_fixed - waic_gc_null]
gc_nullwAIC[,diff_waic_gc_null_gc_null_v_fixed:= waic_gc_null_v_fixed - waic_gc_null]
gc_nullwAIC[,diff_waic_gc_null_gc_null_both_fixed:= waic_gc_null_both_fixed - waic_gc_null]

gc_nullwAIC[,diff_waic_gc_null_v_fixed_gc_null_centre_fixed:= waic_gc_null_centre_fixed - waic_gc_null_v_fixed]
gc_nullwAIC[,diff_waic_gc_null_centre_gc_null_v_fixed:= waic_gc_null_v_fixed - waic_gc_null_centre_fixed]

gc_nullwAIC[,diff_waic_gc_null_v_fixed_gc_null_both_fixed:= waic_gc_null_both_fixed - waic_gc_null_v_fixed]
gc_nullwAIC[,diff_waic_gc_null_centre_gc_null_both_fixed:= waic_gc_null_both_fixed - waic_gc_null_centre_fixed]

#
gc_nullwAIC[diff_AIC_gc_null_gc_null_centre_fixed>2&diff_AIC_gc_null_gc_null_v_fixed>2&diff_AIC_gc_null_gc_null_both_fixed>2];#3 loci#
gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_v_fixed>2&diff_AIC_gc_null_centre_gc_null_both_fixed>2&diff_AIC_gc_null_gc_null_centre_fixed < -2];#3 loci#
gc_nullwAIC[diff_AIC_gc_null_v_fixed_gc_null_centre_fixed>2&diff_AIC_gc_null_v_fixed_gc_null_both_fixed>2&diff_AIC_gc_null_gc_null_v_fixed < -2];#2 loci#
#
gc_nullwAIC[diff_waic_gc_null_gc_null_centre_fixed>2&diff_waic_gc_null_gc_null_v_fixed>2&diff_waic_gc_null_gc_null_both_fixed>2];#2 loci#
gc_nullwAIC[diff_waic_gc_null_centre_gc_null_v_fixed>2&diff_waic_gc_null_centre_gc_null_both_fixed>2&diff_waic_gc_null_gc_null_centre_fixed < -2];#3 loci#
gc_nullwAIC[diff_waic_gc_null_v_fixed_gc_null_centre_fixed>2&diff_waic_gc_null_v_fixed_gc_null_both_fixed>2&diff_waic_gc_null_gc_null_v_fixed < -2];#5 loci#


###

breaks=seq(-5,16,length=21)

#
par(mfrow=c(2,2))
#
plot(gc_nullwAIC[,diff_AIC_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[,diff_AIC_gc_null_v_fixed_gc_null_both_fixed],type="n",xlim=c(-5,16),ylim=c(-5,16))
points(gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed<0&diff_AIC_gc_null_v_fixed_gc_null_both_fixed<0,diff_AIC_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed<0&diff_AIC_gc_null_v_fixed_gc_null_both_fixed<0,diff_AIC_gc_null_v_fixed_gc_null_both_fixed],col="green",cex=0.5,pch=16)
points(gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=0|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=0,diff_AIC_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=0|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=0,diff_AIC_gc_null_v_fixed_gc_null_both_fixed],col="lightblue",cex=0.5,pch=16)
points(gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=2|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=2,diff_AIC_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=2|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=2,diff_AIC_gc_null_v_fixed_gc_null_both_fixed],col="blue",cex=0.5,pch=16)
points(gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=4|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=4,diff_AIC_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed>=4|diff_AIC_gc_null_v_fixed_gc_null_both_fixed>=4,diff_AIC_gc_null_v_fixed_gc_null_both_fixed],col="red",cex=0.5,pch=16)
#
plot(gc_nullwAIC[,diff_waic_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[,diff_waic_gc_null_v_fixed_gc_null_both_fixed],type="n",xlim=c(-5,16),ylim=c(-5,16))
points(gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed<0&diff_waic_gc_null_v_fixed_gc_null_both_fixed<0,diff_waic_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed<0&diff_waic_gc_null_v_fixed_gc_null_both_fixed<0,diff_waic_gc_null_v_fixed_gc_null_both_fixed],col="green",cex=0.5,pch=16)
points(gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=0|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=0,diff_waic_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=0|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=0,diff_waic_gc_null_v_fixed_gc_null_both_fixed],col="lightblue",cex=0.5,pch=16)
points(gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=2|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=2,diff_waic_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=2|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=2,diff_waic_gc_null_v_fixed_gc_null_both_fixed],col="blue",cex=0.5,pch=16)
points(gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=4|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=4,diff_waic_gc_null_centre_gc_null_both_fixed]~gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed>=4|diff_waic_gc_null_v_fixed_gc_null_both_fixed>=4,diff_waic_gc_null_v_fixed_gc_null_both_fixed],col="red",cex=0.5,pch=16)
#
hist(gc_nullwAIC[,diff_AIC_gc_null_gc_null_both_fixed],breaks=breaks,xlim=c(-5,16),ylim=c(0,400))
abline(v=c(0,2,4),col=c("lightblue","blue","red"))
#
hist(gc_nullwAIC[,diff_waic_gc_null_gc_null_both_fixed],breaks=breaks,xlim=c(-5,16),ylim=c(0,400))
abline(v=c(0,2,4),col=c("lightblue","blue","red"))
#


gc_nullwAIC[order(diff_AIC_gc_null_gc_null_both_fixed)]
gc_nullwAIC[order(diff_waic_gc_null_gc_null_both_fixed)]




###############################
###############################
#Try something else############
###############################
###############################

#These are the values presented in the paper**********

gc_nullwAIC[diff_AIC_gc_null_v_fixed_gc_null_both_fixed >= 2,.N]#43 = 4.3%
gc_nullwAIC[diff_AIC_gc_null_centre_gc_null_both_fixed >= 2,.N]#51 = 5.1%

gc_nullwAIC[diff_waic_gc_null_v_fixed_gc_null_both_fixed >= 2,.N]#43 = 4.3%
gc_nullwAIC[diff_waic_gc_null_centre_gc_null_both_fixed >= 2,.N]#51 = 5.1%




gc_nullwAIC[diff_AIC_gc_null_v_fixed_gc_null_both_fixed > 0|diff_AIC_gc_null_centre_gc_null_both_fixed > 0]#324

gc_nullwAIC[diff_AIC_gc_null_v_fixed_gc_null_both_fixed > 0&diff_AIC_gc_null_centre_gc_null_both_fixed > 0]#23


(43/1000)*100


