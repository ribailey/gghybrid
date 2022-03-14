###########################################################
#2.3.3c Grid of centre & v values, testing all 4 models####
###########################################################

###########################################################
#Richard Bailey 20 Feb 2022################################
###########################################################
library(gghybrid)
library(ggplot2)
library(cowplot)
###########################################################

#centre (0.000001 to 0.999999) and ln(v) (-4 to 4) 

#Evenly distributed centres and widths#

vs=seq(-4,4,length=49)
expvs=exp(vs)

centres=seq(0.000001,0.999999,length=49)

centres==0.5

simparam=expand.grid(expvs,centres)

head(simparam)

simparam=data.table(simparam);setnames(simparam,c("Var1","Var2"),c("v","centre"))

simparam[centre>0.49&centre<0.51]
simparam[v==1]

########################################
#Plot the true log(v) and centre values#
########################################

simparam[,logv:=log(v)]

base <- ggplot() + xlim(0, 1) + theme_bw()
base + geom_point(aes(x=centre,y=logv),data=simparam,cex=0.1)




##########################################
#Regularly spaced hybrid indices##########
##########################################

hybind=data.table(h_posterior_mode=seq(0,1,length=1000))#Assume this is a set of test individuals

hybind[,INDLABEL:=paste("A",seq(.N),sep="")]


#########################################
#What about parental allele frequencies?#
#########################################

################
#All loci fixed#
################

simparam[,S0_af1.2:=0][,S1_af1.2:=1]


#Add locus names#

simparam[,index:=seq(1,.N)]
simparam[,locus:=paste("L",index,sep="")]




############################################################
############################################################
#Now simulate the individual genotypes, assuming no linkage#
############################################################
############################################################


#######################################
#Now I need to create a data.prep file#
#######################################


simparam[,index:=NULL]#Only used to create the locus name#

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(simparam,locus);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]

setkey(prepdata,locus)
prepdata=simparam[prepdata]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,u:=qlogis(centre)*v]

prepdata[,Source_allele:=rbinom(.N,1,prob=S0_af1.2 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1_af1.2 - S0_af1.2))]

prepdata[is.na(Source_allele)]#4 failures: all log(v)=4, hindex 0 or 1, centre one or other extreme. Remove these#

prepdata=prepdata[is.na(Source_allele)==FALSE]

setnames(prepdata,c("S0_af1.2","S1_af1.2"),c("S0.prop_1","S1.prop_1"))

prepdata[,v:=NULL][,centre:=NULL][,logv:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#

prepdata[,Source:="TEST"]

#################################

hindex=list()
hindex$hi=hybind
hindex$hi[,Source:="TEST"]
hindex$test.subject="INDLABEL"


simparam[v<5&v>2&centre>0.3&centre<0.7]

#####################################################################
#Run genomic clines full model on all 2401 loci and 1000 individuals#
#####################################################################

gc=ggcline(
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
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("L1255","L965"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-2, 6),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)), 
  prior.logv=c(0,sqrt(10)), 
  nitt=5000,                              #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                            #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)

gc
##########################################################################


#
par(mfrow=c(1,2))
#
plot(gc$gc[,-log10(v_pvalue + 1e-200)]~gc$gc[,mean_log_v],pch=16,col="blue",cex=0.1)
abline(h=-log10(0.05),v=0,col="red",lty=2)
#
plot(gc$gc[,-log10(centre_pvalue + 1e-200)]~gc$gc[,invlogit_mean_logit_centre],pch=16,col="blue",cex=0.1)
abline(h=-log10(0.05),v=0,col="red",lty=2)
#



gc$gc[v_pvalue==0]#None
gc$gc[centre_pvalue==0]#655


setkey(simparam,locus);setkey(gc$gc,locus)

gcsim=simparam[gc$gc]

###

#fdrlevel=0.05

#gcsim[,width_lower:=plogis(mean_logit_centre - (1/mean_log_v))]
#gcsim[,width_upper:=plogis(mean_logit_centre + (1/mean_log_v))]
#gcsim[,v_fdrSig:=qvalue(p = gcsim[,v_pvalue],fdr.level=fdrlevel, pfdr = TRUE)$significant]            #REMEMBER SOME WILL BE NA#
#Error in smooth.spline(lambda, pi0, df = smooth.df) : 
#  missing or infinite values in inputs are not allowed
#gcsim[,centre_fdrSig:=qvalue(p = gcsim[,centre_pvalue],fdr.level=fdrlevel, pfdr = TRUE)$significant]  #REMEMBER SOME WILL BE NA#




######################################################################################
######################################################################################
######################################################################################

#############################################################
#Make plots for parameter accuracy and parameter uncertainty#
#############################################################

gcsim[,logitcentre:=qlogis(centre)]

#gcsim[,v_range:=v_upper_95 - v_lower_95]
#gcsim[,centre_range:=centre_upper_95 - centre_lower_95]

#######################
#Parameter uncertainty#
#######################

#I DON'T THINK THIS IS USEFUL**********************

###
par(mfrow=c(1,2))
###LATENT PARAMETER TRUE VALUES AND POSTERIOR VARIANCES###
plot(var_logit_centre~logitcentre,data=gcsim,type="n")
points(var_logit_centre~logitcentre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(var_logit_centre~logitcentre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(var_logit_centre~logitcentre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
plot(var_log_v~logv,data=gcsim,type="n")
points(var_log_v~logv,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(var_log_v~logv,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(var_log_v~logv,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")
###DATA SCALE PARAMETER TRUE VALUES AND 95% POSTERIOR CI RANGE###
plot(centre_range~centre,data=gcsim,type="n")
points(centre_range~centre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(centre_range~centre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(centre_range~centre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
plot(v_range~v,data=gcsim,type="n")
points(v_range~v,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(v_range~v,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(v_range~v,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")
###


####################
#Parameter accuracy#
####################

#Include both (1) the error in the best posterior estimate and (2) the quantile of the true value wrt to the posterior distribution.

#Also include plots of actual values against estimated values.

gcsim[,logitcentre_squarederror:=(logitcentre - mean_logit_centre)^2]
gcsim[,logv_squarederror:=(logv - mean_log_v)^2]

###
par(mfrow=c(2,2))
###SQUARED ERROR FOR LATENT PARAMETER ESTIMATES###
#plot(logitcentre_squarederror~logitcentre,data=gcsim,type="n")
points(logitcentre_squarederror~logitcentre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(logitcentre_squarederror~logitcentre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(logitcentre_squarederror~logitcentre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
#plot(logv_squarederror~logv,data=gcsim,type="n")
points(logv_squarederror~logv,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(logv_squarederror~logv,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(logv_squarederror~logv,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")
###TRUE AGAINST ESTIMATED LATENT PARAMETER VALUES##
plot(mean_logit_centre~logitcentre,data=gcsim,type="n")
abline(a=0,b=1,col="grey",lty=2,lwd=2)
points(mean_logit_centre~logitcentre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(mean_logit_centre~logitcentre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(mean_logit_centre~logitcentre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
plot(mean_log_v~logv,data=gcsim,type="n")
abline(a=0,b=1,col="grey",lty=2,lwd=2)
points(mean_log_v~logv,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(mean_log_v~logv,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(mean_log_v~logv,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")
###TRUE AGAINST ESTIMATED DATA SCALE PARAMETER VALUES##
plot(invlogit_mean_logit_centre~centre,data=gcsim,type="n")
abline(a=0,b=1,col="grey",lty=2,lwd=2)
points(invlogit_mean_logit_centre~centre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(invlogit_mean_logit_centre~centre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(invlogit_mean_logit_centre~centre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
plot(exp_mean_log_v~v,data=gcsim,type="n")
abline(a=0,b=1,col="grey",lty=2,lwd=2)
points(exp_mean_log_v~v,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(exp_mean_log_v~v,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(exp_mean_log_v~v,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")


###################
#Plot of quantiles#
###################

gcsim[mean_logit_centre<=logitcentre,
    centre_true_quant:=(pnorm(logitcentre,mean=mean_logit_centre,sd=sqrt(var_logit_centre),lower.tail=FALSE))*2][mean_logit_centre>logitcentre,
    centre_true_quant:=(pnorm(logitcentre,mean=mean_logit_centre,sd=sqrt(var_logit_centre),lower.tail=TRUE))*2];

gcsim[mean_log_v<=logv,
    v_true_quant:=(pnorm(logv,mean=mean_log_v,sd=sqrt(var_log_v),lower.tail=FALSE))*2][mean_log_v>logv,
    v_true_quant:=(pnorm(logv,mean=mean_log_v,sd=sqrt(var_log_v),lower.tail=TRUE))*2];




###
par(mfrow=c(2,2))
###QUANTILE OF TRUE VALUE WRT POSTERIOR DISTRIBUTION###
plot(-log10(centre_true_quant)~logitcentre,data=gcsim,type="n")
abline(h=-log10(0.05),col="grey",lty=2,lwd=2)
points(-log10(centre_true_quant)~logitcentre,data=gcsim[logv >= -2&logv <= 2],pch=16,cex=0.5)
points(-log10(centre_true_quant)~logitcentre,data=gcsim[logv < -2],pch=16,cex=0.5,col="blue")
points(-log10(centre_true_quant)~logitcentre,data=gcsim[logv > 2],pch=16,cex=0.5,col="red")
###
plot(-log10(v_true_quant)~logv,data=gcsim,type="n")
abline(h=-log10(0.05),col="grey",lty=2,lwd=2)
points(-log10(v_true_quant)~logv,data=gcsim[logitcentre>= -5&logitcentre<= 5],pch=16,cex=0.5)
points(-log10(v_true_quant)~logv,data=gcsim[logitcentre< -5],pch=16,cex=0.5,col="orange")
points(-log10(v_true_quant)~logv,data=gcsim[logitcentre> 5],pch=16,cex=0.5,col="magenta")
###
hist(gcsim[,centre_true_quant],breaks=21)
abline(v=0.05,col="red",lty=2,lwd=2)
###
hist(gcsim[,v_true_quant],breaks=21)
abline(v=0.05,col="red",lty=2,lwd=2)
###






######################################
######################################
#ggplot for manuscript################
######################################
######################################

#Include the percentage in the bottom quantile (see the two histograms above) in the results text.

###
p1=ggplot()
p1=p1 + geom_abline(intercept=0,slope=1,lty=2,col="grey",lwd=1) +
 geom_point(data=gcsim[logv >= -2&logv <= 2],aes(y=invlogit_mean_logit_centre,x=centre),col="black",cex=0.9) +
 geom_point(data=gcsim[logv < -2],aes(y=invlogit_mean_logit_centre,x=centre),col="magenta",cex=0.9, alpha = 0.5) +
 geom_point(data=gcsim[logv > 2],aes(y=invlogit_mean_logit_centre,x=centre),col="red",cex=0.9, alpha = 0.5) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 #labs(x="true centre", y = "best posterior centre", tag = "A")
 labs(x=NULL,y = "best posterior centre", tag = "A")

#p1
###
p2=ggplot()
p2=p2 + geom_abline(intercept=0,slope=1,lty=2,col="grey",lwd=1) +
 geom_point(data=gcsim[logitcentre>= -5&logitcentre<= 5],aes(y=exp_mean_log_v,x=v),col="black",cex=0.9) +
 geom_point(data=gcsim[logitcentre < -5],aes(y=exp_mean_log_v,x=v),col="purple",cex=0.9) +
 geom_point(data=gcsim[logitcentre > 5],aes(y=exp_mean_log_v,x=v),col="orange",cex=0.9) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x=NULL, y = "best posterior v", tag = "B")

#p2
###
p3=ggplot()
p3=p3 + geom_hline(yintercept= -log10(0.05),lty=2,col="grey",lwd=1) +
 geom_point(data=gcsim[logv >= -2&logv <= 2],aes(y=-log10(centre_true_quant),x=centre),col="black",cex=0.9) +
 geom_point(data=gcsim[logv < -2],aes(y=-log10(centre_true_quant),x=centre),col="magenta",cex=0.9, alpha = 0.5) +
 geom_point(data=gcsim[logv > 2],aes(y=-log10(centre_true_quant),x=centre),col="red",cex=0.9, alpha = 0.5) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true centre", y = "-log10(true centre quantile)", tag = "C")

#p3
###
p4=ggplot()
p4=p4 + geom_hline(yintercept= -log10(0.05),lty=2,col="grey",lwd=1) +
 geom_point(data=gcsim[logitcentre>= -5&logitcentre<= 5],aes(y=-log10(v_true_quant),x=v),col="black",cex=0.9) +
 geom_point(data=gcsim[logitcentre < -5],aes(y=-log10(v_true_quant),x=v),col="purple",cex=0.9) +
 geom_point(data=gcsim[logitcentre > 5],aes(y=-log10(v_true_quant),x=v),col="orange",cex=0.9) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true v", y = "-log10(true v quantile)", tag = "D")

#p4
###

plot_grid(p1,p2,p3,p4,nrow=2)






####################################################################################################
#Now run again with each parameter fixed to null value, and both parameters fixed, for waic and aic#
####################################################################################################

locuslist=gcsim$locus

gc_centre_fixed=ggcline(
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
  fix.value.centre = 0.5,
  #plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.ylim = c(-2, 3),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
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


gc_v_fixed=ggcline(
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
  fix.value.v = 1,
  fix.subject.centre = FALSE,
  fix.value.centre,
  #plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.ylim = c(-2, 3),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
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



gc_both_fixed=ggcline(
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
  fix.value.v = 1,
  fix.subject.centre = locuslist,
  fix.value.centre = 0.5,
  #plot.test.subject = c("L279","L529"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.ylim = c(-2, 3),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
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


##############################################
##############################################
##############################################

aic_gc=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc,
 test.subject="locus")

aic_gc_centre_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_centre_fixed,
 test.subject="locus")

aic_gc_v_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_v_fixed,
 test.subject="locus")

aic_gc_both_fixed=calc_AIC(
 data.prep.object=prepdata,
 esth.object=hindex,
 ggcline.object=gc_both_fixed,
 test.subject="locus")



comp_gc_gc_centre_fixed=compare.models(
 ggcline.object1=gc,
 ggcline.object2=gc_centre_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_gc_v_fixed=compare.models(
 ggcline.object1=gc,
 ggcline.object2=gc_v_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_gc_both_fixed=compare.models(
 ggcline.object1=gc,
 ggcline.object2=gc_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_centre_fixed_gc_both_fixed=compare.models(
 ggcline.object1=gc_centre_fixed,
 ggcline.object2=gc_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_gc_v_fixed_gc_both_fixed=compare.models(
 ggcline.object1=gc_v_fixed,
 ggcline.object2=gc_both_fixed,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)



#Combine the AIC results#
setnames(aic_gc,"AIC","AIC_gc")
setnames(aic_gc_centre_fixed,"AIC","AIC_gc_centre_fixed")
setnames(aic_gc_v_fixed,"AIC","AIC_gc_v_fixed")
setnames(aic_gc_both_fixed,"AIC","AIC_gc_both_fixed")

setkey(aic_gc,locus);setkey(aic_gc_centre_fixed,locus);setkey(aic_gc_v_fixed,locus);setkey(aic_gc_both_fixed,locus)
gcwAIC=aic_gc[,.(locus,AIC_gc)][aic_gc_centre_fixed[,.(locus,AIC_gc_centre_fixed)]][aic_gc_v_fixed[,.(locus,AIC_gc_v_fixed)]][aic_gc_both_fixed[,.(locus,AIC_gc_both_fixed)]]

#Add p values from the full model#
setkey(gc$gc,locus)
gcwAIC=gc$gc[,.(locus,v_pvalue,centre_pvalue)][gcwAIC]

#Add the true v, logv, centre, logitcentre#
setkey(gcsim,locus);setkey(gcwAIC,locus)
gcwAIC=gcsim[,.(locus,v,logv,centre,logitcentre)][gcwAIC]


#Now add the waic values#
setnames(comp_gc_gc_centre_fixed,"waic1_npar1_AICscale","waic_gc")
setnames(comp_gc_gc_centre_fixed,"waic2_npar2_AICscale","waic_gc_centre_fixed")
setnames(comp_gc_gc_v_fixed,"waic2_npar2_AICscale","waic_gc_v_fixed")
setnames(comp_gc_gc_both_fixed,"waic2_npar2_AICscale","waic_gc_both_fixed")

setkey(comp_gc_gc_centre_fixed,locus);setkey(comp_gc_gc_v_fixed,locus);setkey(comp_gc_gc_both_fixed,locus);setkey(gcwAIC,locus)
gcwAIC=gcwAIC[comp_gc_gc_centre_fixed[,.(locus,waic_gc,waic_gc_centre_fixed)]][comp_gc_gc_v_fixed[,.(locus,waic_gc_v_fixed)]][comp_gc_gc_both_fixed[,.(locus,waic_gc_both_fixed)]]

################################################
################################################
################################################

#AIC comparisons#
gcwAIC[AIC_gc_both_fixed<AIC_gc&AIC_gc_both_fixed<AIC_gc_centre_fixed&AIC_gc_both_fixed<AIC_gc_v_fixed]#none with null as best model#
gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed]#407 with centre fixed as best model#
gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed]#101 with v fixed as best model#
gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed]#1893 with full model as best model#


#waic comparisons#
gcwAIC[waic_gc_both_fixed<waic_gc&waic_gc_both_fixed<waic_gc_centre_fixed&waic_gc_both_fixed<waic_gc_v_fixed]#none with null as best model#
gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed]#164 with centre fixed as best model#
gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed]#100 with v fixed as best model#
gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed]#2137 with full model as best model#

#############################################################################################################
#Plots#######################################################################################################
#############################################################################################################

#For information criteria, these just show which is the optimal model, even if only marginally#

#
par(mfrow=c(1,3))
#p values#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")
points(gcwAIC[v_pvalue>0.05&centre_pvalue>0.05,logv]~gcwAIC[v_pvalue>0.05&centre_pvalue>0.05,centre],col="grey",pch=16)
points(gcwAIC[v_pvalue<=0.05&centre_pvalue>0.05,logv]~gcwAIC[v_pvalue<=0.05&centre_pvalue>0.05,centre],col="blue",pch=16)
points(gcwAIC[v_pvalue>0.05&centre_pvalue<=0.05,logv]~gcwAIC[v_pvalue>0.05&centre_pvalue<=0.05,centre],col="green",pch=16)
points(gcwAIC[v_pvalue<=0.05&centre_pvalue<=0.05,logv]~gcwAIC[v_pvalue<=0.05&centre_pvalue<=0.05,centre],col="red",pch=16)

#AIC#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")
points(gcwAIC[AIC_gc_both_fixed<AIC_gc&AIC_gc_both_fixed<AIC_gc_centre_fixed&AIC_gc_both_fixed<AIC_gc_v_fixed,logv]~gcwAIC[AIC_gc_both_fixed<AIC_gc&AIC_gc_both_fixed<AIC_gc_centre_fixed&AIC_gc_both_fixed<AIC_gc_v_fixed,centre],col="grey",pch=16)
points(gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed,logv]~gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed,centre],col="green",pch=16)
points(gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed,logv]~gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed,centre],col="blue",pch=16)
points(gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed,logv]~gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed,centre],col="red",pch=16)

#waic#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")
points(gcwAIC[waic_gc_both_fixed<waic_gc&waic_gc_both_fixed<waic_gc_centre_fixed&waic_gc_both_fixed<waic_gc_v_fixed,logv]~gcwAIC[waic_gc_both_fixed<waic_gc&waic_gc_both_fixed<waic_gc_centre_fixed&waic_gc_both_fixed<waic_gc_v_fixed,centre],col="grey",pch=16)
points(gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed,logv]~gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed,centre],col="green",pch=16)
points(gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed,logv]~gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed,centre],col="blue",pch=16)
points(gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed,logv]~gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed,centre],col="red",pch=16)




###################
###################
###################

#positive values mean the first model is favoured#

gcwAIC[,diff_AIC_gc_gc_centre_fixed:= AIC_gc_centre_fixed - AIC_gc]
gcwAIC[,diff_AIC_gc_gc_v_fixed:= AIC_gc_v_fixed - AIC_gc]
gcwAIC[,diff_AIC_gc_gc_both_fixed:= AIC_gc_both_fixed - AIC_gc]

gcwAIC[,diff_AIC_gc_v_fixed_gc_centre_fixed:= AIC_gc_centre_fixed - AIC_gc_v_fixed]
gcwAIC[,diff_AIC_gc_centre_gc_v_fixed:= AIC_gc_v_fixed - AIC_gc_centre_fixed]

gcwAIC[,diff_AIC_gc_v_fixed_gc_both_fixed:= AIC_gc_both_fixed - AIC_gc_v_fixed]
gcwAIC[,diff_AIC_gc_centre_gc_both_fixed:= AIC_gc_both_fixed - AIC_gc_centre_fixed]


setkey(gcwAIC,diff_AIC_gc_gc_centre_fixed)

gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2]


plot(gcwAIC[,logv]~gcwAIC[,centre],col="grey",pch=16)
points(gcwAIC[diff_AIC_gc_gc_centre_fixed>2&diff_AIC_gc_gc_v_fixed>2&diff_AIC_gc_gc_both_fixed>2,logv]~gcwAIC[diff_AIC_gc_gc_centre_fixed>2&diff_AIC_gc_gc_v_fixed>2&diff_AIC_gc_gc_both_fixed>2,centre],col="red",pch=16)
points(gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2,logv]~gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2,centre],col="blue",pch=16)
points(gcwAIC[diff_AIC_gc_v_fixed_gc_centre_fixed>2&diff_AIC_gc_v_fixed_gc_both_fixed>2&diff_AIC_gc_gc_v_fixed < -2,logv]~gcwAIC[diff_AIC_gc_v_fixed_gc_centre_fixed>2&diff_AIC_gc_v_fixed_gc_both_fixed>2&diff_AIC_gc_gc_v_fixed < -2,centre],col="green",pch=16)


#In this final plot, grey means there's no clear winner in terms of estimating v only, centre only, or both#


gcwAIC[,diff_waic_gc_gc_centre_fixed:= waic_gc_centre_fixed - waic_gc]
gcwAIC[,diff_waic_gc_gc_v_fixed:= waic_gc_v_fixed - waic_gc]
gcwAIC[,diff_waic_gc_gc_both_fixed:= waic_gc_both_fixed - waic_gc]

gcwAIC[,diff_waic_gc_v_fixed_gc_centre_fixed:= waic_gc_centre_fixed - waic_gc_v_fixed]
gcwAIC[,diff_waic_gc_centre_gc_v_fixed:= waic_gc_v_fixed - waic_gc_centre_fixed]

gcwAIC[,diff_waic_gc_v_fixed_gc_both_fixed:= waic_gc_both_fixed - waic_gc_v_fixed]
gcwAIC[,diff_waic_gc_centre_gc_both_fixed:= waic_gc_both_fixed - waic_gc_centre_fixed]

################################
#Putting all the plots together#
################################

par(mfrow=c(2,3))

#AIC#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")
abline(v=0.5,h=0,lty=2)
points(gcwAIC[AIC_gc_both_fixed<AIC_gc&AIC_gc_both_fixed<AIC_gc_centre_fixed&AIC_gc_both_fixed<AIC_gc_v_fixed,logv]~gcwAIC[AIC_gc_both_fixed<AIC_gc&AIC_gc_both_fixed<AIC_gc_centre_fixed&AIC_gc_both_fixed<AIC_gc_v_fixed,centre],col="grey",pch=16)
points(gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed,logv]~gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed,centre],col="green",pch=16)
points(gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed,logv]~gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed,centre],col="blue",pch=16)
points(gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed,logv]~gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed,centre],col="red",pch=16)

#waic#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")
abline(v=0.5,h=0,lty=2)
points(gcwAIC[waic_gc_both_fixed<waic_gc&waic_gc_both_fixed<waic_gc_centre_fixed&waic_gc_both_fixed<waic_gc_v_fixed,logv]~gcwAIC[waic_gc_both_fixed<waic_gc&waic_gc_both_fixed<waic_gc_centre_fixed&waic_gc_both_fixed<waic_gc_v_fixed,centre],col="grey",pch=16)
points(gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed,logv]~gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed,centre],col="green",pch=16)
points(gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed,logv]~gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed,centre],col="blue",pch=16)
points(gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed,logv]~gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed,centre],col="red",pch=16)


#
#p values#
plot(gcwAIC[,logv]~gcwAIC[,centre],type="n")

points(gcwAIC[v_pvalue>0.01&centre_pvalue>0.01,logv]~gcwAIC[v_pvalue>0.01&centre_pvalue>0.01,centre],col="grey",pch=16)
points(gcwAIC[v_pvalue<=0.01&centre_pvalue>0.01,logv]~gcwAIC[v_pvalue<=0.01&centre_pvalue>0.01,centre],col="blue",pch=16)
points(gcwAIC[v_pvalue>0.01&centre_pvalue<=0.01,logv]~gcwAIC[v_pvalue>0.01&centre_pvalue<=0.01,centre],col="green",pch=16)
points(gcwAIC[v_pvalue<=0.01&centre_pvalue<=0.01,logv]~gcwAIC[v_pvalue<=0.01&centre_pvalue<=0.01,centre],col="red",pch=16)
abline(v=0.5,h=0,lty=2)

#
plot(gcwAIC[,logv]~gcwAIC[,centre],col="grey",pch=16)
abline(v=0.5,h=0,lty=2)
points(gcwAIC[diff_AIC_gc_gc_centre_fixed>2&diff_AIC_gc_gc_v_fixed>2&diff_AIC_gc_gc_both_fixed>2,logv]~gcwAIC[diff_AIC_gc_gc_centre_fixed>2&diff_AIC_gc_gc_v_fixed>2&diff_AIC_gc_gc_both_fixed>2,centre],col="red",pch=16)
points(gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2,logv]~gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2,centre],col="blue",pch=16)
points(gcwAIC[diff_AIC_gc_v_fixed_gc_centre_fixed>2&diff_AIC_gc_v_fixed_gc_both_fixed>2&diff_AIC_gc_gc_v_fixed < -2,logv]~gcwAIC[diff_AIC_gc_v_fixed_gc_centre_fixed>2&diff_AIC_gc_v_fixed_gc_both_fixed>2&diff_AIC_gc_gc_v_fixed < -2,centre],col="green",pch=16)
#
plot(gcwAIC[,logv]~gcwAIC[,centre],col="grey",pch=16)
abline(v=0.5,h=0,lty=2)
points(gcwAIC[diff_waic_gc_gc_centre_fixed>2&diff_waic_gc_gc_v_fixed>2&diff_waic_gc_gc_both_fixed>2,logv]~gcwAIC[diff_waic_gc_gc_centre_fixed>2&diff_waic_gc_gc_v_fixed>2&diff_waic_gc_gc_both_fixed>2,centre],col="red",pch=16)
points(gcwAIC[diff_waic_gc_centre_gc_v_fixed>2&diff_waic_gc_centre_gc_both_fixed>2&diff_waic_gc_gc_centre_fixed < -2,logv]~gcwAIC[diff_waic_gc_centre_gc_v_fixed>2&diff_waic_gc_centre_gc_both_fixed>2&diff_waic_gc_gc_centre_fixed < -2,centre],col="blue",pch=16)
points(gcwAIC[diff_waic_gc_v_fixed_gc_centre_fixed>2&diff_waic_gc_v_fixed_gc_both_fixed>2&diff_waic_gc_gc_v_fixed < -2,logv]~gcwAIC[diff_waic_gc_v_fixed_gc_centre_fixed>2&diff_waic_gc_v_fixed_gc_both_fixed>2&diff_waic_gc_gc_v_fixed < -2,centre],col="green",pch=16)




#################################
#################################
#Try it in ggplot################
#################################
#################################

palette.colors(palette = "Okabe-Ito")





###AIC###
p1=ggplot()
p1=p1 + geom_point(data=gcwAIC[AIC_gc_v_fixed<AIC_gc&AIC_gc_v_fixed<AIC_gc_both_fixed&AIC_gc_v_fixed<AIC_gc_centre_fixed],aes(x=centre,y=logv),col="#009E73",cex=1,pch=15) + #bluish green#
 geom_point(data=gcwAIC[AIC_gc_centre_fixed<AIC_gc&AIC_gc_centre_fixed<AIC_gc_both_fixed&AIC_gc_centre_fixed<AIC_gc_v_fixed],aes(x=centre,y=logv),col="#0072B2",cex=1,pch=15) + #blue#
 geom_point(data=gcwAIC[AIC_gc<AIC_gc_v_fixed&AIC_gc<AIC_gc_centre_fixed&AIC_gc<AIC_gc_both_fixed],aes(x=centre,y=logv),col="#CC79A7",cex=1,pch=15) + #reddish purple#
 geom_hline(yintercept=0,lty=2,col="black") +
 geom_vline(xintercept=0.5,lty=2,col="black") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="", y = "true ln(v)", tag = "A")

#p1
###waic###
p2=ggplot()
p2=p2 + geom_point(data=gcwAIC[waic_gc_v_fixed<waic_gc&waic_gc_v_fixed<waic_gc_both_fixed&waic_gc_v_fixed<waic_gc_centre_fixed],aes(x=centre,y=logv),col="#009E73",cex=1,pch=15) +
 geom_point(data=gcwAIC[waic_gc_centre_fixed<waic_gc&waic_gc_centre_fixed<waic_gc_both_fixed&waic_gc_centre_fixed<waic_gc_v_fixed],aes(x=centre,y=logv),col="#0072B2",cex=1,pch=15) +
 geom_point(data=gcwAIC[waic_gc<waic_gc_v_fixed&waic_gc<waic_gc_centre_fixed&waic_gc<waic_gc_both_fixed],aes(x=centre,y=logv),col="#CC79A7",cex=1,pch=15) +
 geom_hline(yintercept=0,lty=2,col="black") +
 geom_vline(xintercept=0.5,lty=2,col="black") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="", y = "true ln(v)", tag = "C")

#p2
###


###AIC>2###
p11=ggplot()
p11=p11 + geom_point(data=gcwAIC,aes(x=centre,y=logv),col="#999999",cex=1,pch=15) + #grey#
 geom_point(data=gcwAIC[diff_AIC_gc_v_fixed_gc_centre_fixed>2&diff_AIC_gc_v_fixed_gc_both_fixed>2&diff_AIC_gc_gc_v_fixed < -2],aes(x=centre,y=logv),col="#009E73",cex=1,pch=15) + #bluish green#
 geom_point(data=gcwAIC[diff_AIC_gc_centre_gc_v_fixed>2&diff_AIC_gc_centre_gc_both_fixed>2&diff_AIC_gc_gc_centre_fixed < -2],aes(x=centre,y=logv),col="#0072B2",cex=1,pch=15) + #blue#
 geom_point(data=gcwAIC[diff_AIC_gc_gc_centre_fixed>2&diff_AIC_gc_gc_v_fixed>2&diff_AIC_gc_gc_both_fixed>2],aes(x=centre,y=logv),col="#CC79A7",cex=1,pch=15) + #reddish purple#
 geom_hline(yintercept=0,lty=2,col="black") +
 geom_vline(xintercept=0.5,lty=2,col="black") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="", y = "", tag = "B")

#p11
###waic>2###
p22=ggplot()
p22=p22 + geom_point(data=gcwAIC,aes(x=centre,y=logv),col="#999999",cex=1,pch=15) + 
 geom_point(data=gcwAIC[diff_waic_gc_v_fixed_gc_centre_fixed>2&diff_waic_gc_v_fixed_gc_both_fixed>2&diff_waic_gc_gc_v_fixed < -2],aes(x=centre,y=logv),col="#009E73",cex=1,pch=15) +
 geom_point(data=gcwAIC[diff_waic_gc_centre_gc_v_fixed>2&diff_waic_gc_centre_gc_both_fixed>2&diff_waic_gc_gc_centre_fixed < -2],aes(x=centre,y=logv),col="#0072B2",cex=1,pch=15) +
 geom_point(data=gcwAIC[diff_waic_gc_gc_centre_fixed>2&diff_waic_gc_gc_v_fixed>2&diff_waic_gc_gc_both_fixed>2],aes(x=centre,y=logv),col="#CC79A7",cex=1,pch=15) +
 geom_hline(yintercept=0,lty=2,col="black") +
 geom_vline(xintercept=0.5,lty=2,col="black") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true centre", y = "", tag = "D")

#p22
###
p33=ggplot()
p33=p33 + geom_point(data=gcwAIC[v_pvalue>0.05&centre_pvalue>0.05],aes(x=centre,y=logv),col="#999999",cex=1,pch=15) +
 geom_point(data=gcwAIC[v_pvalue<=0.05&centre_pvalue>0.05],aes(x=centre,y=logv),col="#0072B2",cex=1,pch=15) +

 geom_point(data=gcwAIC[v_pvalue>0.05&centre_pvalue<=0.05],aes(x=centre,y=logv),col="#009E73",cex=1,pch=15) +
 geom_point(data=gcwAIC[v_pvalue<=0.05&centre_pvalue<=0.05],aes(x=centre,y=logv),col="#CC79A7",cex=1,pch=15) +
 geom_hline(yintercept=0,lty=2,col="black") +
 geom_vline(xintercept=0.5,lty=2,col="black") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true centre", y = "true ln(v)", tag = "E")
###



plot_grid(p1,p11,p2,p22,p33,nrow=3)











