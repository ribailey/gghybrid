################################################################################
#Statistical properties of ggcline tests########################################
################################################################################

################################################################################
#Richard Bailey 19 Feb 2022#####################################################
################################################################################
library(gghybrid)
library(ggplot2)
library(cowplot)
################################################################################


#####################################################################################################################################
#####################################################################################################################################
#2.3.3b Comparison of waic versus AIC versus p value for testing for true deviations from null. 
#Make new simulated datasets with either v or centre having null true value, and compare 1 parameter fixed to both parameters fixed.
#####################################################################################################################################
#####################################################################################################################################

#Given that we are only varying, and only testing, one parameter in each dataset, these results can be used to compare AIC, waic 
#and p value for (a) choosing loci based on a threshold and (b) ranking loci according to deviation from null.

#Make a custom esth object.

hybind=data.table(h_posterior_mode=seq(0,1,length=1000))#Assume this is a set of test individuals
hybind[,INDLABEL:=paste("A",seq(.N),sep="")]
hindex=list()
hindex$hi=hybind
hindex$hi[,Source:="TEST"]
hindex$test.subject="INDLABEL"


#Fixed parental allele frq#

#Choose evenly spaced v values on the log scale#

#Choose evenly spaced centre on the data scale#0.000001 to 0.999999 to match those depicted in Figure 1#


#Include the null value for each parameter#

vs=rep(seq(-4,4,length=49),20)#This goes well beyond the range depicted in figure 1 a-c#
expvs=exp(vs)
centres=rep(seq(0.000001,0.999999,length=49),20)

####################################################################################################################################################################################
#Create SNP datasets with either centre=0.5 (null value) across all loci and v values created above (simparam1) or v=1 (null value) across all loci and centre values created above#
####################################################################################################################################################################################

simparam1=data.table(v=expvs,centre=0.5)
simparam2=data.table(v=1,centre=centres)
###

simparam1[,S0_af1.2:=0][,S1_af1.2:=1]
simparam1[,index:=seq(1,.N)]
simparam1[,locus:=paste("L",index,sep="")][,index:=NULL]

simparam2[,S0_af1.2:=0][,S1_af1.2:=1]
simparam2[,index:=seq(1,.N)]
simparam2[,locus:=paste("L",index,sep="")][,index:=NULL]

#Diploid loci, 1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam1$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))
indloc[,index:=seq(1,.N)]

setkey(hybind,INDLABEL);setkey(simparam1,locus);setkey(indloc,INDLABEL)
prepdata4=hybind[indloc]
setkey(prepdata4,locus)
prepdata4=simparam1[prepdata4]

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam2$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))
indloc[,index:=seq(1,.N)]

setkey(hybind,INDLABEL);setkey(simparam2,locus);setkey(indloc,INDLABEL)
prepdata5=hybind[indloc]
setkey(prepdata5,locus)
prepdata5=simparam2[prepdata5]


##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata4[,u:=qlogis(centre)*v]
prepdata4[,Source_allele:=rbinom(.N,1,prob=S0_af1.2 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1_af1.2 - S0_af1.2))]
prepdata4[is.na(Source_allele)]#No missing#

setnames(prepdata4,c("S0_af1.2","S1_af1.2"),c("S0.prop_1","S1.prop_1"))

prepdata4[,v:=NULL][,centre:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata4[,Source:="TEST"]
###
prepdata5[,u:=qlogis(centre)*v]
prepdata5[,Source_allele:=rbinom(.N,1,prob=S0_af1.2 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1_af1.2 - S0_af1.2))]
prepdata5[is.na(Source_allele)]#No missing#

setnames(prepdata5,c("S0_af1.2","S1_af1.2"),c("S0.prop_1","S1.prop_1"))

prepdata5[,v:=NULL][,centre:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata5[,Source:="TEST"]

###


###############################################################################################################################
#Run 1-parameter versus no-parameter models for each dataset, in each case fixing the other parameter to the (true) null value#
###############################################################################################################################


locuslist=simparam1$locus
###
gc4=ggcline(
  data.prep.object=prepdata4,              #Needs an entry#
  esth.object=hindex,                      #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                  #***NEED THIS TRUE***#
  read.data.precols="INDLABEL",            #Needs an entry#
  #fix.subject.v = FALSE,
  #fix.value.v,
  fix.subject.centre = locuslist,          #***FIX CENTRE FOR ALL LOCI***#
  fix.value.centre=0.5,                    #***FIXING CENTRE TO THE TRUE VALUE***#
  plot.test.subject = c("L1","L25"),       #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),           #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-0.5, 1.5),                    #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),             #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=5000,                               #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                             #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
###
gc4
##########################################################################

simparam2[locus%in%c("L1","L25")]


locuslist=simparam2$locus
###
gc5=ggcline(
  data.prep.object=prepdata5,              #Needs an entry#
  esth.object=hindex,                      #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                  #***NEED THIS TRUE***#
  read.data.precols="INDLABEL",            #Needs an entry#
  fix.subject.v = locuslist,               #***FIX V FOR ALL LOCI***#
  fix.value.v = 1,                         #***FIXING V TO THE TRUE VALUE***#
  fix.subject.centre = FALSE,          
  fix.value.centre,                    
  plot.test.subject = c("L1","L25"),       #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),           #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-11, 1),                    #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),             #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=5000,                               #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=2000,                             #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
###
gc5
##########################################################################


######################################
#Now the null model for both datasets#
######################################

locuslist=simparam1$locus
###
gc4.1=ggcline(
  data.prep.object=prepdata4,              #Needs an entry#
  esth.object=hindex,                      #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                  #***NEED THIS TRUE***#
  read.data.precols="INDLABEL",            #Needs an entry#
  fix.subject.v = locuslist,               #***FIX V FOR ALL LOCI***#
  fix.value.v = 1,                         #***FIXING V TO THE TRUE VALUE***#
  fix.subject.centre = locuslist,          #***FIX CENTRE FOR ALL LOCI***#
  fix.value.centre=0.5,                    #***FIXING CENTRE TO THE TRUE VALUE***#
  #plot.test.subject = c("L1","L15"),       #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),           #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.ylim = c(-0.5, 1.5),                    #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),             #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=2,                               #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=0,                             #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
###


locuslist=simparam2$locus
###
gc5.1=ggcline(
  data.prep.object=prepdata5,              #Needs an entry#
  esth.object=hindex,                      #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = FALSE,                  #Default is FALSE#
  return.likmeans = TRUE,                  #***NEED THIS TRUE***#
  read.data.precols="INDLABEL",            #Needs an entry#
  fix.subject.v = locuslist,               #***FIX V FOR ALL LOCI***#
  fix.value.v = 1,                         #***FIXING V TO THE TRUE VALUE***#
  fix.subject.centre = locuslist,          #***FIX CENTRE FOR ALL LOCI***#
  fix.value.centre=0.5,                    #***FIXING CENTRE TO THE TRUE VALUE***#
  #plot.test.subject = c("L1","L15"),       #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.col = c("orange","cyan"),           #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.ylim = c(-0.5, 1.5),                    #Remove plots to speed up runtime, they are just to check the MCMC is working#
  #plot.pch.v.centre = c(1, 3),             #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre=c(0,sqrt(50)),         #Preferred default#
  prior.logv=c(0,sqrt(10)),                #Preferred default#
  nitt=2,                               #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=0,                             #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 50
)
###


#########################
#Compare models##########
#########################

comp_centre_fixed=compare.models(
 ggcline.object1=gc4,
 ggcline.object2=gc4.1,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_centre_fixed


comp_v_fixed=compare.models(
 ggcline.object1=gc5,
 ggcline.object2=gc5.1,
 ggcline.pooled=FALSE              #No pooling so use the default FALSE#
)

comp_v_fixed

calc1=calc_AIC(
 data.prep.object=prepdata5,
 esth.object=hindex,
 ggcline.object=gc5,
 test.subject="locus")

calc2=calc_AIC(
 data.prep.object=prepdata5,
 esth.object=hindex,
 ggcline.object=gc5.1,
 test.subject="locus")

####

calc4=calc_AIC(
 data.prep.object=prepdata4,
 esth.object=hindex,
 ggcline.object=gc4,
 test.subject="locus")

calc5=calc_AIC(
 data.prep.object=prepdata4,
 esth.object=hindex,
 ggcline.object=gc4.1,
 test.subject="locus")


########################################
########################################
#Plotting waic diff against true values#
########################################
########################################


#p values against logv and centre

#waicdiff against logv and centre

#AICdiff against logv and centre

#Direct comparison of all three




######################
######################
######################

setkey(simparam1,locus);setkey(gc4$gc,locus)
gcsim1=simparam1[gc4$gc];gcsim1[,logv:=log(v)]
setkey(comp_centre_fixed,locus);setkey(gcsim1,locus);
gcsim1=comp_centre_fixed[gcsim1];


setkey(gcsim1,locus);setkey(calc4,locus);setkey(calc5,locus)
gcsim1=gcsim1[calc4][calc5]
#
setkey(simparam2,locus);setkey(gc5$gc,locus)
gcsim2=simparam2[gc5$gc];
setkey(comp_v_fixed,locus);setkey(gcsim2,locus);
gcsim2=comp_v_fixed[gcsim2];

setkey(gcsim2,locus);setkey(calc1,locus);setkey(calc2,locus)
gcsim2=gcsim2[calc1][calc2]
#

###
par(mfrow=c(2,3))
###
plot(gcsim1[,-log10(v_pvalue+1e-200)]~gcsim1[,logv])#true values#
abline(h=-log10(0.05),v=0,col="red",lty=2)
###
plot(gcsim1[,waicdiff_npar_AICscale]~gcsim1[,logv])#true values#
#points(gcsim1[,AIC - i.AIC]~gcsim1[,logv],col="orange",pch=3)
abline(h=2,v=0,col="red",lty=2)
###
plot(gcsim1[,waicdiff_npar_AICscale]~gcsim1[,AIC - i.AIC])
abline(a=0,b=1,col="red",lty=2)
###
plot(gcsim2[,-log10(centre_pvalue+1e-200)]~gcsim2[,centre])#the added value is because one p value was 0#
abline(h=-log10(0.05),v=0.5,col="red",lty=2)
###
plot(gcsim2[,waicdiff_npar_AICscale]~gcsim2[,centre])
#points(gcsim2[,AIC - i.AIC]~gcsim2[,centre],col="orange",pch=3)
abline(h=2,v=0.5,col="red",lty=2)
###
plot(gcsim2[,waicdiff_npar_AICscale]~gcsim2[,AIC - i.AIC])
abline(a=0,b=1,col="red",lty=2)
###

gcsim2[order(AIC)]

######################
######################
######################
plot(gcsim1[,-log10(v_pvalue)]~gcsim1[,waicdiff_npar_AICscale])
points(gcsim1[,-log10(v_pvalue)]~gcsim1[,AIC - i.AIC],col="orange",pch=3)
abline(h=-log10(0.05),v=2,col="red",lty=2)

plot(gcsim2[,-log10(centre_pvalue+1e-200)]~gcsim2[,waicdiff_npar_AICscale])
points(gcsim2[,-log10(centre_pvalue+1e-200)]~gcsim2[,AIC - i.AIC],col="orange",pch=3)
abline(h=-log10(0.05),v=2,col="red",lty=2)
#####################################################################################################
#Conclusion: waic model comparison is way better behaved than the p value!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#####################################################################################################

#Recommend that p value or FDR can be used as a cutoff, but not for ranking loci. Best option for ranking is AIC/waic, 
#but if only running the full model, use a p value cutoff then rank by parameter estimates.

###########################################
###########################################
#ggplot version############################
###########################################
###########################################


###
p1=ggplot()
p1=p1 + geom_point(data=gcsim1[v_pvalue <=0.05],aes(x=logv,y=-log10(v_pvalue+1e-200)),col="blue",cex=0.5) +
 geom_point(data=gcsim1[v_pvalue >0.05],aes(x=logv,y=-log10(v_pvalue+1e-200)),col="red",cex=0.5) +
 geom_hline(yintercept=-log10(0.05),lty=2,col="grey") +
 geom_vline(xintercept=0,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="ln(true v)", y = "-log10 (v p value + 1e-200)", tag = "A")

p1
###
p2=ggplot()
p2=p2 + geom_point(data=gcsim2[centre_pvalue <=0.05],aes(x=centre,y=-log10(centre_pvalue+1e-200)),col="blue",cex=0.5) +
 geom_point(data=gcsim2[centre_pvalue >0.05],aes(x=centre,y=-log10(centre_pvalue+1e-200)),col="red",cex=0.5) +
 geom_hline(yintercept=-log10(0.05),lty=2,col="grey") +
 geom_vline(xintercept=0.5,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true centre", y = "-log10 (centre p value + 1e-200)", tag = "B")

p2
###
p3=ggplot()
p3=p3 + geom_point(data=gcsim1[-waicdiff_npar_AICscale >= 2],aes(x=logv,y=-waicdiff_npar_AICscale),col="blue",cex=0.5) +
 geom_point(data=gcsim1[-waicdiff_npar_AICscale < 2],aes(x=logv,y=-waicdiff_npar_AICscale),col="red",cex=0.5) +
 geom_hline(yintercept=2,lty=2,col="grey") +
 geom_vline(xintercept=0,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="ln(true v)", y = "waic difference", tag = "C")

p3
###
p4=ggplot()
p4=p4 + geom_point(data=gcsim2[-waicdiff_npar_AICscale >= 2],aes(x=centre,y=-waicdiff_npar_AICscale),col="blue",cex=0.5) +
 geom_point(data=gcsim2[-waicdiff_npar_AICscale < 2],aes(x=centre,y=-waicdiff_npar_AICscale),col="red",cex=0.5) +
 geom_hline(yintercept=2,lty=2,col="grey") +
 geom_vline(xintercept=0.5,lty=2,col="grey") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true centre", y = "waic difference", tag = "D")

p4
###
p5=ggplot()
p5=p5 + geom_point(data=gcsim1,aes(x=(i.AIC - AIC),y=-waicdiff_npar_AICscale),cex=0.5) +
 geom_abline(intercept=0,slope=1,lty=2,col="red") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="AIC difference", y = "waic difference", tag = "E")

p5
###
p6=ggplot()
p6=p6 + geom_point(data=gcsim2,aes(x=(i.AIC - AIC),y=-waicdiff_npar_AICscale),cex=0.5) +
 geom_abline(intercept=0,slope=1,lty=2,col="red") +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="AIC difference", y = "waic difference", tag = "F")

p6
###


plot_grid(p1,p2,p3,p4,p5,p6,nrow=3)










