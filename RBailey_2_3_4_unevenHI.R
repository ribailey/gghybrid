############################################################
#2.3.4 Testing hybrid indexsampling effects#################
############################################################
library(gghybrid)
library(ggplot2)
library(cowplot)
############################################################

############################################################
############################################################
#Choose uneven distribution of hybrid indices###############
############################################################
############################################################

#Making a custom esth object.

hindlabel=list()
hindlabel$test.subject="INDLABEL"

hindlabel$hi=data.table(h_posterior_mode=
c(
rep(0,50),
rep(1,10),
seq(0.1,0.2,length=200),
rep(0.25,100),
rep(0.35,100),
0.4,0.45,0.55,0.7,0.9,
seq(0.91,0.99,length=100)
)
)

hist(hindlabel$hi[,h_posterior_mode],breaks=21)

#Add INDLABEL column#

hindlabel$hi[,index:=seq(1,.N)]
hindlabel$hi[,INDLABEL:=paste("IND",index,sep="")]
hindlabel$hi[,index:=NULL]




###############################################################
#Sample genotypes with evenly distributed centres and v=30#####
###############################################################

Nloci=1000                                    #Number of loci to be simulated#
sim_v=rep(30,Nloci)                            #Cline steepness (I chose the same value for all loci)#
sim_centre=seq(0.0001,0.9999,length=Nloci)    #Cline centre (I chose evenly spaced centres)#
###############################################################

simparam=data.table(sim_v=sim_v,sim_centre=sim_centre)

#Histogram of cline centres.

hist(simparam$sim_centre, breaks=21)


##############################################
#All loci fixed in parental reference samples#
##############################################

#This could be altered, but I think it would make the results harder to interpret.

simparam[,S0_af1.2:=0][,S1_af1.2:=1]#Parental reference freqency of S1 allele#


#Add locus names#

simparam[,index:=seq(1,.N)]
simparam[,locus:=paste("L",index,sep="")]

simparam

############################################################
############################################################
#Now simulate the individual genotypes, assuming no linkage#
############################################################
############################################################

#######################################
#Now I need to create a data.prep file#
#######################################


simparam[,index:=NULL]#Only used to create the locus name#

#Make a table with 1 row per allele copy for a diploid set of loci.

indloc=expand.grid(hindlabel$hi$INDLABEL,simparam$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))#This makes them all diploid#
indloc[,index:=seq(1,.N)]

indloc


#Join the tables.

setkey(hindlabel$hi,INDLABEL);setkey(simparam,locus);setkey(indloc,INDLABEL)
prepdata=hindlabel$hi[indloc]

setkey(prepdata,locus)
prepdata=simparam[prepdata]

prepdata[,Source:="TEST"]


########################################################
#Now simulate genotypes using the genomic cline formula#
########################################################

prepdata[,sim_u:=qlogis(sim_centre)*sim_v]

prepdata[,Source_allele:=rbinom(.N,1,prob=S0_af1.2 + (h_posterior_mode^sim_v/(h_posterior_mode^sim_v + (1 - h_posterior_mode)^sim_v*exp(sim_u)))*(S1_af1.2 - S0_af1.2))]

prepdata[is.na(Source_allele)]#No missing, so it worked. Source_allele is the name of the column containing the allele IDs used in esth and ggcline#

setnames(prepdata,c("S0_af1.2","S1_af1.2"),c("S0.prop_1","S1.prop_1"))#I should have named them properly at the start, but whatever#

prepdata=prepdata[,.(S0.prop_1,S1.prop_1,locus,Source,INDLABEL,POPID,index,Source_allele)]#Saving memory#


prepdata



###########################################
###########################################
#Run genomic clines########################
###########################################
###########################################

#First on the complete sample##############


gc_uneven=ggcline(
  data.prep.object=prepdata,
  esth.object=hindlabel,
  #include.Source=TRUE,                    #All are TEST#
  read.data.precols=c("INDLABEL"),
  plot.test.subject = c("L200","L600"),      #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-10,40),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  nitt=5000,
  burnin=2000
)

gc_uneven




##################################################
##################################################
#Repeat everything with v=1#######################
##################################################
##################################################

hindlabel=list()
hindlabel$test.subject="INDLABEL"

hindlabel$hi=data.table(h_posterior_mode=
c(
rep(0,50),
rep(1,10),
seq(0.1,0.2,length=200),
rep(0.25,100),
rep(0.35,100),
0.4,0.45,0.55,0.7,0.9,
seq(0.91,0.99,length=100)
)
)

hist(hindlabel$hi[,h_posterior_mode],breaks=21)

#Add INDLABEL column#

hindlabel$hi[,index:=seq(1,.N)]
hindlabel$hi[,INDLABEL:=paste("IND",index,sep="")]
hindlabel$hi[,index:=NULL]



###############################################################
#Sample genotypes with evenly distributed centres and v=1######
###############################################################


Nloci=1000                                    #Number of loci to be simulated#
sim_v=rep(1,Nloci)                            #Cline steepness (I chose the same value for all loci)#
sim_centre=seq(0.0001,0.9999,length=Nloci)    #Cline centre (I chose evenly spaced centres)#
###############################################################

simparam=data.table(sim_v=sim_v,sim_centre=sim_centre)

#Histogram of cline centres.

hist(simparam$sim_centre, breaks=21)


##############################################
#All loci fixed in parental reference samples#
##############################################

#This could be altered, but I think it would make the results harder to interpret.

simparam[,S0_af1.2:=0][,S1_af1.2:=1]#Parental reference freqency of S1 allele#


#Add locus names#

simparam[,index:=seq(1,.N)]
simparam[,locus:=paste("L",index,sep="")]

simparam

############################################################
############################################################
#Now simulate the individual genotypes, assuming no linkage#
############################################################
############################################################

#######################################
#Now I need to create a data.prep file#
#######################################


simparam[,index:=NULL]#Only used to create the locus name#

#Make a table with 1 row per allele copy for a diploid set of loci.

indloc=expand.grid(hindlabel$hi$INDLABEL,simparam$locus)
head(indloc)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc=rbind(indloc,copy(indloc))#This makes them all diploid#
indloc[,index:=seq(1,.N)]

indloc


#Join the tables.

setkey(hindlabel$hi,INDLABEL);setkey(simparam,locus);setkey(indloc,INDLABEL)
prepdata=hindlabel$hi[indloc]

setkey(prepdata,locus)
prepdata=simparam[prepdata]

prepdata[,Source:="TEST"]


########################################################
#Now simulate genotypes using the genomic cline formula#
########################################################

prepdata[,sim_u:=qlogis(sim_centre)*sim_v]

prepdata[,Source_allele:=rbinom(.N,1,prob=S0_af1.2 + (h_posterior_mode^sim_v/(h_posterior_mode^sim_v + (1 - h_posterior_mode)^sim_v*exp(sim_u)))*(S1_af1.2 - S0_af1.2))]

prepdata[is.na(Source_allele)]#No missing, so it worked. Source_allele is the name of the column containing the allele IDs used in esth and ggcline#

setnames(prepdata,c("S0_af1.2","S1_af1.2"),c("S0.prop_1","S1.prop_1"))#I should have named them properly at the start, but whatever#

prepdata=prepdata[,.(S0.prop_1,S1.prop_1,locus,Source,INDLABEL,index,Source_allele)]#Saving memory#


prepdata



###########################################
###########################################
#Run genomic clines########################
###########################################
###########################################

#First on the complete sample##############


gc_uneven2=ggcline(
  data.prep.object=prepdata,
  esth.object=hindlabel,
  #include.Source=TRUE,                    #All are TEST#
  read.data.precols=c("INDLABEL"),
  plot.test.subject = c("L200","L600"),      #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("orange","cyan"),          #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-5,5),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  nitt=5000,
  burnin=2000
)

gc_uneven2



#######################################################################
#######################################################################
#Now join the results to the true centre and v values##################
#######################################################################
#######################################################################


setkey(simparam,locus);setkey(gc_uneven2$gc,locus)
res_all=simparam[gc_uneven2$gc]


###
#v#
###


###
p1=ggplot(data=res_all,aes(x=sim_centre,y=mean_log_v))

p1 =p1 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") +
geom_errorbar(aes(ymin=log(v_lower_95), ymax=log(v_upper_95)),width=0,col="grey") +
geom_point(cex=1) +
#geom_hline(yintercept=log(30),col="red",lty=2) +
geom_hline(yintercept=log(1),col="blue",lty=2) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true cline centre", y = "estimated ln(v)", tag = "A")

p1
###


########
#centre#
########

###
p3=ggplot(data=res_all,aes(x=sim_centre,y=invlogit_mean_logit_centre))

p3=p3 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") +
geom_errorbar(aes(ymin=centre_lower_95, ymax=centre_upper_95),width=0,col="grey") +
geom_point(cex=1) +
geom_abline(intercept=0,slope=1,col="red",lty=2) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true cline centre", y = "estimated cline centre", tag = "B")

p3
###





##############################################################################
#Now the same for v=30########################################################
##############################################################################


setkey(simparam,locus);setkey(gc_uneven$gc,locus)
res_all=simparam[gc_uneven$gc]

fwrite(res_all,"gc_uneven_v30.csv")

#Now back in v1 workspace#

res_all30=fread("gc_uneven_v30.csv")



###
#v#
###


###
p12=ggplot(data=res_all30,aes(x=sim_centre,y=mean_log_v))

p12 =p12 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") +
geom_errorbar(aes(ymin=log(v_lower_95), ymax=log(v_upper_95)),width=0,col="grey") +
geom_point(cex=1) +
geom_hline(yintercept=log(30),col="red",lty=2) +
geom_hline(yintercept=log(1),col="blue",lty=2) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true cline centre", y = "estimated ln(v)", tag = "C")

p12
###


########
#centre#
########

###
p32=ggplot(data=res_all30,aes(x=sim_centre,y=invlogit_mean_logit_centre))

p32=p32 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") +
geom_errorbar(aes(ymin=centre_lower_95, ymax=centre_upper_95),width=0,col="grey") +
geom_point(cex=1) +
geom_abline(intercept=0,slope=1,col="red",lty=2) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true cline centre", y = "estimated cline centre", tag = "B")

p32
###


plot_grid(p1,p3,p12,p32,nrow=2)

plot_grid(p12,p32,nrow=2)









############################################################################################################
#Next, use the threshold criteria for choosing interesting loci with respect to v, and plot the histogram###
############################################################################################################

#All

res_all30[order(v_pvalue)]#Just choose some arbitrary p value - I'm struggling with error messages from qvalue#


#Make the bins#

###
res_all30[,histfac:=cut(invlogit_mean_logit_centre,breaks=41)];#Was 'gc.centre_mean'#

res_all30[,histfacL:=tstrsplit(histfac,"(",fixed=TRUE, keep=2L)];#Lower breakpoint only#
res_all30[,histfacL:=tstrsplit(histfacL,",",fixed=T,keep=1L)];#Lower breakpoint only#

res_all30[,histfacU:=tstrsplit(histfac,",",fixed=TRUE, keep=2L)];#Upper breakpoint only#
res_all30[,histfacU:=tstrsplit(histfacU,"]",fixed=T,keep=1L)];#Upper breakpoint only#

res_all30[,histfacL:=as.numeric(histfacL)]
res_all30[,histfacU:=as.numeric(histfacU)]

res_all30[,histfacmid:=histfacL + (histfacU - histfacL)/2]#Midpoint of each bin#
###


#True numbers of loci v>20 per bin#NOTE THAT THESE BINS ARE NOT IDENTICAL TO THE BINS BASED ON ESTIMATED CLINE CENTRE#

hist(res_all30$sim_centre, breaks=41)#With 41 breaks, 20 loci per bin#


#Now the estimated numbers#

vcutoff=20
pcutoff=0.05


vN_all=res_all30[v_pvalue < pcutoff &exp_mean_log_v>vcutoff,.N,by=histfacmid]

vN_all[,sum(N)]#Total number passing the criteria#


p7=ggplot(data=vN_all,aes(x=histfacmid,y=N))
p7 = p7 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") + 
 geom_col(fill="blue") +
 geom_hline(yintercept=25,col="red",lty=2) +  #NOT SURE IF 20 IS CORRECT HERE#
 ylim(0,100) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="estimated cline centre binned", y = "N loci above threshold", tag = "B")

 #labs(title="N significantly steep loci with v > 20 by cline centre", x="Cline centre", y = "N loci")


p7

#So a lot of potential biases here!!







#################################
#Same plots, against true centre#
#################################

###
res_all30[,histfac:=cut(sim_centre,breaks=41)];#Was 'gc.centre_mean'#

res_all30[,histfacL:=tstrsplit(histfac,"(",fixed=TRUE, keep=2L)];#Lower breakpoint only#
res_all30[,histfacL:=tstrsplit(histfacL,",",fixed=T,keep=1L)];#Lower breakpoint only#

res_all30[,histfacU:=tstrsplit(histfac,",",fixed=TRUE, keep=2L)];#Upper breakpoint only#
res_all30[,histfacU:=tstrsplit(histfacU,"]",fixed=T,keep=1L)];#Upper breakpoint only#

res_all30[,histfacL:=as.numeric(histfacL)]
res_all30[,histfacU:=as.numeric(histfacU)]

res_all30[,histfacmid:=histfacL + (histfacU - histfacL)/2]#Midpoint of each bin#
###

vcutoff=20
pcutoff=0.05


vN_all=res_all30[v_pvalue < pcutoff &exp_mean_log_v>vcutoff,.N,by=histfacmid]

vN_all[,sum(N)]#Total number passing the criteria#


p72=ggplot(data=vN_all,aes(x=histfacmid,y=N))
p72 = p72 + geom_vline(xintercept=hindlabel$hi$h_posterior_mode,col="lightblue") + 
 geom_col(fill="blue") +
 geom_hline(yintercept=25,col="red",lty=2) +
 ylim(0,100) +
theme_classic() +
 theme(legend.position="none",axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
 labs(x="true cline centre binned", y = "N loci above threshold", tag = "A")

 #labs(title="N significantly steep loci with v > 20 by cline centre", x="Cline centre", y = "N loci")




plot_grid(p72,p7,nrow=2)



#################################################################################################################################################
#However, it should hopefully be okay to compare subsets (e.g. fertility genes) to the global number per bin, to identify overrepresented groups#
#################################################################################################################################################














