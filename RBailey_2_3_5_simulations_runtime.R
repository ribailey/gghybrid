####################################################################################################
#2.3.5 Code for calculating runtime with different sample sizes#####################################
####################################################################################################

####################################################################################################
#Richard Ian Bailey 23 February 2022################################################################
####################################################################################################
library(gghybrid)
library(ggplot2)
library(cowplot)
####################################################################################################

#Specs for my laptop: HP ZBook Create G7 Notebook PC, with 32Gb RAM and an Intel(R) Core(TM) i9-10885H CPU @ 2.40GHz processor

################################################################################################
#100 haploid individuals with an even distribution of hybrid indices; change the number of loci#
################################################################################################

##########################################
#Regularly spaced hybrid indices##########
##########################################

#To be used for all ggcline runs.

Nind=100

#Make a custom esth object.

hybind=data.table(h_posterior_mode=seq(0,1,length=Nind))#Assume this is a set of test individuals
hybind[,INDLABEL:=paste("A",seq(.N),sep="")]
hindex=list()
hindex$hi=hybind
hindex$hi[,Source:="TEST"]
hindex$test.subject="INDLABEL"



##########################################
##########################################
#Nloci=1##################################
##########################################
##########################################
Nloci=1

simparam=data.table(index=seq(Nloci))
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]

#######################################
#Now I need to create a data.prep file#
#######################################

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]
prepdata[,v:=1][,u:=0][,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]
prepdata[,v:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata[,Source:="TEST"]

prepdata

#################################

set.seed(12345)

gc_L1_ygenomsamp=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    plot.test.subject=c("L1"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,2),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50))

gc_L1_ygenomsamp2=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    nitt = 5000, burnin = 2000, print.k = 50))


Nallelecopies_L1 = prepdata[,.N]

gc_L1_ygenomsamp
#   user  system elapsed 
# 137.63    2.84  142.99 

gc_L1_ygenomsamp2
#   user  system elapsed 
# 130.06    2.11  134.50 



##########################################
##########################################
#Nloci=10##################################
##########################################
##########################################
Nloci=10



simparam=data.table(index=seq(Nloci))
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]


#######################################
#Now I need to create a data.prep file#
#######################################

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]
prepdata[,v:=1][,u:=0][,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]
prepdata[,v:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata[,Source:="TEST"]

prepdata

#################################



set.seed(12345)

gc_L10_ygenomsamp=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    plot.test.subject=c("L10"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,2),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50))

gc_L10_ygenomsamp2=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    nitt = 5000, burnin = 2000, print.k = 50))


Nallelecopies_L10 = prepdata[,.N]


gc_L10_ygenomsamp
#   user  system elapsed 
# 152.09    6.57  162.19

gc_L10_ygenomsamp2
#   user  system elapsed 
# 149.97    6.97  157.14 




##########################################
##########################################
#Nloci=100##################################
##########################################
##########################################
Nloci=100



simparam=data.table(index=seq(Nloci))
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]


#######################################
#Now I need to create a data.prep file#
#######################################

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]
prepdata[,v:=1][,u:=0][,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]
prepdata[,v:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata[,Source:="TEST"]

prepdata

#################################



set.seed(12345)

gc_L100_ygenomsamp=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    plot.test.subject=c("L100"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,2),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50))

gc_L100_ygenomsamp2=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    nitt = 5000, burnin = 2000, print.k = 50))


Nallelecopies_L100 = prepdata[,.N]


gc_L100_ygenomsamp
#   user  system elapsed 
# 271.33   41.19  286.14

gc_L100_ygenomsamp2
#   user  system elapsed 
# 274.03   43.06  285.15




##########################################
##########################################
#Nloci=1000##################################
##########################################
##########################################
Nloci=1000



simparam=data.table(index=seq(Nloci))
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]


#######################################
#Now I need to create a data.prep file#
#######################################

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]
prepdata[,v:=1][,u:=0][,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]
prepdata[,v:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata[,Source:="TEST"]

prepdata

#################################



set.seed(12345)

gc_L1000_ygenomsamp=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    plot.test.subject=c("L1000"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,2),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50))

gc_L1000_ygenomsamp2=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    nitt = 5000, burnin = 2000, print.k = 50))


Nallelecopies_L1000 = prepdata[,.N]


gc_L1000_ygenomsamp
#   user  system elapsed 
#1466.97  224.13 1193.51

gc_L1000_ygenomsamp2
#   user  system elapsed 
#1488.85  203.52 1186.68



##########################################
##########################################
#Nloci=10000##################################
##########################################
##########################################
Nloci=10000



simparam=data.table(index=seq(Nloci))
simparam[,locus:=paste("L",index,sep="")]
simparam[,index:=NULL]


#######################################
#Now I need to create a data.prep file#
#######################################

#1 row per allele copy#
indloc=expand.grid(hybind$INDLABEL,simparam$locus)
indloc=data.table(indloc);setnames(indloc,c("Var1","Var2"),c("INDLABEL","locus"));setkey(indloc,INDLABEL,locus)
indloc[,index:=seq(1,.N)]


setkey(hybind,INDLABEL);setkey(indloc,INDLABEL)
prepdata=hybind[indloc]
prepdata[,v:=1][,u:=0][,S0.prop_1:=0][,S1.prop_1:=1]

##########################################
#Simulate for fixed parental allele freqs#
##########################################

prepdata[,Source_allele:=rbinom(.N,1,prob=S0.prop_1 + (h_posterior_mode^v/(h_posterior_mode^v + (1 - h_posterior_mode)^v*exp(u)))*(S1.prop_1 - S0.prop_1))]
prepdata[,v:=NULL][,h_posterior_mode:=NULL][,u:=NULL]#Saving memory#
prepdata[,Source:="TEST"]

prepdata

#################################



set.seed(12345)

#gc_L10000_ygenomsamp=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    plot.test.subject=c("L10000"),	#Optional
    plot.col=c("orange"),		#Optional
    plot.ylim=c(-1,2),			#Optional
    nitt = 5000, burnin = 2000, print.k = 50))

gc_L10000_ygenomsamp2=system.time(ggcline(
    data.prep.object=prepdata,
    esth.object=hindex,
    read.data.precols=c("INDLABEL"),
    nitt = 5000, burnin = 2000, print.k = 50))


Nallelecopies_L10000 = prepdata[,.N]


#gc_L10000_ygenomsamp#I didn't run this one - decided to only summarize results for runs without plots#


gc_L10000_ygenomsamp2
#    user   system  elapsed 
#14131.67   952.20  9383.67


##############################################################
#Putting it all together######################################
##############################################################


names(gc_L10000_ygenomsamp2)

gc_L10000_ygenomsamp2[1:3]


timeres=data.table(
rbind(
t(data.table(gc_L1_ygenomsamp2)),
t(data.table(gc_L10_ygenomsamp2)),
t(data.table(gc_L100_ygenomsamp2)),
t(data.table(gc_L1000_ygenomsamp2)),
t(data.table(gc_L10000_ygenomsamp2))
),
keep.rownames=T);

timeres[,V4:=NULL];
timeres[,V5:=NULL];

setnames(timeres,c("filename","user.self","sys.self","elapsed"));

timeres[,nloci:=c(1,10,100,1000,10000)]

timeres[,elapsed_min:=elapsed/60]


########
#logged#
########

p1=ggplot()
p1=p1 + geom_point(data=timeres,aes(x=nloci,y=log10(elapsed_min)),col=c("blue","red","green","purple","cyan"))


p1

##########
#unlogged#
##########

p2=ggplot()
p2=p2 + geom_point(data=timeres,aes(x=log10(nloci),y=elapsed_min),col=c("blue","red","green","purple","cyan"))


p2



p3=ggplot()
p3=p3 + geom_point(data=timeres,aes(x=nloci,y=elapsed_min),col=c("blue","red","green","purple","cyan"))


p3

#

p4=ggplot()
p4=p4 + geom_point(data=timeres,aes(x=log10(nloci),y=log10(elapsed_min)),col=c("blue","red","green","purple","cyan"))


p4



plot_grid(p3,p4,p1,p2,nrow=2)


#####################################################
#Calculation of extra seconds per extra locus (ESEL)#
#####################################################

timeres[,timeperlocus_sec:=(elapsed - 134.5)/(nloci - 1)]




















