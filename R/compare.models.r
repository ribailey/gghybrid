#' Statistical comparison of pairs of either \code{esth} or \code{ggcline} objects with different pooling or fixing of parameters, using the 
#'   Bayesian widely applicable information criterion (waic).
#'
#' @param ggcline.object1 Name of a \code{ggcline} object.
#' @param ggcline.object2 Name of a \code{ggcline} object.
#' @param ggcline.pooled Logical. If \code{FALSE}, each row of the \code{ggcline$gc} objects is compared separately. If 
#'   \code{TRUE}, a single comparison is made across all rows (to be used when one or both parameters is pooled across 
#'   all samples in at least one of the two models being compared).
#' @param esth.object1 Name of an \code{esth} object.
#' @param esth.object2 Name of an \code{esth} object.
#' @details \code{compare.models} While the calculation of the widely applicable information criterion (waic) normally 
#'   involves an estimated effective number of parameters, testing in \code{gghybrid} suggests that this value is 
#'   systematically overestimated. \code{compare.models} presents the estimated effective number of parameters and the 
#'   resulting waic and difference in waic, in columns including \sQuote{effpar} in the heading, but always produces a 
#'   warning to instead use the equivalent \sQuote{npar} columns, which use the actual number of estimated parameters 
#'   in the calculations. This produces results closely comparable with AIC.
#' 
#' Unlike the \code{calc_Aic} function, which calculates AIC for a single model, \code{compare.models} is set up to 
#'   compare two models, and therefore requires two \code{esth} or two \code{ggcline} objects as input.
#' 
#' waic is calculated as \code{-2*lppd + 2*npar}, and is therefore on the same scale as AIC, for direct comparison.
#' @return A \code{data.table} and \code{data.frame} with the name of the test subject (when there is more than one) plus following fields:
#'   \item{lppd1_sum}{Sum of the mean posterior likelihoods for the individual data points for the first model.}
#'   \item{effpar1_sum}{Sum of the effective number of parameters for the individual data points for the first model.}
#'   \item{npar1}{True number of parameters for the first model.}
#'   \item{waic1_npar1_AICscale}{widely applicable information criterion for the first model using the true number of parameters.}
#'   \item{waic1_effpar1_AICscale}{widely applicable information criterion for the first model using the effective number of parameters.}
#'   \item{lppd2_sum}{Sum of the mean posterior likelihoods for the individual data points for the second model.}
#'   \item{effpar1_sum}{Sum of the effective number of parameters for the individual data points for the second model.}
#'   \item{npar2}{number of parameters for the second model.}
#'   \item{waic2_npar2_AICscale}{widely applicable information criterion for the second model using the true number of parameters.}
#'   \item{waic2_effpar2_AICscale}{widely applicable information criterion for the second model using the effective number of parameters.}
#'   \item{waicdiff_npar_AICscale}{Difference in waic between the two models using the true number of parameters, with negative values favouring the more complex model.}
#'   \item{waicdiff_effpar_AICscale}{Difference in waic between the two models using the effective number of parameters, with negative values favouring the more complex model.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' ####################################################################################################
#' ####################################################################################################
#' #Comparison of hybrid index estimated at two different resolutions, individual and population level#
#' ####################################################################################################
#' ####################################################################################################
#' 
#' #Read and prepare the sparrow data, and run esth with INDLABEL and POPID respectively as test.subject.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);
#' 
#' #Get the individual level hybrid index estimates.
#' hindlabel=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  test.subject = "INDLABEL",                  #default#
#'  test.subject.compare = "POPID",             #***the test.subject you want to compare with the current test.subject***#
#'  include.Source=TRUE,
#'  return.likmeans=TRUE,                       #Must be set to TRUE before using compare.models#
#'  nitt=3000,
#'  burnin=1000
#'  )
#' 
#' #Run esth again, with test.subject="POPID". No need to declare 'test.subject.compare' when running the lower resolution model.
#' hindlabel2=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  test.subject = "POPID",                     #***Changed to "POPID"***#
#'  include.Source=TRUE,
#'  return.likmeans=TRUE,                       #Must be set to TRUE before using compare.models#
#'  nitt=3000,
#'  burnin=1000
#'  )
#' 
#' #Compare waic between the two models for each POPID.
#' comp1=compare.models(esth.object1=hindlabel,esth.object2=hindlabel2);
#' 
#' #Order by the difference in waic between the two models. Negative values for 'waicdiff_npar_AICscale' favour the individual level model.
#' comp1[order(waicdiff_npar_AICscale)]
#' 
#' ################################
#' ################################
#' #Genomic cline model comparison#
#' ################################
#' ################################
#' 
#' ##############################################################
#' #Parameter v pooled versus unpooled across a whole chromosome#
#' ##############################################################
#' 
#' #Calculate a single waic for the sparrow Z chromosome for two ggcline models run in the examples under ggcline, 
#' #one with v pooled and one unpooled.
#' 
#' #A negative value for 'waicdiff_npar_AICscale' indicates support for the more complex model, with v estimated independently per locus.
#' gclinecomp=compare.models(
#'  ggcline.object1=gc_unpooled_z,
#'  ggcline.object2=gc_poolv_Z,
#'  ggcline.pooled=TRUE              #Default is FALSE. This must be set to TRUE when one of the ggcline runs includes pooling#
#' )
#' 
#' #############################################################################################################
#' #At the individual locus level, compare a model with v and centre estimated versus a model with v fixed at 1#
#' #############################################################################################################
#' 
#' #Compare the ggcline objects gc (no fixing) and gc_v_fixed created in the ggcline examples.
#' 
#' gclinecomp_vfixed=compare.models(
#'  ggcline.object1=gc_v_fixed,
#'  ggcline.object2=gc,
#'  ggcline.pooled=FALSE              #No pooling so use the default FALSE#
#' )
#' 
#' #This time we get a separate estimate of difference in waic per locus. Negative values for 'waicdiff_npar_AICscale' support deviation from v=1.
#' gclinecomp_vfixed[order(waicdiff_npar_AICscale)]
#' 
#' #########################################################
#' #Comparing differing resolutions at the lower resolution#
#' #########################################################
#' 
#' #Compare the ggcline objects gc (1 cline per locus) and gc_by_chr (1 cline per chromosome) created in the ggcline examples.
#' 
#' #If the column "chrom" is absent from the 'gc' results object and the prepdata$data.prep object, it must first be added.
#' chrom=fread("markerchr.csv")#This contains only the loci retained in prepdata$data.prep#
#' 
#' setkey(prepdata$data.prep,locus);setkey(chrom,locus);
#' prepdata$data.prep=prepdata$data.prep[chrom]
#' 
#' setkey(gc$gc,locus);setkey(chrom,locus);setkey(gc$likmeans,locus);
#' gc$gc=gc$gc[chrom];
#' gc$likmeans=gc$likmeans[chrom]
#' 
#' #This produces one model comparison per chromosome.
#' gclinecomp_chr=compare.models(
#'  ggcline.object1=gc,
#'  ggcline.object2=gc_by_chr,
#'  ggcline.pooled=FALSE              #No pooling so use the default FALSE#
#' )
#' 
#' #Again, negative values for 'waicdiff_npar_AICscale' favour the more complex model (independent cline parameters per 
#' #locus within each chromosome).
#' gclinecomp_chr[order(waicdiff_npar_AICscale)]
#' }
#' @export
compare.models=function(ggcline.object1=NULL,ggcline.object2=NULL,
    ggcline.pooled=FALSE,
    esth.object1=NULL,esth.object2=NULL){

    if(is.null(ggcline.object1)==F){
        if(ggcline.pooled){
            npar1=unique(ggcline.object1$gc[,sum_npar]);
            MET.waic1=ggcline.object1$likmeans;

            MET.1=data.table(lppd1_sum=MET.waic1[,sum(log(postmean.exp.loglik), na.rm=T)],
			    effpar1_sum=MET.waic1[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T)]);#added 14 Dec 2021. Using sum rather than mean so is on the same scale as lppd1 and AIC#
            MET.1[,
                npar1:=npar1][,
                waic1_npar1_AICscale:=-2*lppd1_sum + 2*npar1][,
                waic1_effpar1_AICscale:=-2*lppd1_sum + 2*effpar1_sum];

            npar2=unique(ggcline.object2$gc[,sum_npar]);
            MET.waic2=ggcline.object2$likmeans;

            MET.2=data.table(lppd2_sum=MET.waic2[,sum(log(postmean.exp.loglik), na.rm=T)],
			    effpar2_sum=MET.waic2[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T)]);#added 14 Dec 2021. Using sum rather than mean so is on the same scale as lppd1 and AIC#
            MET.2[,
                npar2:=npar2][,
                waic2_npar2_AICscale:=-2*lppd2_sum + 2*npar2][,
                waic2_effpar2_AICscale:=-2*lppd2_sum + 2*effpar2_sum];

#setkey(MET.waic1,index);setkey(MET.waic2,index);
#MET.waic=MET.waic1[MET.waic2];
#se_diff_mean=MET.waic[,sqrt(var(log(postmean.exp.loglik) - log(i.postmean.exp.loglik)))];

            MET.both=cbind(MET.1,MET.2);
            MET.both[npar1 >= npar2,
			    waicdiff_npar_AICscale:=waic1_npar1_AICscale - waic2_npar2_AICscale][npar1 < npar2,
			    waicdiff_npar_AICscale:=waic2_npar2_AICscale - waic1_npar1_AICscale][effpar1_sum >= effpar2_sum,
				waicdiff_effpar_AICscale:=waic1_effpar1_AICscale - waic2_effpar2_AICscale][effpar1_sum < effpar2_sum,
			    waicdiff_effpar_AICscale:=waic2_effpar2_AICscale - waic1_effpar1_AICscale];
			
warning("While the effective number of parameters ('effpar') is calculated and presented, and used in some of the waic calculations, testing suggests it is systematically overestimated and should not be trusted. Unless the estimated effpar values are the same as or less than the npar values (the actual number of parameters), use the calculations based on npar");
            return(MET.both);
                          }else{ 
			
			if(unique(ggcline.object1$gc$sum_npar) > unique(ggcline.object2$gc$sum_npar)){

#setkey(ggcline.object1$likmeans,index);setkey(ggcline.object2$likmeans,index);
#MET.waic=ggcline.object1$likmeans[ggcline.object2$likmeans];
#se_diff_mean=MET.waic[,sqrt(var(log(postmean.exp.loglik) - log(i.postmean.exp.loglik))),by=c(ggcline.object2$test.subject)];
#setkeyv(se_diff_mean,ggcline.object2$test.subject);
#setnames(se_diff_mean,"V1","se_diff_mean");
			
			    npar1=ggcline.object1$gc[,sum(npar),by=c(ggcline.object2$test.subject)];
			    npar2=ggcline.object2$gc[,sum(npar),by=c(ggcline.object2$test.subject)];

				lppd1_sum=ggcline.object1$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(ggcline.object2$test.subject)];
				lppd2_sum=ggcline.object2$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(ggcline.object2$test.subject)];

                effpar1_sum=ggcline.object1$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(ggcline.object2$test.subject)];
                effpar2_sum=ggcline.object2$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(ggcline.object2$test.subject)];			
				
				setnames(npar1,"V1","npar1");
				setnames(npar2,"V1","npar2");
				setnames(lppd1_sum,"V1","lppd1_sum");
				setnames(lppd2_sum,"V1","lppd2_sum");							
				setnames(effpar1_sum,"V1","effpar1_sum");
				setnames(effpar2_sum,"V1","effpar2_sum");
				
				setkeyv(npar1,ggcline.object2$test.subject);
				setkeyv(npar2,ggcline.object2$test.subject);
				setkeyv(lppd1_sum,ggcline.object2$test.subject);
				setkeyv(lppd2_sum,ggcline.object2$test.subject);
				setkeyv(effpar1_sum,ggcline.object2$test.subject);
				setkeyv(effpar2_sum,ggcline.object2$test.subject);				
				
				MET1=npar1[lppd1_sum][effpar1_sum];
				
				MET1[,waic1_npar1_AICscale:=-2*lppd1_sum + 2*npar1][,
                waic1_effpar1_AICscale:=-2*lppd1_sum + 2*effpar1_sum];
				
###

				MET2=npar2[lppd2_sum][effpar2_sum];
				
				MET2[,waic2_npar2_AICscale:=-2*lppd2_sum + 2*npar2][,
                waic2_effpar2_AICscale:=-2*lppd2_sum + 2*effpar2_sum];

				setkeyv(MET1,ggcline.object2$test.subject);				
				setkeyv(MET2,ggcline.object2$test.subject);
                MET=MET1[MET2];				
                MET[,
			        waicdiff_npar_AICscale:=waic1_npar1_AICscale - waic2_npar2_AICscale][,
				    waicdiff_effpar_AICscale:=waic1_effpar1_AICscale - waic2_effpar2_AICscale];
				
				#setkeyv(MET,ggcline.object2$test.subject);
				#MET=MET[se_diff_mean];
				#MET[,waicdiff_effpar_mean_quantzero:=pnorm(0,mean=abs(waicdiff_effpar_mean),sd=se_diff_mean,lower.tail=TRUE)];#This is the quantile of zero difference, equivalent to a one-tailed p value#
				#Should result in negative value if more complex model is favoured#

warning("While the effective number of parameters ('effpar') is calculated and presented, and used in some of the waic calculations, testing suggests it is systematically overestimated and should not be trusted. Unless the estimated effpar values are the same as or less than the npar values (the actual number of parameters), use the calculations based on npar");				
				return(MET);

			                                                                             }else{
																						 
#setkey(ggcline.object1$likmeans,index);setkey(ggcline.object2$likmeans,index);
#MET.waic=ggcline.object1$likmeans[ggcline.object2$likmeans];
#se_diff_mean=MET.waic[,sqrt(var(log(postmean.exp.loglik) - log(i.postmean.exp.loglik))),by=c(ggcline.object1$test.subject)];
#setkeyv(se_diff_mean,ggcline.object1$test.subject);
#setnames(se_diff_mean,"V1","se_diff_mean");
			
			    npar1=ggcline.object1$gc[,sum(npar),by=c(ggcline.object1$test.subject)];
			    npar2=ggcline.object2$gc[,sum(npar),by=c(ggcline.object1$test.subject)];

				lppd1_sum=ggcline.object1$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(ggcline.object1$test.subject)];
				lppd2_sum=ggcline.object2$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(ggcline.object1$test.subject)];

                effpar1_sum=ggcline.object1$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(ggcline.object1$test.subject)];
                effpar2_sum=ggcline.object2$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(ggcline.object1$test.subject)];			
				
				setnames(npar1,"V1","npar1");
				setnames(npar2,"V1","npar2");
				setnames(lppd1_sum,"V1","lppd1_sum");
				setnames(lppd2_sum,"V1","lppd2_sum");							
				setnames(effpar1_sum,"V1","effpar1_sum");
				setnames(effpar2_sum,"V1","effpar2_sum");
				
				setkeyv(npar1,ggcline.object1$test.subject);
				setkeyv(npar2,ggcline.object1$test.subject);
				setkeyv(lppd1_sum,ggcline.object1$test.subject);
				setkeyv(lppd2_sum,ggcline.object1$test.subject);
				setkeyv(effpar1_sum,ggcline.object1$test.subject);
				setkeyv(effpar2_sum,ggcline.object1$test.subject);				
				
				MET1=npar1[lppd1_sum][effpar1_sum];
				
				MET1[,waic1_npar1_AICscale:=-2*lppd1_sum + 2*npar1][,
                waic1_effpar1_AICscale:=-2*lppd1_sum + 2*effpar1_sum];
				
###

				MET2=npar2[lppd2_sum][effpar2_sum];
				
				MET2[,waic2_npar2_AICscale:=-2*lppd2_sum + 2*npar2][,
                waic2_effpar2_AICscale:=-2*lppd2_sum + 2*effpar2_sum];

				setkeyv(MET1,ggcline.object1$test.subject);				
				setkeyv(MET2,ggcline.object1$test.subject);
                MET=MET1[MET2];				
                MET[,
			        waicdiff_npar_AICscale:=waic2_npar2_AICscale - waic1_npar1_AICscale][,
				    waicdiff_effpar_AICscale:=waic2_effpar2_AICscale - waic1_effpar1_AICscale];
				
				#setkeyv(MET,ggcline.object1$test.subject);
				#MET=MET[se_diff_mean];
				#MET[,waicdiff_effpar_mean_quantzero:=pnorm(0,mean=abs(waicdiff_effpar_mean),sd=se_diff_mean,lower.tail=TRUE)];					

warning("While the effective number of parameters ('effpar') is calculated and presented, and used in some of the waic calculations, testing suggests it is systematically overestimated and should not be trusted. Unless the estimated effpar values are the same as or less than the npar values (the actual number of parameters), use the calculations based on npar");				
				return(MET);																						 
																						 
																						      };						  
						 
                               };
                                   }else{
        if(sum(names(esth.object1$hi)=="npar_compare")==0){
            MET1=esth.object1$hi[,c(esth.object1$test.subject),with=F];#Lower resolution analysis#1 row per lower resolution test.subject + npar#
            MET1[,npar1:=1];
            MET2=unique(esth.object2$hi[,
                c(esth.object1$test.subject,"npar_compare"),with=F]);   #Higher resolution analysis#1 row per lower resolution test.subject + npar#
            MET2[,npar2:=npar_compare][,npar_compare:=NULL];
			
			setkeyv(MET1,esth.object1$test.subject);
			setkeyv(MET2,esth.object1$test.subject);
			
			
			
#########################*********
#setkey(esth.object1$likmeans,index);setkey(esth.object2$likmeans,index);
#MET.waic=esth.object1$likmeans[esth.object2$likmeans];
#se_diff_mean=MET.waic[,sqrt(var(log(postmean.exp.loglik) - log(i.postmean.exp.loglik))),by=c(esth.object1$test.subject)];
#setkeyv(se_diff_mean,esth.object1$test.subject);
#setnames(se_diff_mean,"V1","se_diff_mean");

				lppd1_sum=esth.object1$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(esth.object1$test.subject)];
				lppd2_sum=esth.object2$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(esth.object1$test.subject)];

                effpar1_sum=esth.object1$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(esth.object1$test.subject)];
                effpar2_sum=esth.object2$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(esth.object1$test.subject)];			

				setnames(lppd1_sum,"V1","lppd1_sum");
				setnames(lppd2_sum,"V1","lppd2_sum");							
				setnames(effpar1_sum,"V1","effpar1_sum");
				setnames(effpar2_sum,"V1","effpar2_sum");
				
				setkeyv(lppd1_sum,esth.object1$test.subject);
				setkeyv(lppd2_sum,esth.object1$test.subject);
				setkeyv(effpar1_sum,esth.object1$test.subject);
				setkeyv(effpar2_sum,esth.object1$test.subject);								
				
				MET1=MET1[lppd1_sum][effpar1_sum];
				
				MET1[,waic1_npar1_AICscale:=-2*lppd1_sum + 2*npar1][,
                waic1_effpar1_AICscale:=-2*lppd1_sum + 2*effpar1_sum];
				
###

				MET2=MET2[lppd2_sum][effpar2_sum];
				
				MET2[,waic2_npar2_AICscale:=-2*lppd2_sum + 2*npar2][,
                waic2_effpar2_AICscale:=-2*lppd2_sum + 2*effpar2_sum];

				setkeyv(MET1,esth.object1$test.subject);				
				setkeyv(MET2,esth.object1$test.subject);
                MET=MET1[MET2];				
                MET[,
			        waicdiff_npar_AICscale:=waic2_npar2_AICscale - waic1_npar1_AICscale][,
				    waicdiff_effpar_AICscale:=waic2_effpar2_AICscale - waic1_effpar1_AICscale];
				
				#setkeyv(MET,esth.object1$test.subject);
				#MET=MET[se_diff_mean];
				#MET[,waicdiff_effpar_mean_quantzero:=pnorm(0,mean=abs(waicdiff_effpar_mean),sd=se_diff_mean,lower.tail=TRUE)];

warning("While the effective number of parameters ('effpar') is calculated and presented, and used in some of the waic calculations, testing suggests it is systematically overestimated and should not be trusted. Unless the estimated effpar values are the same as or less than the npar values (the actual number of parameters), use the calculations based on npar");				
				return(MET);														
#########################*********						

                                                          }else{
														  
            MET1=esth.object2$hi[,c(esth.object2$test.subject),with=F];#Lower resolution analysis#1 row per lower resolution test.subject + npar#
            MET1[,npar1:=1];
            MET2=unique(esth.object1$hi[,
                c(esth.object2$test.subject,"npar_compare"),with=F]);   #Higher resolution analysis#1 row per lower resolution test.subject + npar#
            MET2[,npar2:=npar_compare][,npar_compare:=NULL];
			
			setkeyv(MET1,esth.object2$test.subject);
			setkeyv(MET2,esth.object2$test.subject);						
			
#########################*********
#setkey(esth.object1$likmeans,index);setkey(esth.object2$likmeans,index);
#MET.waic=esth.object1$likmeans[esth.object2$likmeans];
#se_diff_mean=MET.waic[,sqrt(var(log(postmean.exp.loglik) - log(i.postmean.exp.loglik))),by=c(esth.object2$test.subject)];
#setkeyv(se_diff_mean,esth.object2$test.subject);
#setnames(se_diff_mean,"V1","se_diff_mean");

				lppd1_sum=esth.object2$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(esth.object2$test.subject)];
				lppd2_sum=esth.object1$likmeans[,sum(log(postmean.exp.loglik),na.rm=T),by=c(esth.object2$test.subject)];

                effpar1_sum=esth.object2$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(esth.object2$test.subject)];
                effpar2_sum=esth.object1$likmeans[,sum(postmean.loglik.sqrd - postmean.loglik^2, na.rm=T),by=c(esth.object2$test.subject)];			

				setnames(lppd1_sum,"V1","lppd1_sum");
				setnames(lppd2_sum,"V1","lppd2_sum");							
				setnames(effpar1_sum,"V1","effpar1_sum");
				setnames(effpar2_sum,"V1","effpar2_sum");
				
				setkeyv(lppd1_sum,esth.object2$test.subject);
				setkeyv(lppd2_sum,esth.object2$test.subject);
				setkeyv(effpar1_sum,esth.object2$test.subject);
				setkeyv(effpar2_sum,esth.object2$test.subject);								
				
				MET1=MET1[lppd1_sum][effpar1_sum];
				
				MET1[,waic1_npar1_AICscale:=-2*lppd1_sum + 2*npar1][,
                waic1_effpar1_AICscale:=-2*lppd1_sum + 2*effpar1_sum];

				MET2=MET2[lppd2_sum][effpar2_sum];
				
				MET2[,waic2_npar2_AICscale:=-2*lppd2_sum + 2*npar2][,
                waic2_effpar2_AICscale:=-2*lppd2_sum + 2*effpar2_sum];

				setkeyv(MET1,esth.object2$test.subject);				
				setkeyv(MET2,esth.object2$test.subject);
                MET=MET1[MET2];				
                MET[,
			        waicdiff_npar_AICscale:=waic2_npar2_AICscale - waic1_npar1_AICscale][,
				    waicdiff_effpar_AICscale:=waic2_effpar2_AICscale - waic1_effpar1_AICscale];
				
				#setkeyv(MET,esth.object2$test.subject);
				#MET=MET[se_diff_mean];
				#MET[,waicdiff_effpar_mean_quantzero:=pnorm(0,mean=abs(waicdiff_effpar_mean),sd=se_diff_mean,lower.tail=TRUE)];

warning("While the effective number of parameters ('effpar') is calculated and presented, and used in some of the waic calculations, testing suggests it is systematically overestimated and should not be trusted. Unless the estimated effpar values are the same as or less than the npar values (the actual number of parameters), use the calculations based on npar");				
				return(MET);																  														  													  

                                                               };
                                        };

                                        }
