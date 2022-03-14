#' Function to calculate AIC from the best posterior parameter estimate(s) for an esth or ggcline object.
#' 
#' @param data.prep.object Name of the \code{data.prep} object produced by the \code{data.prep} function.
#' @param esth.object Name of the \code{esth} object.
#' @param esth.colname Name of the column in the \code{esth.object} with the best posterior hybrid index estimates. 
#'   Default \sQuote{h_posterior_mode}.
#' @param ggcline.object Name of the \code{ggcline} object. Only needed when calculating AIC for a genomic cline model.
#' @param ggcline.pooled Logical. If set to \option{TRUE}, calculates a single AIC for the whole dataset. Useful when 
#'   a comparison is being made between a model with one of the genomic cline parameters pooled across the whole dataset 
#'   and another model. Default \option{FALSE}.
#' @param test.subject The test subject for which you wish to make a comparison between models.
#' @details \code{calc_AIC} calculates the AIC for one or more test subjects for a single model, using the likelihoods for 
#'   the best posterior parameter estimates.
#' 
#' For hybrid index model comparison, \code{calc_AIC} is set up to compare models of differing resolution, for 
#'   example, for each of a set of populations, comparing AIC for a single hybrid index estimated at the population level 
#'   to separate hybrid index estimates per individual within the population.
#' 
#' For genomic cline model comparisons there are further options. Comparisons can be made between models run at the same resolution 
#'   but with differing combinations of parameters fixed versus estimated, or parameters fixed at different values. When one or both 
#'   cline parameters are pooled across the dataset in one model, it can be compared at the level of the whole dataset with models 
#'   with no or different pooling.
#' @return  A \code{data.table} and \code{data.frame} with the name of the test subject (when there is more than one) plus following fields:
#'   \item{ml}{The model maximum likelihood for each test subject.}
#'   \item{npar}{The number of parameters for each test subject.}
#'   \item{AIC}{AIC for each test subject.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' ############################################################################################################
#' #Hybrid index example: compare the fit between running esth at the individual level versus population level#
#' ############################################################################################################
#' 
#' #This will result in a separate model comparison for each POPID (the level with lower resolution).
#' 
#' #Use the sparrow dataset, as in other examples.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);
#' 
#' #First run esth with test.subject="INDLABEL" (default). 
#' #'return.likmeans=TRUE' is not needed to calculate AIC; it is only necessary to calculate waic using 'compare.models'.
#' #For the higher resolution test subject only, we also need to declare what lower resolution test subject we're going to compare with.
#' hindlabel=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  test.subject = "INDLABEL",                                  #default#
#'  test.subject.compare = "POPID",                             #***the test.subject you want to compare with the current test.subject***#
#'  include.Source=TRUE,
#'  nitt=3000,
#'  burnin=1000
#' )
#' 
#' #Then run esth again, with test.subject="POPID". No need to declare 'test.subject.compare' when running the lower resolution model.
#' hindlabel2=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  test.subject = "POPID",                                     #***Changed to "POPID"***#
#'  include.Source=TRUE,
#'  nitt=3000,
#'  burnin=1000
#' )
#' 
#' #Now compare models using AIC.
#' #AIC is calculated per declared test subject, so for direct comparison use the lower-resolution test subject to calculate AIC for both models.
#' aic_INDLABEL=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel,
#'  test.subject="POPID"              #Make this the same for both models for direct comparison of fit#
#' );
#' 
#' aic_POPID=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel2,
#'  test.subject="POPID"              #Make this the same for both models for direct comparison of fit#
#' )
#' 
#' #Now join the tables and calculate the difference in AIC.
#' setkey(aic_INDLABEL,POPID);setkey(aic_POPID,POPID);
#' aic_both=aic_INDLABEL[aic_POPID];#The columns named 'i...' refer to the second table, aic_POPID#
#' aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (a single pooled hybrid index for each populations)#
#' 
#' #Show the results in order of aicdiff for each POPID. Those at the top of the ordered table more strongly favour the individual-level model.
#' aic_both[order(aicdiff)]
#' 
#' #################################
#' #################################
#' #Genomic cline model comparisons#
#' #################################
#' #################################
#' 
#' ##############################################################
#' #Parameter v pooled versus unpooled across a whole chromosome#
#' ##############################################################
#' 
#' #Calculate a single AIC for the sparrow Z chromosome for two ggcline models run in the examples under ggcline, 
#' #one with v pooled and one unpooled.
#' aic_unpooled_z=calc_AIC(
#'  data.prep.object=prepZ,
#'  esth.object=hindlabel,
#'  ggcline.object=gc_unpooled_z,
#'  test.subject="locus",
#'  ggcline.pooled=TRUE
#' )
#' 
#' aic_poolv_Z=calc_AIC(
#'  data.prep.object=prepZ,
#'  esth.object=hindlabel,
#'  ggcline.object=gc_poolv_Z,
#'  test.subject="locus",
#'  ggcline.pooled=TRUE
#' )
#' 
#' #Calculate the difference in AIC. Negative favours the more complex model, with v unpooled.
#' aic_unpooled_z$AIC - aic_poolv_Z$AIC
#' 
#' #############################################################################################################
#' #At the individual locus level, compare a model with v and centre estimated versus a model with v fixed at 1#
#' #############################################################################################################
#' 
#' #Compare the ggcline objects gc (no fixing) and gc_v_fixed created in the ggcline examples.
#' aic_v_fixed=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel,
#'  ggcline.object=gc_v_fixed,
#'  test.subject="locus",
#'  ggcline.pooled=FALSE
#' )
#' 
#' aic_no_fixed=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel,
#'  ggcline.object=gc,
#'  test.subject="locus",
#'  ggcline.pooled=FALSE
#' )
#' 
#' #Now join the tables and calculate the difference in AIC for each locus.
#' setkey(aic_no_fixed,locus);setkey(aic_v_fixed,locus);
#' aic_both=aic_no_fixed[aic_v_fixed];#The columns named 'i...' refer to the second table, aic_v_fixed#
#' aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (v fixed at 1)#
#' 
#' #Show the results in order of aicdiff for each POPID. Those at the top of the ordered table more strongly favour the individual-level model.
#' aic_both[order(aicdiff)]
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
#' #Now run AIC for each gc object.
#' aic_gc=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel,
#'  ggcline.object=gc,
#'  test.subject="chrom",           #Use the lower resolution test subject#
#'  ggcline.pooled=FALSE
#' )
#' 
#' aic_gc_by_chr=calc_AIC(
#'  data.prep.object=prepdata$data.prep,
#'  esth.object=hindlabel,
#'  ggcline.object=gc_by_chr,
#'  test.subject="chrom",           #Use the lower resolution test subject#
#'  ggcline.pooled=FALSE
#' )
#' 
#' #Now join the tables and calculate the difference in AIC for each chromosome.
#' setkey(aic_gc,chrom);setkey(aic_gc_by_chr,chrom);
#' aic_both=aic_gc[aic_gc_by_chr];#The columns named 'i...' refer to the second table, aic_gc_by_chr#
#' aic_both[,aicdiff:=AIC - i.AIC]#Positive values favour the less complex model (single cline for the whole chromosome)#
#' 
#' #Show the results in order of aicdiff for each POPID. Those at the top of the ordered table more strongly favour the individual locus level model.
#' aic_both[order(aicdiff)]
#' }
#' @export
calc_AIC=function(data.prep.object,
    esth.object=NULL,
    esth.colname="h_posterior_mode",
    ggcline.object=NULL,             #LEAVE BLANK IF YOU WANT AIC FOR ESTH RESULTS#
    ggcline.pooled=FALSE,
    test.subject
         ){

    #if(length(esth.object) > 1){
     #   stop("AIC can only be calculated for one results object at a time")
      #                                                      };

    #if(length(ggcline.object) > 1){
     #   stop("AIC can only be calculated for one results object at a time")
      #                                                      };
###############
    if(is.null(ggcline.object)==F){

        if(ggcline.pooled){
            npar=unique(ggcline.object$gc[,sum_npar]);

###MAX LIKELIHOOD CALCULATION###
setkeyv(ggcline.object$gc,test.subject);
setkeyv(data.prep.object,test.subject);
obj=ggcline.object$gc[,c(test.subject,"exp_mean_log_v","mean_logit_centre"),with=FALSE][data.prep.object];

setkeyv(obj,esth.object$test.subject);
setkeyv(esth.object$hi,esth.object$test.subject);
obj=esth.object$hi[,c(esth.object$test.subject,esth.colname),with=FALSE][obj];

setnames(obj,c("exp_mean_log_v","mean_logit_centre",esth.colname),c("v","centre","h"));

obj[,u:=centre*v];

setkey(obj,Source_allele);
obj[.(1),
    loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_1 + 
        (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_1)][.(0),
    loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*(1 - S1.prop_1) + 
        (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*(1 - S0.prop_1))];

sum.loglik=data.table(obj[,sum(loglik.gc, na.rm=T)]);
setnames(sum.loglik,"V1","ml");
################################

        sum.loglik[,
            npar:=npar][,
            AIC:=-2*ml + 2*npar];			

        return(sum.loglik);
                          }else{

        npar=ggcline.object$gc[,sum(npar),by=c(test.subject)];
        setnames(npar,"V1","npar");
        setkeyv(npar,test.subject);

###MAX LIKELIHOOD CALCULATION###
setkeyv(ggcline.object$gc,ggcline.object$test.subject);#test.subject);
setkeyv(data.prep.object,ggcline.object$test.subject);#test.subject);
obj=ggcline.object$gc[,c(test.subject,ggcline.object$test.subject,"exp_mean_log_v","mean_logit_centre"),with=FALSE][data.prep.object];#***PROBLEM LINE***#

setkeyv(obj,esth.object$test.subject);
setkeyv(esth.object$hi,esth.object$test.subject);
obj=esth.object$hi[,c(esth.object$test.subject,esth.colname),with=FALSE][obj];

setnames(obj,c("exp_mean_log_v","mean_logit_centre",esth.colname),c("v","centre","h"));

obj[,u:=centre*v];

setkey(obj,Source_allele);
obj[.(1),
    loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_1 + 
        (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_1)][.(0),
    loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*(1 - S1.prop_1) + 
        (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*(1 - S0.prop_1))];

sum.loglik=obj[,sum(loglik.gc, na.rm=T),by=test.subject];
setnames(sum.loglik,"V1","ml");
################################

        setkeyv(sum.loglik,test.subject);
        sum.loglik=sum.loglik[npar];

        sum.loglik[,AIC:=-2*ml + 2*npar];			

        return(sum.loglik);

                               };
}else{#ESTH ONLY***

###MAX LIKELIHOOD CALCULATION###
setkeyv(data.prep.object,esth.object$test.subject);
setkeyv(esth.object$hi,esth.object$test.subject);
obj=esth.object$hi[,c(esth.object$test.subject,test.subject,esth.colname),with=FALSE][data.prep.object];

setnames(obj,c(esth.colname),c("h"));
setkey(obj,Source_allele);
obj[.(1),loglik:=log(h*S1.prop_1 + (1-h)*S0.prop_1)];
obj[.(0),loglik:=log(h*(1 - S1.prop_1) + (1-h)*(1 - S0.prop_1))];

sum.loglik=obj[,sum(loglik, na.rm=T),by=test.subject];
setnames(sum.loglik,"V1","ml");
###MAX LIKELIHOOD CALCULATION###

    npar=esth.object$hi[,.N,by=test.subject];
    setnames(npar,"N","npar");
    setkeyv(npar,test.subject);
    setkeyv(sum.loglik,test.subject);
    sum.loglik=sum.loglik[npar];

    sum.loglik[,AIC:=-2*ml + 2*npar];

    return(sum.loglik);

};


}#END OF FUNCTION#