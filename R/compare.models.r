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
#' @details \code{compare.models} uses the actual number of parameters for each model, rather than an estimate,
#'   which would be applicable to more complex hierarchichal models.
#' @return A \code{data.table} and \code{data.frame} with the name of the test subject plus following fields:*************
#'   \item{npar1}{number of parameters for the first model.}
#'   \item{lppd1}{Sum of the mean posterior likelihoods for the individual data points for the first model.}
#'   \item{waic1}{widely applicable information criterion for the first model.}
#'   \item{se1}{standard error of waic1.}
#'   \item{npar2}{number of parameters for the second model.}
#'   \item{lppd2}{Sum of the mean posterior likelihoods for the individual data points for the second model.}
#'   \item{waic2}{widely applicable information criterion for the second model.}
#'   \item{se2}{standard error of waic2.}
#'   \item{diff}{difference in waic between the first and second models. A positive value favours the first model.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' #A comparison of individual-level hybrid index estimates with population-level hybrid index estimates, to see if there is 
#' #evidence of individuals within populations differing in hybrid index:
#' esthcomp=compare.models(hest.object1=hindlabel,
#' hest.object2=hpopid)
#' #A comparison of genomic cline estimates versus likelihoods of genomic clines with the two parameters fixed to null values (an 
#' #alternative way of testing significance to the p values presented in \code{ggcline}):
#' gclinecomp=compare.models(ggcline.object1=gc1,
#' ggcline.object2=gcAllFixed)
#' }
#' @export
compare.models=function(ggcline.object1=NULL,ggcline.object2=NULL,
    ggcline.pooled=FALSE,
    esth.object1=NULL,esth.object2=NULL){

    if(is.null(ggcline.object1)==F){
        if(ggcline.pooled){
            npar1=unique(ggcline.object1$gc[,sum_npar]);
            MET.waic1=ggcline.object1$likmeans;

            MET.1=data.table(lppd1=MET.waic1[,sum(log(postmean.exp.loglik), na.rm=T)]);
            MET.1[,
                npar1:=npar1][,
                waic1:=-2*lppd1 + 2*npar1][,
                se1:=sqrt(MET.waic1[,.N*var(log(postmean.exp.loglik), na.rm=T)])];

            npar2=unique(ggcline.object2$gc[,sum_npar]);
            MET.waic2=ggcline.object2$likmeans;

            MET.2=data.table(lppd2=MET.waic2[,sum(log(postmean.exp.loglik), na.rm=T)]);
            MET.2[,
                npar2:=npar2][,
                waic2:=-2*lppd2 + 2*npar2][,
                se2:=sqrt(MET.waic2[,.N*var(log(postmean.exp.loglik), na.rm=T)])];

            MET.both=cbind(MET.1,MET.2);
            MET.both[,diff:=waic2 - waic1];

            print(MET.both);

            return(MET.both);
                          }else{ 
            setkeyv(ggcline.object1$gc,ggcline.object1$test.subject);
            setkeyv(ggcline.object2$gc,ggcline.object2$test.subject);
            setkeyv(ggcline.object1$likmeans,ggcline.object1$test.subject);
            setkeyv(ggcline.object2$likmeans,ggcline.object2$test.subject);
            MET.waic1=ggcline.object1$likmeans[,
                c(ggcline.object1$test.subject,"postmean.exp.loglik"),with=F];
            MET.waic2=ggcline.object2$likmeans[,
                c(ggcline.object1$test.subject,"postmean.exp.loglik"),with=F];

            setkeyv(MET.waic1,ggcline.object1$test.subject);
            MET.1=data.table(ggcline.object1$gc[,c(ggcline.object1$test.subject),with=F]);
            MET.1[,
                npar1:=ggcline.object1$gc[,npar]][,
                lppd1:=MET.waic1[,sum(log(postmean.exp.loglik), 
                na.rm=T),by=c(ggcline.object1$test.subject)]$V1][,
                waic1:=-2*lppd1 + 2*npar1][,
                se1:=sqrt(MET.waic1[,.N*var(log(postmean.exp.loglik), 
                na.rm=T),by=c(ggcline.object1$test.subject)]$V1)];

            setkeyv(MET.waic2,ggcline.object2$test.subject);
            MET.2=data.table(ggcline.object2$gc[,c(ggcline.object2$test.subject),with=F]);
            MET.2[,
                npar2:=ggcline.object2$gc[,npar]][,
                lppd2:=MET.waic2[,sum(log(postmean.exp.loglik), 
			        na.rm=T),by=c(ggcline.object2$test.subject)]$V1][,
                waic2:=-2*lppd2 + 2*npar2][,
                se2:=sqrt(MET.waic2[,.N*var(log(postmean.exp.loglik), 
                    na.rm=T),by=c(ggcline.object2$test.subject)]$V1)];

            setkeyv(MET.1,ggcline.object1$test.subject);
            setkeyv(MET.2,ggcline.object2$test.subject);
            MET.both=MET.1[MET.2];
            MET.both[,diff:=waic2 - waic1];

            print(MET.both);

            return(MET.both);
                               };
                                   }else{
        if(sum(names(esth.object1$hi)=="npar_compare")==0){
            MET.1=esth.object1$hi[,c(esth.object1$test.subject),with=F];
            MET.1[,npar1:=1];
            MET.2=unique(esth.object2$hi[,
                c(esth.object1$test.subject,"npar_compare"),with=F]);
            MET.2[,npar2:=npar_compare][,npar_compare:=NULL];

            MET.1[,
			    lppd1:=esth.object1$likmeans[,sum(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object1$test.subject)]$V1][,
                waic1:=-2*lppd1 + 2*npar1][,
                se1:=sqrt(esth.object1$likmeans[,.N*var(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object1$test.subject)]$V1)];

            MET.2[,
			    lppd2:=esth.object2$likmeans[,sum(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object1$test.subject)]$V1][,
                waic2:=-2*lppd2 + 2*npar2][,
                se2:=sqrt(esth.object2$likmeans[,.N*var(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object1$test.subject)]$V1)];

            setkeyv(MET.1,esth.object1$test.subject);
            setkeyv(MET.2,esth.object1$test.subject);
            MET.both=MET.1[MET.2];
            MET.both[,diff:=waic2 - waic1];

            print(MET.both);

            return(MET.both);

                                                          }else{
            MET.1=esth.object2$hi[,c(esth.object2$test.subject),with=F];
            MET.1[,npar1:=1];
            MET.2=unique(esth.object1$hi[,
                c(esth.object2$test.subject,"npar_compare"),with=F]);
            MET.2[,npar2:=npar_compare][,npar_compare:=NULL];

            MET.1[,
			    lppd1:=esth.object2$likmeans[,sum(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object2$test.subject)]$V1][,
                waic1:=-2*lppd1 + 2*npar1][,
                se1:=sqrt(esth.object2$likmeans[,.N*var(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object2$test.subject)]$V1)];

            MET.2[,
			    lppd2:=esth.object1$likmeans[,sum(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object2$test.subject)]$V1][,
                waic2:=-2*lppd2 + 2*npar2][,
                se2:=sqrt(esth.object1$likmeans[,.N*var(log(postmean.exp.loglik), 
                    na.rm=T),by=c(esth.object2$test.subject)]$V1)];

            setkeyv(MET.1,esth.object2$test.subject);
            setkeyv(MET.2,esth.object2$test.subject);
            MET.both=MET.2[MET.1];
            MET.both[,diff:=waic2 - waic1];

            print(MET.both);

            return(MET.both);
                                                               };
                                        };
                                        }