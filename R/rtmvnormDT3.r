#' Function to sample from the posterior of estimated genomic cline parameters.
#' 
#' @param nsamp Numeric. Number of random posterior samples. Default \kbd{1}.
#' @param meanlogv Name of the column in a \code{data.table} containing posterior mean values for ln(\code{v}), or another variable from a bivariate normal distribution. 
#'   Entry should not include quotation marks. No default.
#' @param meanlogitcentre Name of the column in a \code{data.table} containing mean values for logit(\code{centre}), or another variable from a 
#'   bivariate normal distribution. Entry should not include quotation marks. No default.
#' @param varlogv Name of the column in a \code{data.table} containing values for the posterior variance of ln(\code{v}), or another variable from a 
#'   bivariate normal distribution. Entry should not include quotation marks. No default.
#' @param varlogitcentre Name of the column in a \code{data.table} containing values for the posterior variance of logit(\code{centre}), or another variable from a 
#'   bivariate normal distribution. Entry should not include quotation marks. No default.
#' @param covlogvlogitcentre Name of the column in a \code{data.table} containing the posterior covariance between ln(\code{v}) and logit(\code{centre}), or 
#'   any other pair of bivariate normally distributed variables. Entry should not include quotation marks. No default.
#' @param lower Numeric vector. If truncation is desired, the lower cutoff for the two variables.
#' @param upper Numeric vector. If truncation is desired, the upper cutoff for the two variables.
#' @details \code{rtmvnormDT3} is a generalized version of two internal \code{gghybrid} functions that themselves are designed for sampling from the parameter 
#'   proposal distributions for ln(\code{v}) and logit(\code{centre}). Its purpose is to take samples from the joint posterior of the two parameters for one 
#'   or more test subjects simultaneously. While the option names are specific to \code{gghybrid}, the column names are left open and so the function can be 
#'   used to sample from any set of bivariate normal distributions, including with truncation.
#' 
#' The function must be run from within a \code{data.table}, which can be subsetted by row if desired, as shown in the example.
#' @return A \code{data.table} and \code{data.frame} with one row per sample per test subject.
#' 
#'   \item{locus}{the name of the chosen test.subject. If a test subject other than locus is chosen, the column will be named accordingly.}
#'   \item{log_v}{posterior samples of ln(v), or another chosen variable.}
#'   \item{logit_centre}{posterior samples of logit(centre), or another chosen variable.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' ##########################################################################
#' #Sample from the genomic cline parameter posteriors of a subset of 5 loci#
#' ##########################################################################
#' 
#' #The function must be run from within a data.table, after the comma. Subsetting by row can be carried out.
#' #In this case I'm sampling from the 'gc' object created in the documentation examples for 'ggcline'.
#' 
#' gcsamp=gc$gc[locus%in%c("ACO1","CLTA","LNPEP","HSDL2","MCCC2"), #subsetting the data.table by row to only include these 5 loci#
#'  rtmvnormDT3(
#'   nsamp=500,                                 #default is 1 sample per locus#
#'   meanlogv=mean_log_v,                       #no default, column must be specified#
#'   meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
#'   varlogv=var_log_v,                         #no default, column must be specified#
#'   varlogitcentre=var_logit_centre,           #no default, column must be specified#
#'   covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
#'   lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
#'   upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
#'             ),                               #end of rtmvnormDT3 function#
#'  by="locus"];                                #indicate a grouping variable for the sampling. This will be the name of the first column in the resulting table#
#' }
#' @export
rtmvnormDT3=function(
 nsamp=1,#Number of samples#
 meanlogv,
 meanlogitcentre,
 varlogv,
 varlogitcentre,
 covlogvlogitcentre,
 lower=c(-Inf,-Inf), upper=c(Inf,Inf)){

for(i in 1:nsamp){

    retval2=numeric();
    sigma=corpcor::make.positive.definite(matrix(c(varlogv,covlogvlogitcentre,covlogvlogitcentre,varlogitcentre),2,2));
    R=chol(sigma, pivot=TRUE);
    pivot=attr(R, "pivot");
    R=R[, order(pivot)];
    while(length(retval2) < 2){
        retval=as.numeric((rnorm(2) %*% R) + c(meanlogv,meanlogitcentre));
        logv=retval[1];logitcentre=retval[2];
        logv=logv[logv >=lower[1] & logv <=upper[1]];
        logitcentre=logitcentre[logitcentre >=lower[2] & logitcentre <=upper[2]];
        retval2 <- c(logv,logitcentre);
                              };
    stopifnot(length(retval2)==2);
    out<- data.table(retval2[1], retval2[2]);

    if(i==1){outM=out}else{outM=rbind(outM,out)};

};#endfor i#

    setnames(outM,c("V1","V2"),c("log_v","logit_centre"));

    return(outM)
                                     };#End of function####