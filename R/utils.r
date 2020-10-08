#' Internal functions not to be exported.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}

rtnorm=function (n, mean=0, sd=1, lower=-Inf, upper=Inf){
    ret=numeric();
    if (length(n) > 1)
        n=length(n);
    while (length(ret) < n) {
        y=rnorm(n - length(ret), mean, sd);
        y=y[y >=lower & y <=upper];
        ret=c(ret, y)
                            };
    stopifnot(length(ret)==n);
    ret
                                                        };
##ggcline##
rtlnorm=function (n, mean=0, sd=1, lower=-Inf, upper=Inf){
    ret=numeric();
    if (length(n) > 1)
        n=length(n);
    while (length(ret) < n) {
        y=rlnorm(n - length(ret), mean, sd);
        y=y[y >=lower & y <=upper];
        ret=c(ret, y)
                            };
    stopifnot(length(ret)==n);
    ret
                                                         };

rtmvnormDT=function(mean1,mean2,M11,M12,M22,
    lower=c(-7,-Inf), upper=c(7,Inf)){
    retval2=numeric();
    sigma=corpcor::make.positive.definite(matrix(c(M11,M12,M12,M22),2,2));
    R=chol(sigma, pivot=TRUE);
    pivot=attr(R, "pivot");
    R=R[, order(pivot)];
    while(length(retval2) < 2){
        retval=as.numeric((rnorm(2) %*% R) + c(mean1,mean2));
        logv=retval[1];centre=retval[2];
        logv=logv[logv >=lower[1] & logv <=upper[1]];
        centre=centre[centre >=lower[2] & centre <=upper[2]];
        retval2 <- c(logv,centre);
                              };
    stopifnot(length(retval2)==2);
    out<- list(exp(retval2[1]), retval2[2]);
    return(out)
                                     };

rtmvnormDT2=function(x,lenv,lenu,
    lower=c(-7,-Inf), upper=c(7,Inf)){
    retval2=numeric();
    mean1=x[1:lenv,log(v)];
    mean2=x[1:lenu,centre];
    sigma=diag(c(x[1:lenv,indM11],
        x[1:lenu,indM22]));
    if(lenv<lenu){
        sigma[row(sigma)>col(sigma) & col(sigma)<=lenv|
        row(sigma)<col(sigma) & row(sigma)<=lenv]=
        x[1:max(lenv,lenu),indM12];
                 }else{
        sigma[row(sigma)>col(sigma) & row(sigma)==lenv+1|
        row(sigma)<col(sigma) & col(sigma)==lenv+1]=
        x[1:max(lenv,lenu),indM12];
                      };
    sigma=corpcor::make.positive.definite(sigma);
    R=chol(sigma, pivot=TRUE);
    pivot=attr(R, "pivot");
    R=R[, order(pivot)];
    while(length(retval2) < lenv + lenu){
        retval=as.numeric((rnorm(lenv + lenu) %*% R) + c(mean1,mean2));
        logv=retval[1:lenv];
        centre=retval[lenv+1:lenu];
        logv=logv[logv >=lower[1] & logv <=upper[1]];
        centre=centre[centre >=lower[2] & centre <=upper[2]];
        retval2 <- c(logv,centre);
                                        };
    stopifnot(length(retval2)==lenv + lenu);
    out<- list(exp(retval2[1:lenv]), 
        retval2[lenv+1:lenu]);
    return(out)
                                     };				
									 
##data.prep##
#meanfun=function(x){mean(x,na.rm=T)*ploidy}#??No longer needed??#

replace.fun.A = function(X){replace(X, 
	X==sort(unique(na.omit(X)))[1], 
	"A")}

replace.fun = function(X){replace(X, 
	X !="A", 
	"B")}
	
##read.data##
unique.fun=function(X){sort(unique(na.omit(X)))};
