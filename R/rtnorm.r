#' Internal function that samples from a truncated normal distribution.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @export															 
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