#' Bayesian hybrid index estimation.
#'
#' @param data.prep.object Name of the \code{data.prep} object produced by the \code{data.prep} function.
#' @param read.data.precols A \code{precols} object produced by \code{read.data}, or a custom character vector with 
#'   names of all pre-marker columns in the \code{data}.
#' @param test.subject Character string. Name of the field identifying subjects for hybrid index estimation. Default \dQuote{INDLABEL}.
#' @param test.subject.compare Character string. Name of a different field identifying subjects for 
#'   hybrid index estimation, to be run through \code{esth} separately and compared with the current run 
#'   downstream using \code{compare.models}.
#' @param include.Source Logical. Whether to estimate hybrid indices for test subjects declared as parental reference. Default is \code{TRUE}.
#' @param return.likmeans Logical. Whether to return a table of likelihood calculations. Required for \code{compare.models}.
#' @param fix.subject Logical default, otherwise a character vector. List of test subjects (normally individual references 
#'   in the \dQuote{INDLABEL} field) for which you wish to fix rather than estimate the hybrid index. Potentially useful for model 
#'   comparison with \code{compare.models}.
#' @param fix.value Logical default, otherwise a numeric vector of length 1 or the same length as \code{fix.subject}. 
#'   The fixed hybrid index value(s) for the test subjects listed in \code{fix.subject}.
#' @param plot.ind Character vector. List of up to 6 test subjects for which you wish to plot the posterior distribution values in real time. 
#'   Default \code{NULL}.
#' @param plot.col Character vector. List of up to 6 colours to plot the \code{plot.ind} test subjects. Default \code{NULL}.
#' @param nitt Numeric. The total number of MCMC iterations including burnin. At least 6000 is recommended.
#' @param burnin Numeric. The number of burnin MCMC iterations. At least 3000 is recommended.
#' @param init.var Numeric. Starting value for the variance of the parameter proposal distribution. For data sets with 
#'   very large numbers of loci it may be useful to reduce this number as uncertainty in hybrid index will be very small. Default 
#'   is 0.002, which works well with tested data sets.
#' @param prior Numeric vector. Values of \code{shape1} and \code{shape2} parameters for the beta-distributed prior. 
#'   Default is \code{shape1=shape2=0.5} (Jeffrey's prior).
#' @param print.k Numeric. The iteration is printed to screen on multiples of this number. Default is \code{50}.
#' @details \code{esth} estimates hybrid index using the likelihood formulae of Buerkle (2005), with the addition of a
#'   prior. The default is Jeffrey's prior (beta distribution with \code{shape1=shape2=0.5}), which testing suggests 
#'   is an improvement over a uniform prior (\code{shape1=shape2=1}).
#'
#' \code{esth} excludes any pre-marker columns from the output with finer resolution than the declared \code{test.subject}, 
#'   and hence has the same number of rows as the number of unique \code{test.subject} values.
#'
#' Set \code{return.likmeans=TRUE} if you intend to carry out model comparison.
#' @return list containing \code{hi}, \code{test.subject} and \code{likmeans} (optional). \code{hi} is a \code{data.table} 
#'   and \code{data.frame} containing all pre-marker columns that do not have finer resolution 
#'   than the declared \code{test.subject}, 
#'   as well as the following fields:
#'   \item{Source}{whether the sample is from source 0 ('S0'), source 1 ('S1') or is a test individual ('TEST').}
#'   \item{h_posterior_mode}{the hybrid index estimate (mode of the beta-distributed posterior).}
#'   \item{h_cred_int_lower}{lower 95 percent credible interval (2.5 percent quantile of the posterior beta distribution).}
#'   \item{h_cred_int_upper}{upper 95 percent credible interval (97.5 percent quantile of the posterior beta distribution).}
#'   \item{beta_mean}{mean of the posterior beta distribution (included for convenience).}
#'   \item{beta_var}{variance of the posterior beta distribution (included for convenience).}
#'   \item{beta_shape1}{beta shape1 parameter estimate for the posterior beta distribution (included for convenience).}
#'   \item{beta_shape2}{beta shape2 parameter estimate for the posterior beta distribution (included for convenience).}
#'   \item{npar_compare}{see compare.models.}
#'   \item{npar}{number of parameters. see compare.models.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @references
#'   Buerkle, C. A. (2005). Maximum likelihood estimation of a hybrid index based on molecular markers. 
#'   Molecular Ecology Notes, 5(3), 684-687.
#' @examples
#'
#' \dontrun{
#' hindlabel= esth(data.prep.object = prepdata$data.prep,
#' read.data.precols = dat$precols,
#' test.subject="INDLABEL",
#' test.subject.compare="POPID",
#' include.Source = TRUE,
#' return.likmeans = TRUE,
#' plot.ind = c("indref1","indref2","indref3","indref4","indref5","indref6"),
#' plot.col = c("blue","green","cyan","purple","magenta","red"),
#' nitt=15000,burnin=5000)
#' }
#' @export
esth= function(data.prep.object,read.data.precols,
    test.subject="INDLABEL",test.subject.compare="INDLABEL",
    include.Source = TRUE,
    return.likmeans = FALSE,fix.subject = FALSE,fix.value = FALSE,
    plot.ind = NULL,plot.col = NULL,nitt,burnin,
	init.var=0.002,
##NEW##
init.var2=0.00005,
#######
prior=c(0.5,0.5),
    print.k=50){

    if(return.likmeans==TRUE){
        data.prep.object[,index:=seq(.N)];
        setkey(data.prep.object,index);
        d3=data.table(index=data.prep.object[,index]);
        setkey(d3,index);
        d3=d3[data.prep.object[,c("index",read.data.precols,test.subject,
            "locus","Source_allele"),with=F]];
        d3[,loglik:=0][,postmean.exp.loglik:=0][,postmean.exp.loglik.1:=0][,
            postmean.loglik:=0][,postmean.loglik.1:=0][,
            postvar.true.loglik:=0][,postvar.loglik:=0][,postvar.loglik.1:=0];
                           };

    if(!include.Source){
        setkey(data.prep.object,Source);
        data.prep.object = data.prep.object[Source == "TEST"];
                   };

 #   read.data.precols.red=character(0);
 #   for(i in 1:length(read.data.precols)){
 #       if(uniqueN(data.prep.object[,c(test.subject,read.data.precols[i]),with=F]) <= 
 #           uniqueN(data.prep.object[,test.subject,with=F])){
 #               read.data.precols.red=c(read.data.precols.red,read.data.precols[i]);
 #                                                     };
 #                                        };
	
    read.data.precols.red=test.subject;#character(0);#New edit 14 Oct 2021#
    for(i in 1:length(read.data.precols)){
        if(uniqueN(data.prep.object[,c(test.subject,read.data.precols[i]),with=F]) <= 
            uniqueN(data.prep.object[,test.subject,with=F])){
                read.data.precols.red=c(read.data.precols.red,read.data.precols[i]);
                                                      };
                                         };
					 
    read.data.precols.red=unique(read.data.precols.red);#New line 14 Oct 2021#

    MET = unique(data.prep.object[,c("Source",read.data.precols.red),with=F])[,
        LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
        TEMPF:=exp(DF)][,h:=runif(.N)][,h1:=h][,
        indM:=init.var][,m:=(2.4/sqrt(1))^2][,
        accept := 0][,r:=0.99][,
        Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
        Astar:=0.44][,q:=2000];

    if(fix.subject!=FALSE){
        setkeyv(MET,test.subject);
        MET[fix.subject,h:=fix.value][fix.subject,h1:=h];
                             };

    if(is.null(test.subject.compare)==F){
        MET[,npar_compare:=.N,by=test.subject.compare];
                                        };

    for(k in 1:nitt){
        if(k==1){cat(paste("Estimating hybrid index..."),fill=1); flush.console();
		        };

        if(tester::is_multiple(k,print.k)){cat(paste("\t","\t","\t","Iteration",k,"; acceptance=",unique(MET[,A]),"; indM",unique(MET[,indM])),fill=1); 
            flush.console();     
			                      };

        setkey(MET,h);
        MET[,prior.h:=dbeta(h,shape1=prior[1],shape2=prior[2],log=TRUE)];

        setkeyv(MET,test.subject); setkeyv(data.prep.object,test.subject);
        d2 = data.prep.object[MET];

        setkey(d2,Source_allele);
        d2[.(1),loglik:=log(h*S1.prop_1 + (1-h)*S0.prop_1)];
        d2[.(0),loglik:=log(h*S1.prop_0 + (1-h)*S0.prop_0)];
        setkey(d2,loglik);
        d2[.(-Inf),loglik:= -50000];

        if(return.likmeans==TRUE){
            setkey(d2,index);
            d3[k>=(burnin+1),loglik:=d2[,loglik]][k==(burnin+1),
                postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
                postvar.loglik.1:=0][k==(burnin+1),
                postmean.loglik.1:=loglik];

        d3[k>(burnin+1),
            postmean.loglik:= postmean.loglik.1 + 
            (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
            postmean.exp.loglik:= postmean.exp.loglik.1 + 
            (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
            postvar.loglik:= postvar.loglik.1 + 
            (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
            postvar.true.loglik:= postvar.loglik/(k-(burnin+1))][k>(burnin+1),
            postmean.loglik.1:= postmean.loglik][k>(burnin+1),
            postvar.loglik.1:= postvar.loglik][k>(burnin+1),
            postmean.exp.loglik.1:= postmean.exp.loglik];
                               };

        setkeyv(d2,c(test.subject,"loglik"));
        MET[,LL:=prior.h + d2[,sum(loglik, na.rm=T),by=test.subject]$V1][,DF:=LL-LF][,
            UF:=runif(.N)];

        setkey(MET,DF);
        MET[DF>10,DF:=10]; MET[,TEMPF:=exp(DF)][,DF0:=DF>0];
        setkey(MET,UF,TEMPF);
        MET[,UFTF:=UF<TEMPF];
        setkey(MET,DF0,UFTF);
        MET[!.(FALSE,FALSE),LF:=LL][!.(FALSE,FALSE),accept := 1];
        MET[.(FALSE,FALSE),h:=h1][.(FALSE,FALSE),accept := 0];

        MET[k==as.integer(burnin*0.4)+1,
	        meanh.1:=h
			][k==as.integer(burnin*0.4)+1,
			var1.1:=0];

        MET[k>as.integer(burnin*0.4)+1 & k<burnin+1, 
            meanh:= meanh.1 + (h - meanh.1)/(k-as.integer(burnin*0.4))
            ][k>as.integer(burnin*0.4)+1 & k<burnin+1,
            var1:= var1.1 + (h - meanh.1)*(h - meanh)
			][k>as.integer(burnin*0.4)+1 & k<burnin+1,
            varh:= var1/(k-(as.integer(burnin*0.4)+1))
			][k>as.integer(burnin*0.4)+1 & k<burnin+1,
            meanh.1:=meanh
			][k>as.integer(burnin*0.4)+1 & k<burnin+1,
			var1.1:=var1];

        MET[k==(burnin+1),
		    postmean.h.1:= h
			][k==(burnin+1),
            postvar.h.1:=0
			][k==(burnin+1),
			indM:= varh*m];

##NEW##
MET[k==101,
    indM:= init.var2];
#######

        MET[k>(burnin+1),
            postmean.h:= postmean.h.1 + (h - postmean.h.1)/(k-burnin)
			][k>(burnin+1),
            postvar.h:= postvar.h.1 + (h - postmean.h.1)*(h - postmean.h)
			][k>(burnin+1),
            postvar.true.h:= postvar.h/(k-(burnin+1))
			][k>(burnin+1),
            postmean.h.1:= postmean.h
			][k>(burnin+1),
            postvar.h.1:= postvar.h];

        if(!is.null(plot.ind)){
            setkeyv(MET,test.subject);

            if(k==1){
			    plot(MET[plot.ind[1],h]~k, xlim=c(1,nitt), ylim=c(0,1), 
                    type="p", cex=0.5, col=plot.col[1], xlab="Iteration",
                    ylab="Hybrid index");
                abline(h=c(0,0.5,1),v=c(burnin,nitt),col="grey",lty=2);

                if(length(plot.ind)>1){
				    points(MET[plot.ind[2],h]~k, cex=0.5, col=plot.col[2]);
                                      };
                if(length(plot.ind)>2){
				    points(MET[plot.ind[3],h]~k, cex=0.5, col=plot.col[3]);
                                      };
                if(length(plot.ind)>3){
				    points(MET[plot.ind[4],h]~k, cex=0.5, col=plot.col[4]);
                                      };
                if(length(plot.ind)>4){
				    points(MET[plot.ind[5],h]~k, cex=0.5, col=plot.col[5]);
                                      };
                if(length(plot.ind)>5){
				    points(MET[plot.ind[6],h]~k, cex=0.5, col=plot.col[6]);
                                      };
                    };

            if(k>1){
			    points(MET[plot.ind[1],h]~k, cex=0.5, col=plot.col[1]);

                if(length(plot.ind)>1){
				    points(MET[plot.ind[2],h]~k, cex=0.5, col=plot.col[2]);
                                      };
                if(length(plot.ind)>2){
				    points(MET[plot.ind[3],h]~k, cex=0.5, col=plot.col[3]);
                                      };
                if(length(plot.ind)>3){
				    points(MET[plot.ind[4],h]~k, cex=0.5, col=plot.col[4]);
                                      };
                if(length(plot.ind)>4){
				    points(MET[plot.ind[5],h]~k, cex=0.5, col=plot.col[5]);
                                      };
                if(length(plot.ind)>5){
				    points(MET[plot.ind[6],h]~k, cex=0.5, col=plot.col[6]);
                                      };
                   }
                              };

        MET[k>1, Zk := r*Zk.1 + accept][k>1, Wk := r*Wk.1 + 1][,A := Zk/Wk];
        MET[k==as.integer(burnin*0.4)+1, indM:=indM*(q^(A-Astar))];
        MET[k>1, Zk.1 := Zk][k>1, Wk.1 := Wk];

        setkeyv(MET,test.subject);
        MET[!fix.subject,h1:=h][!fix.subject,h:=truncnorm::rtruncnorm(.N,mean=h,sd=sqrt(indM),a=0,b=1)];
                    };

    setkeyv(MET,test.subject);
    MET[!fix.subject,
        post.ab:= ((postmean.h*(1 - postmean.h))/postvar.true.h) - 1
		][,
        beta_shape1:= post.ab*postmean.h
		][,
        beta_shape2:= post.ab*(1 - postmean.h)
		][beta_shape1>1&beta_shape2>1,
        h_posterior_mode:= (beta_shape1 - 1)/(beta_shape1 + beta_shape2 - 2)
        ][beta_shape1<1&beta_shape2>=1,
        h_posterior_mode:= 0
        ][beta_shape1>=1&beta_shape2<1,
        h_posterior_mode:= 1
		][!fix.subject,
        h_cred_int_lower:= qbeta(0.025,shape1 = beta_shape1,shape2 = beta_shape2),by=test.subject
		][,
        h_cred_int_upper:= qbeta(0.975,shape1 = beta_shape1,shape2 = beta_shape2),by=test.subject];

    MET[!fix.subject,
	    h_posterior_mode:=round(h_posterior_mode,
        digits = nchar(length(unique(data.prep.object[,locus])), type = "bytes"))
		][!fix.subject,
        beta_shape1:=round(beta_shape1,
        digits = nchar(length(unique(data.prep.object[,locus])), type = "bytes"))
		][!fix.subject,
        beta_shape2:=round(beta_shape2,
        digits = nchar(length(unique(data.prep.object[,locus])), type = "bytes"))
		][!fix.subject,
        h_cred_int_lower:=round(h_cred_int_lower,
        digits = nchar(length(unique(data.prep.object[,locus])), type = "bytes"))
		][!fix.subject,
        h_cred_int_upper:=round(h_cred_int_upper,
        digits = nchar(length(unique(data.prep.object[,locus])), type = "bytes"))];

    MET[fix.subject,h_posterior_mode:=fix.value];

    MET[,LL:=NULL][,LF:=NULL][,DF:=NULL][,UF:=NULL][,TEMPF:=NULL][,h:=NULL][,
    h1:=NULL][,indM:=NULL][,m:=NULL][,accept:=NULL][,r:=NULL][,Zk:=NULL][,
    Zk.1:=NULL][,Wk:=NULL][,Wk.1:=NULL][,A:=NULL][,Astar:=NULL][,q:=NULL][,
    post.ab:=NULL][,prior.h:=NULL][,UFTF:=NULL][,meanh.1:=NULL][,var1.1:=NULL][,
    postmean.h.1:=NULL][,postvar.h:=NULL][,postvar.h.1:=NULL][,varh:=NULL][,
    var1:=NULL][,meanh:=NULL][,DF0:=NULL];

    setnames(MET,c("postmean.h","postvar.true.h"),c("beta_mean","beta_var"));

    setcolorder(MET,c("Source",read.data.precols.red,"h_posterior_mode",
        "h_cred_int_lower","h_cred_int_upper","beta_mean","beta_var",
        "beta_shape1","beta_shape2","npar_compare"));

    output=list();
        output$hi=MET;
	output$test.subject=test.subject;

    if(return.likmeans==TRUE){
        d3[,index:=NULL][,loglik:=NULL][,postmean.exp.loglik.1:=NULL][,
            postmean.loglik.1:=NULL][,postvar.loglik:=NULL][,postvar.loglik.1:=NULL];
        setnames(d3,"postvar.true.loglik","postvar.loglik");
        output$likmeans=d3;
                             };

    cat(paste("Done"),fill=1);

	return(output)
               }

