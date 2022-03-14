#' Bayesian hybrid index estimation.
#'
#' @param data.prep.object Name of the \code{data.prep} object produced by the \code{data.prep} function.
#' @param read.data.precols A \code{precols} object produced by \code{read.data}, or a custom character vector with 
#'   names of all pre-marker columns in the \code{data}.
#' @param test.subject Character string. Name of the field identifying subjects for hybrid index estimation. Default \dQuote{INDLABEL}.
#' @param test.subject.compare Character string. Name of another field in the data analysis table, such as \sQuote{POPID}, with which 
#'   you wish to compare the fit of the current test subject using either AIC or waic. Used to compare model fits at differing resolutions, 
#'   such as individual-level versus population-level, with the subsequent comparison being made at the lower resolution. Only needs to be 
#'   declared when running the higher resolution model, and produces a warning if the declared field is at a higher resolution than the test subject. 
#'   A second run of \command{esth} will then be needed, with \code{test.subject.compare} in the current model declared as \code{test.subject}. 
#'   The two models can then be compared using either or both of the downstream functions \command{compare.models} (for waic) or \command{calc_AIC}.
#' @param include.Source Logical. Whether to estimate hybrid indices for test subjects declared as parental reference. Default is \code{TRUE}. 
#'   Ignored if there are no parental reference samples.
#' @param return.likmeans Logical. Whether to return a table of likelihood calculations. Required for \code{compare.models}. Default is \code{FALSE}.
#' @param fix.subject Logical default, or a character vector. Which, if any, test subjects to fix at a specific hybrid index value. If set \code{TRUE}, 
#'   the entry for \code{fix.value} should be a single value to be applied to all test subjects. Alternatively, a character vector including the names 
#'   of all or a subset of test subjects can be entered, in which case the same set of test subjects can be assigned unique fixed values using 
#'   \code{fix.value}. Potentially useful for model comparison with \code{compare.models}.
#' @param fix.value Logical default, otherwise a numeric vector of length 1 or the same length and order as the character vector declared using 
#'   \code{fix.subject}. The fixed hybrid index value(s). Note that if all test subjects are fixed, only \code{nitt=2} and \code{burnin=0} are required.
#' @param plot.ind Character vector. List of up to 6 test subjects for which you wish to plot the posterior distribution values in real time. 
#'   Default \code{NULL}.
#' @param plot.col Character vector. List of up to 6 colours to plot the \code{plot.ind} test subjects. Default \code{NULL}.
#' @param nitt Numeric. The total number of MCMC iterations including burnin. At least 2000, preferably 3000, is recommended.
#' @param burnin Numeric. The number of burnin MCMC iterations. Typically, 1000 is sufficient.
#' @param init.var Numeric. Starting value for the variance of the parameter proposal distribution, applied for the first 100 iterations. Default 
#'   is 0.002, which works well with tested data sets.
#' @param init.var2 Numeric. Variance of the proposal distribution to be applied after the first 100 iterations, by which time the mcmc should be 
#'   sampling the posterior distribution. The default entry is \code{NULL}, in which case the value is calculated internally using the mean 
#'   number of allele copies per test subject. The internal formula is: (1/(mean(N allele copies per test subject)/2))/10, which works well with 
#'   a big variety of tested numbers of loci per test subject, up to several million.
#' @param prior Numeric vector. Values of \code{shape1} and \code{shape2} parameters for the beta-distributed prior. Default is 
#'   \code{shape1=shape2=0.5} (Jeffrey's prior). Doesn't currently allow a unique prior for each test subject. If set to \code{c(0,0)}, the prior 
#'   is excluded from the posterior likelihood calculation, in which case the best posterior estimate is the mean (\sQuote{beta_mean}) rather than the mode 
#'   (\sQuote{h_posterior_mode}), the latter of which is the best estimate in the presence of a prior.
#' @param print.k Numeric. The iteration is printed to screen at multiples of this number. Default is \code{50}.
#' @details \code{esth} estimates hybrid index using the likelihood formulae of Buerkle (2005), with the addition of a
#'   prior. The default is Jeffrey's prior (beta distribution with \code{shape1=shape2=0.5}), which testing suggests 
#'   is an improvement over a uniform prior (\code{shape1=shape2=1}).
#'
#' \code{esth} excludes any pre-marker columns from the output with finer resolution than the declared \code{test.subject}, 
#'   and hence has the same number of rows as the number of unique \code{test.subject} values.
#'
#' Set \code{return.likmeans=TRUE} if you intend to carry out model comparison based on the widely applicable information criterion (waic).
#' @return list containing \code{hi}, \code{test.subject}, \code{init.var2}, and \code{likmeans} (optional). \code{init.var2} is numeric and 
#'   is included in case the user wants to try a different, possibly smaller, value, for test subjects for which mcmc failed. \code{likmeans} 
#'   contains the information needed to calculate waic. \code{hi} is a \code{data.table} and \code{data.frame} containing all pre-marker 
#'   columns that do not have finer resolution than the declared \code{test.subject}, as well as the following fields:
#'   \item{Source}{whether the sample is from source 0 ('S0'), source 1 ('S1') or is a test individual ('TEST').}
#'   \item{h_posterior_mode}{the hybrid index estimate (mode of the beta-distributed posterior).}
#'   \item{h_cred_int_lower}{lower 95 percent credible interval (2.5 percent quantile of the posterior beta distribution).}
#'   \item{h_cred_int_upper}{upper 95 percent credible interval (97.5 percent quantile of the posterior beta distribution).}
#'   \item{beta_mean}{mean of the posterior beta distribution (included for convenience).}
#'   \item{beta_var}{variance of the posterior beta distribution (included for convenience).}
#'   \item{beta_shape1}{beta shape1 parameter estimate for the posterior beta distribution (included for convenience).}
#'   \item{beta_shape2}{beta shape2 parameter estimate for the posterior beta distribution (included for convenience).}
#'   \item{npar_compare}{see compare.models. Absent unless test.subject.compare is set.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @references
#'   Buerkle, C. A. (2005). Maximum likelihood estimation of a hybrid index based on molecular markers. 
#'   Molecular Ecology Notes, 5(3), 684-687.
#' @examples
#'
#' \dontrun{
#' #Use the sparrow data with loci filtered by parental allele frequency credible interval overlap.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569)
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE)
#' 
#' ###############################################
#' #In the presence of parental reference samples#
#' ###############################################
#' hindlabel=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  include.Source=TRUE,	         #Leave at default TRUE if you want hybrid indices for the parental reference individuals#
#'  plot.ind = c("PD07-254","PD07-160","PD07-159","PI07-110","PI08-498","PH08-285"),  #posterior sampling for up to 6 test subjects can be plotted in real time#
#'   plot.col = c("blue","green","cyan","purple","magenta","red"),
#' nitt=3000,                                                  #Testing suggests nitt=3000,burnin=1000 are plenty for accurate posterior estimation#
#'  burnin=1000
#' )
#' 
#' #There will be a warning if any NAs are produced for the best posterior estimate (h_posterior_mode).
#'
#' ###################################################
#' #esth in the absence of parental reference samples#
#' ###################################################
#' 
#' #This requires a marker.info file with correct column names and containing parental allele frequencies. 
#' #First create this file.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",
#'  precols=dat$precols,AF.CIoverlap = FALSE,return.locus.table = TRUE);
#' sparrowMarkers=prepdata$locus.data[,.(locus,S0_allele,S1_allele,S0.prop_0,S1.prop_0)];
#' setnames(sparrowMarkers,c("S0_allele","S1_allele","S0.prop_0","S1.prop_0"),c("refAllele","alternateAllele","S0.prop_r","S1.prop_r"));
#' #Save it to the working directory.
#' fwrite(sparrowMarkers,"sparrowMarkers.csv")#By default data.table's "fwrite" saves as a comma-separated file, whatever the file extension#
#' 
#' #Now run the workflow without specifying parental reference samples.
#' dat=read.data(file="RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,NUMINDS=569,MISSINGVAL=NA,marker.info.file="sparrowMarkers.csv",sourceAbsent=TRUE);
#' 
#' #Run data.prep filtering for allele frequency difference.
#' prepdata= data.prep(
#'  data=dat$data,
#'  loci=dat$loci,
#'  sourceAbsent = TRUE,        #***Required when no parental reference samples are declared***#
#'  marker.info=dat$marker.info,
#'  alleles=dat$alleles,
#'  precols=dat$precols,
#'  min.diff = 0.2
#' )
#' 
#' #esth can now be run as above.
#' hindlabelx=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,nitt=3000,burnin=1000)
#' 
#' #################################################################
#' #Estimating a hybrid index per locus across all TEST individuals#
#' #################################################################
#' 
#' #Also known as locus-specific ancestry.
#' 
#' #Use a dataset including parental reference samples, and then run per-locus hybrid index estimation on the 'TEST' 
#' #individuals only (i.e. excluding S0 and S1). No plots this time.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);
#' hindlabel=esth(
#'  data.prep.object=prepdata$data.prep,
#'  read.data.precols=dat$precols,
#'  test.subject="locus",             #Switch from estimating h-index per individual (INDLABEL, the default) to per locus#
#'  include.Source=FALSE,	          #Exclude the parental reference individuals (which are included by default) to get the h-index per locus for TEST individuals only#
#'  nitt=3000,
#'  burnin=1000
#' )
#' 
#' ##############################################################################
#' #Running esth on one of the files produced by the 'split_data_prep' function##
#' ##############################################################################
#' 
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",precols=dat$precols,AF.CIoverlap = FALSE);
#' 
#' #Produce one data analysis file per individual and save to WD.
#' split_data_prep(data.prep.object=prepdata$data.prep,splitBy="INDLABEL",keepN=1);
#' 
#' #Read in the first sub-file and estimate hybrid index.
#' prep=fread("prepdata_INDLABEL_1.csv");
#' hindlabel=esth(data.prep.object=prep,read.data.precols=dat$precols,nitt=3000,burnin=1000)
#' }
#' @export
esth= function(data.prep.object,
read.data.precols,
    test.subject="INDLABEL",test.subject.compare=NULL,
    include.Source = TRUE,
    return.likmeans = FALSE,fix.subject = FALSE,fix.value = FALSE,
    plot.ind = NULL,plot.col = NULL,nitt,burnin,
	init.var=0.002,
    init.var2=NULL,
    prior=c(0.5,0.5),
    print.k=50){
	
#Just to prevent any errors when people use a different test.subject#
    if(is.null(test.subject.compare)==FALSE){	
	    if(uniqueN(data.prep.object[,test.subject,with=FALSE]) < uniqueN(data.prep.object[,test.subject.compare,with=FALSE])){
	        test.subject.compare=NULL;
			warning("If the resolution (number of unique values) of test.subject is lower than that of the intended model comparison test.subject, leave 'test.subject.compare' as the default, NULL");
	                                                                                                                         };
											};
	
#If init.var2 is set by the user use that value, otherwise set it internally using the mean number of allele copies per test subject#

    if(is.null(init.var2)){	
	    #if(is.null(N_allele_copies)){
		 #   stop("The average number of allele copies per test subject must be declared in 'N_allele_copies' in order to optimise the MCMC.");
			#		      };
		init.var2=(1/(mean(data.prep.object[,.N,by=test.subject]$N)/2))/10;				  
	                      };						 

    if(return.likmeans==TRUE){
        #data.prep.object[,index:=seq(.N)];
        setkey(data.prep.object,index);
        d3=data.table(index=data.prep.object[,index]);
        setkey(d3,index);
        d3=d3[data.prep.object[,c("index",read.data.precols,test.subject,
            "locus","Source_allele"),with=F]];
		setkey(d3,index);
        d3[,loglik:=0][,postmean.exp.loglik:=0][,postmean.exp.loglik.1:=0][,
            postmean.loglik:=0][,postmean.loglik.1:=0][,
            postvar.true.loglik:=0][,postvar.loglik:=0][,postvar.loglik.1:=0][,
			postmean.loglik.sqrd:=0][,postmean.loglik.sqrd.1:=0];#*******NEW 10 DEC 2021*****#;
                           };

        if(!include.Source){
            setkey(data.prep.object,Source);
            data.prep.object = data.prep.object[Source == "TEST"];
                           };
										 
    read.data.precols.red=test.subject;#character(0);#New edit 14 Oct 2021#
    for(i in 1:length(read.data.precols)){
        if(uniqueN(data.prep.object[,c(test.subject,read.data.precols[i]),with=F]) <= 
            uniqueN(data.prep.object[,test.subject,with=F])){
                read.data.precols.red=c(read.data.precols.red,read.data.precols[i]);
                                                      };
                                         };
					 
	if(is.null(test.subject.compare)==FALSE){
        read.data.precols.red=unique(c(read.data.precols.red,test.subject.compare));#New line 27 Feb 2022#
                                            }else{
		read.data.precols.red=unique(read.data.precols.red);
		                                         };		

	if(uniqueN(data.prep.object[,c("Source",read.data.precols.red),with=F])==uniqueN(data.prep.object[,c(read.data.precols.red),with=F])){
        MET = unique(data.prep.object[,c("Source",read.data.precols.red),with=F])[,
            LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
            TEMPF:=exp(DF)][,h:=runif(.N)][,h1:=h][,
            indM:=init.var][,m:=(2.4/sqrt(1))^2][,
            accept := 0][,r:=0.99][,
            Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
            Astar:=0.44][,q:=2000];
		                                                                                                                                  }else{
        MET = unique(data.prep.object[,c(read.data.precols.red),with=F])[,#I can't refer to the 'Source' column when it's absent, i.e. in the absence of parental reference samples#
            LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
            TEMPF:=exp(DF)][,h:=runif(.N)][,h1:=h][,
            indM:=init.var][,m:=(2.4/sqrt(1))^2][,
            accept := 0][,r:=0.99][,
            Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
            Astar:=0.44][,q:=2000];								 
								 

#***NEW 27 FEB 2022*********##########
    if(fix.subject==TRUE){
        if(length(fix.value) > 1){
		    stop("fix.value should be a single number when fix.subject=TRUE. If you wish to assign a unique fixed value to each test subject, the entry for fix.subject should be a character vector of the same length and order as the list of values entered in fix.value.");
	                             };
        MET[,h:=fix.value][,h1:=h];
                         };								                                                                                                                };

    if(is.character(fix.subject)){
	    if(length(fix.subject)!=length(fix.value)){
		    stop("The values in fix.value should be the same length and order as the test subjects identified in fix.subject.");
			                                      };
		fixed_values=data.table(subjects=fix.subject,values=fix.value);
        setkeyv(MET,test.subject);
		setkey(fixed_values,subjects);
		MET=fixed_values[MET];
		setnames(MET,"subjects",test.subject);
        MET[is.na(values)==FALSE,h:=values][is.na(values)==FALSE,h1:=h];
		MET[,values:=NULL];
                                };
#***NEW 27 FEB 2022*********##########

    if(is.null(test.subject.compare)==FALSE){
        MET[,npar_compare:=.N,by=test.subject.compare];
                                        }else{
		MET[,npar_compare:=NA];
		                                     };

    for(k in 1:nitt){
        if(k==1){cat(paste("Estimating hybrid index..."),fill=1); flush.console();
		        };

        #if(tester::is_multiple(k,print.k)){cat(paste("\t","\t","\t","Iteration",k,"; acceptance=",unique(MET[,A]),"; indM",unique(MET[,indM])),fill=1);#For MCMC testing only# 
        if(tester::is_multiple(k,print.k)){
		    cat(paste("\t","\t","\t","Iteration",k),fill=1);				
            flush.console();     
			                              };

        setkey(MET,h);
        MET[,prior.h:=dbeta(h,shape1=prior[1],shape2=prior[2],log=TRUE)];

        setkeyv(MET,test.subject); setkeyv(data.prep.object,test.subject);
        d2 = data.prep.object[MET];

#Edited to only use S0.prop_1 and S1.prop_1#
        setkey(d2,Source_allele);
        d2[.(1),loglik:=log(h*S1.prop_1 + (1-h)*S0.prop_1)];
        #d2[.(0),loglik:=log(h*S1.prop_0 + (1-h)*S0.prop_0)];
		d2[.(0),loglik:=log(h*(1 - S1.prop_1) + (1-h)*(1 - S0.prop_1))];
		
        setkey(d2,loglik);
        d2[.(-Inf),loglik:= -50000];

        if(return.likmeans==TRUE){
            setkey(d2,index);
            d3[k>=(burnin+1),loglik:=d2[,loglik]][k==(burnin+1),
				postmean.loglik.sqrd.1:=loglik^2][k==(burnin+1),        #****NEW 10 DEC 2021****#
                postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
                postvar.loglik.1:=0][k==(burnin+1),
                postmean.loglik.1:=loglik];

            d3[k>(burnin+1),
                postmean.loglik:= postmean.loglik.1 + 
                (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
                    postmean.loglik.sqrd:=postmean.loglik.sqrd.1 + #****NEW 10 DEC 2021****#
                    (loglik^2 - postmean.loglik.sqrd.1)/(k-burnin)][k>(burnin+1),#****NEW 10 DEC 2021****#
                postmean.exp.loglik:= postmean.exp.loglik.1 + 
                (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
                postvar.loglik:= postvar.loglik.1 + 
                (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
                postvar.true.loglik:= postvar.loglik/(k-(burnin+1))][k>(burnin+1),
                postmean.loglik.1:= postmean.loglik][k>(burnin+1),
                    postmean.loglik.sqrd.1:=postmean.loglik.sqrd][k>(burnin+1),#****NEW 10 DEC 2021****#
                postvar.loglik.1:= postvar.loglik][k>(burnin+1),
                postmean.exp.loglik.1:= postmean.exp.loglik];
                               };

        setkeyv(d2,c(test.subject,"loglik"));
        MET[sum(prior)>0,#********NEW 25 FEB 2022********#
		    LL:=prior.h + d2[,sum(loglik, na.rm=T),by=test.subject]$V1][sum(prior)==0,#********NEW 25 FEB 2022********#
		    LL:=d2[,sum(loglik, na.rm=T),by=test.subject]$V1][,
		    DF:=LL-LF][,
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

##Once the MCMC starts sampling the posterior, change the proposal variance according to the number of allele copies##
        MET[k==101,
            indM:= init.var2];
#############################################################################################################

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

	if(uniqueN(data.prep.object[,c("Source",read.data.precols.red),with=F])==uniqueN(data.prep.object[,c(read.data.precols.red),with=F])){
        setcolorder(MET,c("Source",read.data.precols.red,"h_posterior_mode",
            "h_cred_int_lower","h_cred_int_upper","beta_mean","beta_var",
            "beta_shape1","beta_shape2","npar_compare"));
		                                                                                                                                 }else{
        setcolorder(MET,c(read.data.precols.red,"h_posterior_mode",#I can't refer to the 'Source' column when it's absent, i.e. in the absence of parental reference samples#
            "h_cred_int_lower","h_cred_int_upper","beta_mean","beta_var",
            "beta_shape1","beta_shape2","npar_compare"));		
		                                                                                                                                      };
									   
	if(is.null(test.subject.compare)){
	    MET=MET[,-"npar_compare"]
	                                 };
									 
	if(MET[is.na(h_posterior_mode),.N]>0){
	    warning("Some h_posterior_mode estimates returned NA, suggesting failure to properly estimate the posterior for those test subjects. Try increasing the burnin and/or reducing init.var2 (value is in the output from the ggcline run)")
	                                     };

    output=list();
        output$hi=MET;
	    output$test.subject=test.subject;
		output$init.var2=init.var2;

    if(return.likmeans==TRUE){
        d3[,loglik:=NULL][,postmean.exp.loglik.1:=NULL][,postmean.loglik.sqrd.1:=NULL][,#****NEW 10 DEC 2021***#
            postmean.loglik.1:=NULL][,postvar.loglik:=NULL][,postvar.loglik.1:=NULL];
        setnames(d3,"postvar.true.loglik","postvar.loglik");
		setkey(d3,index);
        output$likmeans=d3;
                             };

    cat(paste("Done"),fill=1);

	return(output)
               }

