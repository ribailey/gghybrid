#' Estimate genomic cline parameters, credible intervals and p values.
#'
#' @param data.prep.object Name of the \code{data.prep} object produced by the \code{data.prep} function.
#' @param esth.object Name of the \code{esth} object.
#' @param esth.colname Name of the column in the \code{hi} object within \code{esth} containing the best posterior estimates for 
#'   hybrid index. Default is \sQuote{h_posterior_mode}. If hybrid index was estimated in the absence of a prior, use 
#'   \sQuote{beta_mean} instead.
#' @param test.subject Character string. Name of the field identifying subjects for genomic cline estimation. Default \dQuote{locus}.
#' @param poolv Logical. Whether to pool the rate parameter \code{v} across all test subjects. Default \code{FALSE}.
#' @param poolcentre Logical. Whether to pool the \code{centre} parameter across all test subjects. Default \code{FALSE}.
#' @param include.Source Logical. Whether to include the parental reference individuals in genomic cline 
#'   estimation. Default \code{FALSE}.
#' @param return.likmeans Logical. Whether to return a table of likelihood calculations. Required for applying the widely applicable 
#'   information criterion to model comparison using \code{compare.models}. Default is \code{FALSE}.
#' @param read.data.precols A \code{precols} object produced by \code{read.data}, or a custom character vector with 
#'   names of all pre-marker columns in the \code{data}.
#' @param fix.subject.v Logical default, otherwise a character vector. List of test subjects (normally locus names) 
#'   for which you wish to fix rather than estimate the rate parameter \code{v}. Used for model 
#'   comparison with \code{compare.models}. Set to \code{TRUE} to set all test subjects to the same value. Default is \code{FALSE}.
#' @param fix.value.v Logical default, otherwise a numeric vector of length 1 or the same length as \code{fix.subject.v}. 
#'   The fixed \code{v} value(s) for the test subject(s) listed in \code{fix.subject.v}.
#' @param fix.subject.centre Logical default, otherwise a character vector. List of test subjects (normally locus names) 
#'   for which you wish to fix rather than estimate the \code{centre} parameter . Used for model 
#'   comparison with \code{compare.models}. Set to \code{TRUE} to set all test subjects to the same value. Default is \code{FALSE}.
#' @param fix.value.centre Logical default, otherwise a numeric vector of length 1 or the same length as \code{fix.subject.centre}. 
#'   The fixed \code{centre} value(s) for the test subject(s) listed in \code{fix.subject.centre}.
#' @param plot.test.subject Character vector. List of up to 3 test subjects for which you wish to plot the posterior distribution values in real time. 
#'   Rate parameter \code{v} is plotted with open circles, and logit(\code{centre}) with plus signs, both the same colour for the same 
#'   test subject. Default \code{NULL}.
#' @param plot.col Character vector. List of up to 3 colours to plot the \code{plot.test.subject} test subjects. Default \code{NULL}.
#' @param plot.ylim Numeric vector. Upper and lower plot limits of the y axis. Default \code{c(-5,5)}.
#' @param plot.pch.v.centre Numeric vector. Size of data points for the two parameters. Default \code{c(1,3)}.
#' @param prior.logitcentre Numeric vector. Mean and sd of the normally distributed prior for logit(\code{centre}). Default \code{c(0,sqrt(50))}.
#' @param prior.logv Numeric vector. Mean and sd of the normally distributed prior for ln(\code{v}). Default \code{c(0,sqrt(10))}.
#' @param nitt The total number of MCMC iterations including burnin. 5,000 is recommended, 10,000 if one parameter is pooled.
#' @param burnin Numeric. The number of burnin MCMC iterations. 2,000 is recommended, 5,000 if one parameter is pooled.
#' @param start.v Numeric vector. Optional starting \code{v} values for the MCMC, which are otherwise sampled randomly. Default \code{NULL}.
#' @param start.centre Numeric vector. Optional starting \code{centre} values for the MCMC, which are otherwise sampled randomly. Default \code{NULL}.
#' @param init.var.v Numeric. Optional starting variance of the log(\code{v}) parameter proposal distribution. Internal value works 
#'   well on tested data sets. Default \code{NULL}.
#' @param init.var.centre Optional starting variance of the logit(\code{centre}) parameter proposal distribution. Internal value works 
#'   well on tested data sets. Default \code{NULL}.
#' @param init.cov.vcentre Optional starting covariance of the log(\code{v}) and logit(\code{centre}) multivariate normal proposal distribution. 
#'   Internal value (0) works well on tested data sets. Default \code{NULL}.
#' @param print.k The iteration is printed to screen on multiples of this number. Default is \code{50}.
#' @details \code{ggcline} assumes the posterior of \code{ln(v)} is normally distributed, and the results for this parameter are 
#'   the inverse logs of the estimated mean and upper and lower 95 percent credible intervals of ln(\code{v}). \code{v} is always positive and higher 
#'   values for \code{v} indicate steeper clines, with \code{v=1} being the null value. \code{centre} represents the hybrid index
#'   at which allele frequencies are half way between those of the parents. It ranges between 0 and 1 and the null value is 0.5. It 
#'   is estimated as logit(\code{centre}), and the posterior mean and credible intervals are inverse logit-transformed back to the 
#'   original scale.
#'
#' While Fitzpatrick (2013) describes parameters, \code{v} (the relative cline slope) and \code{u} (the relative cline position), 
#'   the parameter \code{u} is difficult to interpret as its range scales to \code{v}. Noting that for the hybrid index itself,
#'   \code{u=0} and \code{v=1}, the cline \code{centre} (hybrid index value for which allele frequencies are half way between those of
#'   the parents: \code{m} in Fitzpatrick's notation) for individual loci has the relationship \code{logit(centre)=u/v}. \code{centre}
#'   is easier to interpret, and estimating it rather than \code{u} improves MCMC efficiency; hence I estimate \code{centre} 
#'   rather than \code{u}.
#'
#' If both parameters are fixed to the null or other values, only \code{nitt=2} and \code{burnin=0} are required.
#'
#' For \code{poolv} and \code{poolcentre}, the options are no pooling (\code{FALSE}) or pooling across all 
#'   samples in the data set (\code{TRUE}). So if for example you wish to obtain a single pooled estimate of \code{v} across two transects 
#'   for each locus individually (i.e. not the same parameter estimate across all loci), 
#'   the loci must be run individually by subsetting the \code{data.prep.object} to only include one locus in each run, using \command{split_data_prep}.
#' @return list with at least two components, \code{gc}, a \code{data.table} and \code{data.frame} with the estimated 
#'   genomic cline parameters, and \code{test.subject}, a character scalar. The optional output \code{likmeans} is needed downstream 
#'   when model comparison is carried out with \code{compare.models}.
#'
#'   \code{gc} contains the \code{test.subject} and following fields:
#'   \item{locus}{the locus name, or another chosen test.subject.}
#'   \item{S0.prop_1}{allele frequency of the S1 allele in the S0 source population.}
#'   \item{S1.prop_1}{allele frequency of the S1 allele in the S1 source population.}
#'   \item{exp_mean_log_v}{the best posterior estimate of v (inverse-logged mean of the posterior normal distribution for ln(v)).}
#'   \item{v_lower_95}{v lower 95 percent credible interval (inverse-logged 2.5 percent quantile of the posterior normal distribution for ln(v)).}
#'   \item{v_upper_95}{v upper 95 percent credible interval (inverse-logged 97.5 percent quantile of the posterior normal distribution for ln(v)).}
#'   \item{v_pvalue}{Bayesian p value for v (2*Quantile of the null value (log(1)) given the posterior normal distribution for ln(v)).}
#'   \item{invlogit_mean_logit_centre}{best posterior estimate of centre (inverse logit-transformed mean of the posterior normal distribution for logit(centre)).}
#'   \item{centre_lower_95}{centre lower 95 percent credible interval (inverse logit-transformed 2.5 percent quantile of posterior logit(centre)).}
#'   \item{centre_upper_95}{centre upper 95 percent credible interval (inverse logit-transformed 97.5 percent quantile of posterior logit(centre)).}
#'   \item{centre_pvalue}{Bayesian p value for centre (2*Quantile of the null value (logit(0.5)) given the posterior normal distribution for logit(centre)).}
#'   \item{mean_log_v}{posterior mean v on the normally distributed latent scale (ln(v)).}
#'   \item{var_log_v}{posterior variance of v on the normally distributed latent scale (ln(v)).}
#'   \item{mean_logit_centre}{posterior mean centre on the normally distributed latent scale (logit(centre)).}
#'   \item{var_logit_centre}{posterior variance of centre on the normally distributed latent scale (logit(centre)).}
#'   \item{cov_log_v_logit_centre}{posterior covariance between the multivariate normally distributed ln(v) and (logit(centre)).}
#'   \item{npar}{number of parameters for each test subject.}
#'   \item{sum_npar}{sum number of parameters across all test subjects.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @references
#'   Fitzpatrick, B. M. (2013). Alternative forms for genomic clines. Ecology and evolution, 3(7), 1951-1966. 
#' @examples
#'
#' \dontrun{
#' ###########################################################
#' #Read in and prepare data, and run hybrid index estimation#
#' ###########################################################
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569);
#' prepdata=data.prep(data=dat$data,loci=dat$loci,alleles=dat$alleles,S0=c("Kralove","Oslo"),S1="LesinaSPANISH",
#'  precols=dat$precols,AF.CIoverlap = FALSE);
#' hindlabel=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,include.Source=TRUE,
#'  nitt=3000,burnin=1000);
#' 
#' #################################################################################################################
#' #The full set of arguments for ggcline including all defaults (comments indicate which entries are not defaults)#
#' #################################################################################################################
#' gc=ggcline(
#'   data.prep.object=prepdata$data.prep,    #Needs an entry#
#'   esth.object=hindlabel,                  #Needs an entry#
#'   esth.colname="h_posterior_mode",          #Default. Can be replaced with "beta_mean" from the esth object#
#'   test.subject = "locus",
#'   poolv = FALSE,
#'   poolcentre = FALSE,
#'   include.Source = TRUE,                  #Default is FALSE#
#'   return.likmeans = TRUE,                 #Default is FALSE#
#'   read.data.precols=dat$precols,          #Needs an entry#
#'   fix.subject.v = FALSE,
#'   fix.value.v,
#'   fix.subject.centre = FALSE,
#'   fix.value.centre,
#'   plot.test.subject = c("A2ML1_SNP1","GTF2H2"), #Plots are just to check the MCMC is working#
#'   plot.col = c("orange","cyan"),                #Plots are just to check the MCMC is working#
#'   plot.ylim = c(-3, 5),                         #Plots are just to check the MCMC is working#
#'   plot.pch.v.centre = c(1, 3),                  #Plots are just to check the MCMC is working#
#'   prior.logitcentre = c(0, sqrt(50)),         #Default#
#'   prior.logv = c(0, sqrt(10)),                #Default#
#'   nitt=5000,                                     #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
#'   burnin=2000,                                   #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
#'   start.v = NULL,
#'   start.centre = NULL,
#'   init.var.v = NULL,
#'   init.var.centre = NULL,
#'   init.cov.vcentre = NULL,
#'   print.k = 50
#' )
#' 
#' #############################################################################################
#' #Estimate a single v parameter across a whole chromosome, but a unique centre for each locus#
#' #############################################################################################
#' 
#' #First, add a column to prepdata$data.prep indicating the chromosome for each marker.
#' chrom=fread("markerchr.csv");#This contains only the loci retained in prepdata$data.prep after filtering#
#' setkey(prepdata$data.prep,locus);setkey(chrom,locus);
#' prepdata$data.prep=prepdata$data.prep[chrom];
#' 
#' #Make a data.prep object for the Z chromosome markers only.
#' prepZ=prepdata$data.prep[chrom=="Z"]
#' 
#' #First run ggcline for chrZ-only without pooling, for later model comparison.
#' gc_unpooled_z=ggcline(data.prep.object=prepZ,esth.object=hindlabel,include.Source = TRUE,
#'  return.likmeans = TRUE,               #***Important to set this to TRUE waic for model comparison***#
#'  read.data.precols=dat$precols, nitt=5000,burnin=2000)
#' 
#' #Run again but pooling v across all loci. The MCMC is less efficient with pooling, so using longer nitt and burnin.
#' gc_poolv_Z=ggcline(data.prep.object=prepZ,esth.object=hindlabel,              
#'  poolv = TRUE,                         #***Set to TRUE***#
#'  include.Source = TRUE,
#'  return.likmeans = TRUE,               #***Important to set this to TRUE for model comparison***#
#'  read.data.precols=dat$precols,
#'  nitt=10000,                           #Use a longer post-burnin (5000 in this case)#
#'  burnin=5000                           #Use a longer burnin#
#' )
#' 
#' ########################################################################
#' #Set up model comparisons by fixing one or both parameters (no pooling)#
#' ########################################################################
#' 
#' #Run ggcline with both parameters fixed to their null values. This can then be compared with the first run above (gc).
#' gc_v_centre_fixed=ggcline(data.prep.object=prepdata$data.prep,esth.object=hindlabel,include.Source = TRUE,
#'  return.likmeans = TRUE,               #***Important to set this to TRUE for model comparison***#
#'  read.data.precols=dat$precols,
#'  fix.subject.v = gc$gc$locus,          #***Character vector of locus names, here taken from the first ggcline run above***#
#'  fix.value.v = 1,                      #***The null parameter value on the data scale***#
#'  fix.subject.centre = gc$gc$locus,     #***Vector of locus names***#
#'  fix.value.centre = 0.5,               #***The null parameter value on the data scale***#
#'  nitt=2,                               #MCMC estimation not needed when both parameters are fixed for all test subjects, only nitt=2 is necessary to avoid errors#
#'  burnin=0                              #MCMC estimation not needed when both parameters are fixed for all test subjects#
#' )
#' 
#' #The comparison would be cleaner if we only fixed one parameter. For example, v.
#' gc_v_fixed=ggcline(data.prep.object=prepdata$data.prep,esth.object=hindlabel,include.Source = TRUE,
#'  return.likmeans = TRUE,               #***Important to set this to TRUE for model comparison***#
#'  read.data.precols=dat$precols,
#'  fix.subject.v = gc$gc$locus,          #***Vector of locus names***#
#'  fix.value.v = 1,                      #***The null parameter value on the data scale (could also be fixed to some other value)***#
#'  #fix.subject.centre = gc$gc$locus,    #Not fixing centre this time#
#'  #fix.value.centre = 0.5,              #Not fixing centre this time#
#'  nitt=5000,
#'  burnin=2000
#' )
#' 
#' #This can now be compared with either gc or gc_v_centre_fixed.
#' 
#' ##################################################################
#' #Prepare to compare parameters estimated at differing resolutions#
#' ##################################################################
#' 
#' #Run ggcline using chromosome as test.subject rather than locus. We added the column "chrom" to prepdata$data.prep earlier.
#' gc_by_chr=ggcline(data.prep.object=prepdata$data.prep,esth.object=hindlabel,
#'  test.subject = "chrom",               #***Changed from default***#
#'  include.Source = TRUE,return.likmeans = TRUE,read.data.precols=dat$precols,nitt=5000,burnin=2000)
#' 
#' #The column "chrom" is absent from the 'gc' results object created above, and must first be added to both gc$gc and gc$likmeans.
#' chrom=fread("markerchr.csv");#This contains only the loci retained in prepdata$data.prep#
#' setkey(gc$gc,locus);setkey(chrom,locus);setkey(gc$likmeans,locus);
#' gc$gc=gc$gc[chrom];
#' gc$likmeans=gc$likmeans[chrom];
#' 
#' #Now model comparison can be carried out at the chromosome level, between gc and gc_by_chr.
#' }
#' @export
ggcline=function(data.prep.object,
    esth.object,                           #This can now be a custom object, as long as it has the sub-lists $hi and $test.subject#
	esth.colname="h_posterior_mode",
	test.subject="locus",
    poolv=FALSE,poolcentre=FALSE,include.Source=FALSE,
    return.likmeans=FALSE,read.data.precols,
    fix.subject.v=FALSE,fix.value.v,                       #On the data scale#
    fix.subject.centre=FALSE,fix.value.centre,             #On the data scale#
    plot.test.subject=NULL,plot.col=NULL,
    plot.ylim=c(-5,5),plot.pch.v.centre=c(1,3),
    prior.logitcentre=c(0,sqrt(50)),prior.logv=c(0,sqrt(10)), #Latent scale#
    nitt,burnin,
	start.v=NULL,start.centre=NULL,init.var.v=NULL,
	init.var.centre=NULL,init.cov.vcentre=NULL,
    print.k=50){
	
    if(fix.subject.centre[1]!=FALSE){
	    if(fix.value.centre==0){stop("fix.value.centre cannot be 0 or 1")};
	    if(fix.value.centre==1){stop("fix.value.centre cannot be 0 or 1")};		
		                            };

	if(is.null(esth.object$hi)){
	    stop("The declared esth.object must contain lists named $hi (containing a column with the hybrid index value for each test subject and a column named after the test subject, typically INDLABEL) and $test.subject (naming the test subject for the hybrid index values, typically INDLABEL");
			                             };
										 
	if(is.null(esth.object$test.subject)){
	    stop("The declared esth.object must contain lists named $hi (containing a column with the hybrid index value for each test subject and a column named after the test subject, typically INDLABEL) and $test.subject (naming the test subject for the hybrid index values, typically INDLABEL");
			                             };

    test.subject.esth=esth.object$test.subject;
    esth.object=esth.object$hi;
			
    if(!test.subject.esth%in%names(data.prep.object)){
        stop("The declared esth.object$test.subject must be a column name present in the declared data.prep.object so the tables can be joined");
		                                             };

    if(return.likmeans==TRUE){
        if(!include.Source){
            setkey(data.prep.object,Source);
            data.prep.object=data.prep.object[Source=="TEST"];
            setkey(esth.object,Source);
            esth.object=esth.object[Source=="TEST"];
                           };
        #data.prep.object[,index:=seq(.N)];#No longer needed#
        setkey(data.prep.object,index);
        d3=data.table(index=data.prep.object[,index]);
        setkey(d3,index);
        d3=d3[data.prep.object[,c("index",read.data.precols,test.subject,"locus","Source_allele"),with=F]];
        d3[,loglik:=0][,postmean.exp.loglik:=0][,postmean.exp.loglik.1:=0][,
            postmean.loglik:=0][,postmean.loglik.1:=0][,
            postvar.true.loglik:=0][,postvar.loglik:=0][,postvar.loglik.1:=0][,
			postmean.loglik.sqrd:=0][,postmean.loglik.sqrd.1:=0];#*******NEW 10 DEC 2021*****#
                             };

    if(poolcentre + poolv==0){
        setkeyv(data.prep.object,test.subject.esth);
        setkeyv(esth.object,test.subject.esth);
        MELT.DT=data.prep.object[esth.object];

        if(!include.Source){
            setkey(MELT.DT,Source);MELT.DT=MELT.DT[Source=="TEST"];
                           };
        setnames(MELT.DT,esth.colname,"h");

        MET.gc=unique(data.prep.object[,test.subject,with=F]);
        setkeyv(MET.gc,test.subject);

        MET.gc[,LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
            TEMPF:=exp(DF)][,
            m:=(2.4/sqrt(2))^2][,V1:=-5000000][,
            accept:=0][,r:=0.99][,
            Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
            Astar:=0.3875][,q:=2000][,npar:=2];

        if(fix.subject.v[1]!=FALSE){
            MET.gc[fix.subject.v,m:=(2.4/sqrt(1))^2][fix.subject.v,
                Astar:=0.44][fix.subject.v,npar:=npar - 1];
                                   };
        if(fix.subject.centre[1]!=FALSE){
            MET.gc[fix.subject.centre,m:=(2.4/sqrt(1))^2][fix.subject.centre,
                Astar:=0.44][fix.subject.centre,npar:=npar - 1];
                                        };

        MET.gc[,sum_npar:=sum(npar)];

        if(is.null(start.v)){
            MET.gc[,v:=rtlnorm(.N,0,1,lower=0.02,upper=50)][,v1:=v];
                            }else{
            MET.gc[,v:=start.v][,v1:=v];
                                 };

        if(is.null(start.centre)){
            MET.gc[,centre:=gghybrid::rtnorm(.N,mean=0,sd=1,lower=-10,upper=10)][,centre1:=centre];
                            }else{
            MET.gc[,centre:=qlogis(start.centre)][,centre1:=centre];
                                 };

        if(is.null(init.var.v)){
            MET.gc[,indM11:=0.01];
                               }else{
            MET.gc[,indM11:=init.var.v];
                                    };

        if(is.null(init.var.centre)){
            MET.gc[,indM22:=0.02];
                               }else{
            MET.gc[,indM22:=init.var.centre];
                                    };

        if(is.null(init.cov.vcentre)){
            MET.gc[,indM12:=0];
                                }else{
            MET.gc[,indM12:=init.cov.vcentre];
                                     };

        if(fix.subject.v[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.v,v1:=v];
                                   };

        if(fix.subject.centre[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.centre,centre:=qlogis(fix.value.centre)][fix.subject.centre,centre1:=centre];
                                        };

        test=c(test.subject,"loglik.gc");

        for(k in 1:nitt){
            if(k==1){
                cat(paste("Estimating genomic clines..."),fill=1); 
                flush.console();
                    };

            if(tester::is_multiple(k,print.k)){
                cat(paste("\t","\t","\t","Iteration",k),fill=1); 
                flush.console();
                                      };

            MET.gc[,d.centre:=dnorm(centre,mean=prior.logitcentre[1],sd=prior.logitcentre[2],log=TRUE)][,
                d.v:=dnorm(log(v),mean=prior.logv[1],sd=prior.logv[2],log=TRUE)];

            setkeyv(MELT.DT,test.subject);setkeyv(MET.gc,test.subject);
            MELTED.DT=MELT.DT[MET.gc];
			
			MELTED.DT[,u:=centre*v];

#Edited to only use S0.prop_1 and S1.prop_1#
            setkey(MELTED.DT,Source_allele);
            MELTED.DT[.(1),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_1 + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_1)][.(0),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*(1 - S1.prop_1) + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*(1 - S0.prop_1))];
                #loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_0 + 
                #    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_0)];

            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(-Inf),loglik.gc:=-1000];
            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(NaN),loglik.gc:=-1000];

            if(return.likmeans==TRUE){
                setkey(MELTED.DT,index);
                d3[k>=(burnin+1),loglik:=MELTED.DT[,loglik.gc]][k==(burnin+1),
				    postmean.loglik.sqrd.1:=loglik^2][k==(burnin+1),        #****NEW 10 DEC 2021****#
                    postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
                    postvar.loglik.1:=0][k==(burnin+1),
                    postmean.loglik.1:=loglik];

                d3[k>(burnin+1),
                    postmean.loglik:=postmean.loglik.1 + 
                    (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
                    postmean.loglik.sqrd:=postmean.loglik.sqrd.1 + #****NEW 10 DEC 2021****#
                    (loglik^2 - postmean.loglik.sqrd.1)/(k-burnin)][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postmean.exp.loglik:=postmean.exp.loglik.1 + 
                    (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
                    postvar.loglik:=postvar.loglik.1 + 
                    (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
                    postvar.true.loglik:=postvar.loglik/(k-(burnin+1))][k>(burnin+1),
                    postmean.loglik.1:=postmean.loglik][k>(burnin+1),
                    postmean.loglik.sqrd.1:=postmean.loglik.sqrd][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postvar.loglik.1:=postvar.loglik][k>(burnin+1),
                    postmean.exp.loglik.1:=postmean.exp.loglik];
                                     };

            setkeyv(MELTED.DT,test);
            sum.loglik=MELTED.DT[,sum(loglik.gc, na.rm=T),by=test.subject];
            setkeyv(sum.loglik,test.subject);
            MET.gc[,V1:=NULL];
            MET.gc=MET.gc[sum.loglik];
            MET.gc[,LL:=d.centre + d.v + V1][,
                DF:=LL-LF][,UF:=runif(.N)];

            setkey(MET.gc,DF);
            MET.gc[DF>10,DF:=10][,TEMPF:=exp(DF)][,DF0:=DF>0][,UFTF:=UF<TEMPF];
            setkey(MET.gc,DF0,UFTF);
            MET.gc[!.(FALSE,FALSE),
                LF:=LL][!.(FALSE,FALSE),
                accept:=1][.(FALSE,FALSE),
                v:=v1][.(FALSE,FALSE),
                centre:=centre1][.(FALSE,FALSE),
                accept:=0];

            MET.gc[k==as.integer(burnin*0.2)+1,
                meanv.1:=log(v)][k==as.integer(burnin*0.2)+1,
                meancentre.1:=centre][k==as.integer(burnin*0.2)+1,
                var11.1:=0][k==as.integer(burnin*0.2)+1,
                var12.1:=0][k==as.integer(burnin*0.2)+1,
                var22.1:=0][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv:=meanv.1 + (log(v) - meanv.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre:=meancentre.1 + (centre - meancentre.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11:=var11.1 + (log(v) - meanv.1)*(log(v) - meanv)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.v:=var11/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22:=var22.1 + (centre - meancentre.1)*(centre - meancentre)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.centre:=var22/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12:=var12.1 + ((k-as.integer(burnin*0.2+1))/(k-as.integer(burnin*0.2)))*(log(v) - meanv.1)*
                    (centre - meancentre.1)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                cov.vcentre:=var12/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv.1:=meanv][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre.1:=meancentre][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11.1:=var11][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22.1:=var22][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12.1:=var12][k>as.integer(burnin*0.599) & k<burnin+1,
                indM11:=var.v][k>as.integer(burnin*0.599) & k<burnin+1,
                indM12:=cov.vcentre][k>as.integer(burnin*0.599) & k<burnin+1,
                indM22:=var.centre][k==(burnin+1),
                post.mean.logv.1:=log(v)][k==(burnin+1),
                post.mean.centre.1:=centre][k==(burnin+1),
                post.var.logv.1:=0][k==(burnin+1),
                post.var.centre.1:=0][k==(burnin+1),
                var12.1:=0][k>(burnin+1),
                post.mean.logv:=post.mean.logv.1 + (log(v) - post.mean.logv.1)/(k-burnin)][k>(burnin+1),
                post.mean.centre:=post.mean.centre.1 + (centre - post.mean.centre.1)/(k-burnin)][k>(burnin+1),
                post.var.logv:=post.var.logv.1 + (log(v) - post.mean.logv.1)*(log(v) - post.mean.logv)][k>(burnin+1),
                post.var.true.logv:=post.var.logv/(k-(burnin+1))][k>(burnin+1),
                post.var.centre:=post.var.centre.1 + (centre - post.mean.centre.1)*(centre - post.mean.centre)][k>(burnin+1),
                post.var.true.centre:=post.var.centre/(k-(burnin+1))][k>(burnin+1),  #Changes 4 Nov2021 start here#
                var12:=var12.1 + ((k-(burnin+1))/(k-(burnin)))*(log(v) - post.mean.logv.1)*
                    (centre - post.mean.centre.1)][k>(burnin+1),
                cov.vcentre:=var12/(k-(burnin+1))][k>(burnin+1),
                post.mean.logv.1:=post.mean.logv][k>(burnin+1),
                post.mean.centre.1:=post.mean.centre][k>(burnin+1),
                post.var.logv.1:=post.var.logv][k>(burnin+1),
                post.var.centre.1:=post.var.centre][k>(burnin+1),
                var12.1:=var12][k==burnin+2,
				meanv.1:=NULL][k==burnin+2,
				meancentre.1:=NULL][k==burnin+2,
				var11.1:=NULL][k==burnin+2,
				var22.1:=NULL][k==burnin+2,
				meanv:=NULL][k==burnin+2,
				meancentre:=NULL][k==burnin+2,
				var11:=NULL][k==burnin+2,
				var.v:=NULL][k==burnin+2,
				var22:=NULL][k==burnin+2,
				var.centre:=NULL][k>1 & k<as.integer(burnin*0.4)+2,
                Zk:=r*Zk.1 + accept][k>1 & k<as.integer(burnin*0.4)+2, 
                Wk:=r*Wk.1 + 1][k<as.integer(burnin*0.4)+2,
                A:=Zk/Wk][k==as.integer(burnin*0.4)+1,
                indM11:=indM11*(q^(A-Astar))][k==as.integer(burnin*0.4)+1,
                indM12:=indM12*(q^(A-Astar))][k==as.integer(burnin*0.4)+1,
                indM22:=indM22*(q^(A-Astar))][k>1 & k<as.integer(burnin*0.4)+2, 
                Zk.1:=Zk][k>1 & k<as.integer(burnin*0.4)+2, 
                Wk.1:=Wk][k==burnin+1,
                indM11:=indM11*m][k==burnin+1,
                indM12:=indM12*m][k==burnin+1,
                indM22:=indM22*m];

            if(!is.null(plot.test.subject)){
                setkeyv(MET.gc,test.subject);
                if(k==1){
                    plot(MET.gc[plot.test.subject[1],v]~k, col=plot.col[1],xlim=c(1,nitt), 
                        ylim=plot.ylim,pch=plot.pch.v.centre[1],type="p",ylab="Parameter value");
                    abline(v=c((as.integer(burnin*0.4)+1),(as.integer(burnin*0.6)),
                        burnin,nitt),col="grey",lty=2,lwd=c(1,1,2,2));
                    abline(h=c(0,1),col=c("black","grey"),lty=2);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[2]);
                                                   };
                        };

                if(k>1){
                    points(MET.gc[plot.test.subject[1],v]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[1]);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };
                       };
                                           };

            MET.gc[,v1:=v][,centre1:=centre][,
                c("v","centre"):=rtmvnormDT(mean1=log(v),mean2=centre,M11=indM11,M12=indM12,M22=indM22),by=test.subject];

            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.centre,centre:=qlogis(fix.value.centre)];
                        };
                        };

    if(poolcentre + poolv==1){
        lenv=numeric(0);
        lenu=numeric(0);

        setkeyv(data.prep.object,test.subject.esth);
        setkeyv(esth.object,test.subject.esth);
        MELT.DT=data.prep.object[esth.object];

        if(!include.Source){
            setkey(MELT.DT,Source);MELT.DT=MELT.DT[Source=="TEST"];
                           };
        setnames(MELT.DT,esth.colname,"h");

        MET.gc=unique(data.prep.object[,test.subject,with=F]);
        setkeyv(MET.gc,test.subject);
        if(poolv){lenv=1}else{lenv=nrow(MET.gc)};
        if(poolcentre){lenu=1}else{lenu=nrow(MET.gc)};

        MET.gc[,LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
            TEMPF:=exp(DF)][,
            m:=(2.4/sqrt(lenv+lenu))^2][,V1:=-5000000][,
            accept:=0][,r:=0.95][,
            Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
            Astar:=0.44 - 0.21*(1/5)*(lenv + lenu - 1)][Astar<0.23,
            Astar:=0.23][,q:=20][,npar:=NA][,sum_npar:=lenv + lenu];

        if(fix.subject.v[1]!=FALSE){
            MET.gc[fix.subject.v,m:=(2.4/sqrt(lenu))^2][fix.subject.v,
                Astar:=0.44 - 0.21*(1/5)*(lenu - 1)];
            if(lenv>1){
                MET.gc[,sum_npar:=sum_npar - nrow(MET.gc[fix.subject.v])];
                      }else{
                MET.gc[,sum_npar:=sum_npar - 1];
                           };
                                   };
        if(fix.subject.centre[1]!=FALSE){
            MET.gc[fix.subject.centre,m:=(2.4/sqrt(lenv))^2][fix.subject.centre,
                Astar:=0.44 - 0.21*(1/5)*(lenv - 1)];
            if(lenu>1){
                MET.gc[,sum_npar:=sum_npar - nrow(MET.gc[fix.subject.centre])];
                      }else{
                MET.gc[,sum_npar:=sum_npar - 1];
                           };
                                   };

        LL.all=-5000000;
        LF.all=LL.all;
        DF.all=LL.all-LF.all;
        UF.all=runif(1);
        TEMPF.all=exp(DF.all);
        accept.all=0;

        if(is.null(start.v)){
            if(lenv==1){
                MET.gc[,v:=rtlnorm(1,0,1,lower=0.02,upper=50)][,v1:=v];
                       }else{
                MET.gc[,v:=rtlnorm(.N,0,1,lower=0.02,upper=50)][,v1:=v];
                            };
                            }else{
            MET.gc[,v:=start.v][,v1:=v];
                                 };

        if(is.null(start.centre)){
            if(lenu==1){
                MET.gc[,centre:=gghybrid::rtnorm(1,mean=0,sd=1,lower=-10,upper=10)][,centre1:=centre];
                       }else{
                MET.gc[,centre:=gghybrid::rtnorm(.N,mean=0,sd=1,lower=-10,upper=10)][,centre1:=centre];
                            };
                                 }else{
            MET.gc[,centre:=qlogis(start.centre)][,centre1:=centre];
                                      };

        if(is.null(init.var.v)){
            MET.gc[,indM11:=0.01];
                               }else{
            MET.gc[,indM11:=init.var.v];
                                    };

        if(is.null(init.var.centre)){
            MET.gc[,indM22:=0.02];
                               }else{
            MET.gc[,indM22:=init.var.centre];
                                    };

        if(is.null(init.cov.vcentre)){
            MET.gc[,indM12:=0];
                                }else{
            MET.gc[,indM12:=init.cov.vcentre];
                                     };

        MET.gc[,ML:=-500000][,mlv:=v][,mlcentre:=centre];

        if(fix.subject.v[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.v,v1:=v];
                                   };

        if(fix.subject.centre[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.centre,centre:=qlogis(fix.value.centre)][fix.subject.centre,centre1:=centre];
                                   };

        test=c(test.subject,"loglik.gc");

        for(k in 1:nitt){
            if(k==1){
                cat(paste("Estimating genomic clines..."),fill=1); 
                flush.console();
                    };

            if(tester::is_multiple(k,print.k)){
                cat(paste("\t","\t","\t","Iteration",k),fill=1); 
                flush.console();
                                      };

            MET.gc[,d.centre:=dnorm(centre,mean=prior.logitcentre[1],sd=prior.logitcentre[2],log=TRUE)][,
                d.v:=dnorm(log(v),mean=prior.logv[1],sd=prior.logv[2],log=TRUE)];

            setkeyv(MELT.DT,test.subject);setkeyv(MET.gc,test.subject);
            MELTED.DT=MELT.DT[MET.gc];
			
			MELTED.DT[,u:=centre*v];

#Edited to only use S0.prop_1 and S1.prop_1#
            setkey(MELTED.DT,Source_allele);
            MELTED.DT[.(1),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_1 + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_1)][.(0),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*(1 - S1.prop_1) + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*(1 - S0.prop_1))];
                #loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_0 + 
                #    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_0)];

            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(-Inf),loglik.gc:=-1000];
            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(NaN),loglik.gc:=-1000];

            #if(return.likmeans==TRUE){
            #    setkey(MELTED.DT,index);
            #    d3[k>=(burnin+1),loglik:=MELTED.DT[,loglik.gc]][k==(burnin+1),
            #        postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
            #        postvar.loglik.1:=0][k==(burnin+1),
            #        postmean.loglik.1:=loglik];

            #    d3[k>(burnin+1),
            #        postmean.loglik:=postmean.loglik.1 + (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
            #        postmean.exp.loglik:=postmean.exp.loglik.1 + 
            #            (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
            #        postvar.loglik:=postvar.loglik.1 + 
            #            (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
            #        postvar.true.loglik:=postvar.loglik/(k-(burnin+1))][k>(burnin+1),
            #        postmean.loglik.1:=postmean.loglik][k>(burnin+1),
            #        postvar.loglik.1:=postvar.loglik][k>(burnin+1),
            #        postmean.exp.loglik.1:=postmean.exp.loglik];
            #                       };
								   
            if(return.likmeans==TRUE){
                setkey(MELTED.DT,index);
                d3[k>=(burnin+1),loglik:=MELTED.DT[,loglik.gc]][k==(burnin+1),
				    postmean.loglik.sqrd.1:=loglik^2][k==(burnin+1),        #****NEW 10 DEC 2021****#
                    postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
                    postvar.loglik.1:=0][k==(burnin+1),
                    postmean.loglik.1:=loglik];

                d3[k>(burnin+1),
                    postmean.loglik:=postmean.loglik.1 + 
                    (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
                    postmean.loglik.sqrd:=postmean.loglik.sqrd.1 + #****NEW 10 DEC 2021****#
                    (loglik^2 - postmean.loglik.sqrd.1)/(k-burnin)][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postmean.exp.loglik:=postmean.exp.loglik.1 + 
                    (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
                    postvar.loglik:=postvar.loglik.1 + 
                    (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
                    postvar.true.loglik:=postvar.loglik/(k-(burnin+1))][k>(burnin+1),
                    postmean.loglik.1:=postmean.loglik][k>(burnin+1),
                    postmean.loglik.sqrd.1:=postmean.loglik.sqrd][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postvar.loglik.1:=postvar.loglik][k>(burnin+1),
                    postmean.exp.loglik.1:=postmean.exp.loglik];
                                     };

            sum.loglik=MELTED.DT[,sum(loglik.gc, na.rm=T)];
            sum.prior=MET.gc[,sum(d.centre, na.rm=T)] + MET.gc[,sum(d.v, na.rm=T)];
            LL.all=sum.prior + sum.loglik;
            DF.all=LL.all-LF.all;
            UF.all=runif(1);
            MET.gc[,LL:=LL.all];

            setkey(MET.gc,ML,LL);
            MET.gc[ML<LL,mlcentre:=centre][ML<LL,mlv:=v][ML<LL,ML:=LL];

            if(DF.all>10){DF.all=10};
            TEMPF.all=exp(DF.all);
            DF0.all=DF.all>0;
            UFTF.all=UF.all<TEMPF.all;
            if(sum(DF0.all,UFTF.all)>0){
			    LF.all=LL.all;accept.all=1;
			                           }else{
                MET.gc[,v:=v1][,centre:=centre1];
				accept.all=0;
                                            };

            MET.gc[k==as.integer(burnin*0.2)+1,
                meanv.1:=log(v)][k==as.integer(burnin*0.2)+1,
                meancentre.1:=centre][k==as.integer(burnin*0.2)+1,
                var11.1:=0][k==as.integer(burnin*0.2)+1,
                var12.1:=0][k==as.integer(burnin*0.2)+1,
                var22.1:=0][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv:=meanv.1 + (log(v) - meanv.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre:=meancentre.1 + (centre - meancentre.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11:=var11.1 + (log(v) - meanv.1)*(log(v) - meanv)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.v:=var11/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22:=var22.1 + (centre - meancentre.1)*(centre - meancentre)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.centre:=var22/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12:=var12.1 + ((k-as.integer(burnin*0.2+1))/(k-as.integer(burnin*0.2)))*(log(v) - meanv.1)*
                    (centre - meancentre.1)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                cov.vcentre:=var12/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv.1:=meanv][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre.1:=meancentre][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11.1:=var11][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22.1:=var22][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12.1:=var12][k>as.integer(burnin*0.799) & k<burnin+1,
                indM11:=var.v*m][k>as.integer(burnin*0.799) & k<burnin+1,
                indM12:=cov.vcentre*m][k>as.integer(burnin*0.799) & k<burnin+1,
                indM22:=var.centre*m][k==(burnin+1),
                post.mean.logv.1:=log(v)][k==(burnin+1),
                post.mean.centre.1:=centre][k==(burnin+1),
                post.var.logv.1:=0][k==(burnin+1),
                post.var.centre.1:=0][k==(burnin+1),
                var12.1:=0][k>(burnin+1),
                post.mean.logv:=post.mean.logv.1 + (log(v) - post.mean.logv.1)/(k-burnin)][k>(burnin+1),
                post.mean.centre:=post.mean.centre.1 + (centre - post.mean.centre.1)/(k-burnin)][k>(burnin+1),
                post.var.logv:=post.var.logv.1 + (log(v) - post.mean.logv.1)*(log(v) - post.mean.logv)][k>(burnin+1),
                post.var.true.logv:=post.var.logv/(k-(burnin+1))][k>(burnin+1),
                post.var.centre:=post.var.centre.1 + (centre - post.mean.centre.1)*(centre - post.mean.centre)][k>(burnin+1),
                post.var.true.centre:=post.var.centre/(k-(burnin+1))][k>(burnin+1),  #Changes 4 Nov2021 start here#
                var12:=var12.1 + ((k-(burnin+1))/(k-(burnin)))*(log(v) - post.mean.logv.1)*
                    (centre - post.mean.centre.1)][k>(burnin+1),
                cov.vcentre:=var12/(k-(burnin+1))][k>(burnin+1),
                post.mean.logv.1:=post.mean.logv][k>(burnin+1),
                post.mean.centre.1:=post.mean.centre][k>(burnin+1),
                post.var.logv.1:=post.var.logv][k>(burnin+1),
                post.var.centre.1:=post.var.centre][k>(burnin+1),
                var12.1:=var12][k==burnin+2,
				meanv.1:=NULL][k==burnin+2,
				meancentre.1:=NULL][k==burnin+2,
				var11.1:=NULL][k==burnin+2,
				var22.1:=NULL][k==burnin+2,
				meanv:=NULL][k==burnin+2,
				meancentre:=NULL][k==burnin+2,
				var11:=NULL][k==burnin+2,
				var.v:=NULL][k==burnin+2,
				var22:=NULL][k==burnin+2,
				var.centre:=NULL][k>1 & k<burnin+1,
                Zk:=r*Zk.1 + accept][k>1 & k<burnin+1,
                Wk:=r*Wk.1 + 1][k<burnin+1,
                A:=Zk/Wk][k>1 & k<burnin+1,
                Zk.1:=Zk][k>1 & k<burnin+1,
                Wk.1:=Wk][k==burnin+1,
                indM11:=indM11*(q^(A-Astar))][k==burnin+1,
                indM22:=indM22*(q^(A-Astar))][k==burnin+1,
                indM12:=indM12*(q^(A-Astar))];

            if(!is.null(plot.test.subject)){
                setkeyv(MET.gc,test.subject);
                if(k==1){
                    plot(MET.gc[plot.test.subject[1],v]~k, col=plot.col[1],xlim=c(1,nitt),ylab="Parameter value", 
                        ylim=plot.ylim,pch=plot.pch.v.centre[1],type="p");
                    abline(v=c((as.integer(burnin*0.2)+1),(as.integer(burnin*0.8)),
                        burnin,nitt),col="grey",lty=2,lwd=c(1,1,2,2));
                    abline(h=c(0,1),col=c("black","grey"),lty=2);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[2]);
                              };
                        };

                if(k>1){
                    points(MET.gc[plot.test.subject[1],v]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[1]);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };
                       };
                                           };

            MET.gc[k==as.integer(burnin*0.8),
                v:=mlv][k==as.integer(burnin*0.8),
                centre:=mlcentre][,v1:=v][,centre1:=centre][k!=as.integer(burnin*0.8),
                c("v","centre"):=rtmvnormDT2(x=MET.gc,lenv=lenv,lenu=lenu)];

            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.centre,centre:=qlogis(fix.value.centre)];
                        };
                        };

    if(poolcentre + poolv==2){
        lenv=numeric(0);
        lenu=numeric(0);

        setkeyv(data.prep.object,test.subject.esth);
        setkeyv(esth.object,test.subject.esth);
        MELT.DT=data.prep.object[esth.object];

        if(!include.Source){
            setkey(MELT.DT,Source);MELT.DT=MELT.DT[Source=="TEST"];
                   };
        setnames(MELT.DT,esth.colname,"h");

        MET.gc=unique(data.prep.object[,test.subject,with=F]);
        setkeyv(MET.gc,test.subject);

        if(poolv){lenv=1}else{lenv=nrow(MET.gc)};
        if(poolcentre){lenu=1}else{lenu=nrow(MET.gc)};

        MET.gc[,LL:=-5000000][,LF:=LL][,DF:=LL-LF][,UF:=runif(.N)][,
            TEMPF:=exp(DF)][,
            m:=(2.4/sqrt(lenv+lenu))^2][,V1:=-5000000][,
            accept:=0][,r:=0.95][,
            Zk:=1][,Zk.1:=1][,Wk:=r^0][,Wk.1:=r^0][,A:=0][,
            Astar:=0.44 - 0.21*(1/5)*(lenv + lenu - 1)][Astar<0.23,
            Astar:=0.23][,q:=20][,npar:=NA][,sum_npar:=lenv + lenu];

        if(fix.subject.v[1]!=FALSE){
            MET.gc[fix.subject.v,m:=(2.4/sqrt(lenu))^2][fix.subject.v,
                Astar:=0.44 - 0.21*(1/5)*(lenu - 1)][,sum_npar:=sum_npar - 1];
                                   };
        if(fix.subject.centre[1]!=FALSE){
            MET.gc[fix.subject.centre,m:=(2.4/sqrt(lenv))^2][fix.subject.centre,
                Astar:=0.44 - 0.21*(1/5)*(lenv - 1)][,sum_npar:=sum_npar - 1];
                                        };

        LL.all=-5000000;
        LF.all=LL.all;
        DF.all=LL.all-LF.all;
        UF.all=runif(1);
        TEMPF.all=exp(DF.all);
        accept.all=0;

        if(is.null(start.v)){
            if(lenv==1){
                MET.gc[,v:=rtlnorm(1,0,1,lower=0.02,upper=50)][,v1:=v];
                       }else{
                MET.gc[,v:=rtlnorm(.N,0,1,lower=0.02,upper=50)][,v1:=v];
                            };
                            }else{
            MET.gc[,v:=start.v][,v1:=v];
                                 };

        if(is.null(start.centre)){
            if(lenu==1){
                MET.gc[,centre:=gghybrid::rtnorm(1,mean=0,sd=1,lower=-10,upper=10)][,centre1:=centre];		
                       }else{
                MET.gc[,centre:=gghybrid::rtnorm(.N,mean=0,sd=1,lower=-10,upper=10)][,centre1:=centre];
                            };
                            }else{
            MET.gc[,centre:=qlogis(start.centre)][,centre1:=centre];
                                 };

        if(is.null(init.var.v)){
            MET.gc[,indM11:=0.01];
                               }else{
            MET.gc[,indM11:=init.var.v];
                                    };

        if(is.null(init.var.centre)){
            MET.gc[,indM22:=0.02];
                               }else{
            MET.gc[,indM22:=init.var.centre];
                                    };

        if(is.null(init.cov.vcentre)){
            MET.gc[,indM12:=0];
                                }else{
            MET.gc[,indM12:=init.cov.vcentre];
                                     };

        MET.gc[,ML:=-500000][,mlv:=v][,mlcentre:=centre];

        if(fix.subject.v[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.v,v1:=v];
                                   };

        if(fix.subject.centre[1]!=FALSE){
            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.centre,centre:=qlogis(fix.value.centre)][fix.subject.centre,centre1:=centre];
                                        };

        test=c(test.subject,"loglik.gc");

        for(k in 1:nitt){
            if(k==1){
                cat(paste("Estimating genomic clines..."),fill=1); 
                flush.console();
                    };

            if(tester::is_multiple(k,print.k)){
                cat(paste("\t","\t","\t","Iteration",k),fill=1); 
                flush.console();
                                      };

            MET.gc[,d.centre:=dnorm(centre,mean=prior.logitcentre[1],sd=prior.logitcentre[2],log=TRUE)][,
                d.v:=dnorm(log(v),mean=prior.logv[1],sd=prior.logv[2],log=TRUE)];

            setkeyv(MELT.DT,test.subject);setkeyv(MET.gc,test.subject);
            MELTED.DT=MELT.DT[MET.gc];
			
			MELTED.DT[,u:=centre*v];

#Edited to only use S0.prop_1 and S1.prop_1#
            setkey(MELTED.DT,Source_allele);
            MELTED.DT[.(1),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_1 + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_1)][.(0),
                loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*(1 - S1.prop_1) + 
                    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*(1 - S0.prop_1))];
                #loglik.gc:=log((h^v/(h^v + (1 - h)^v*exp(u)))*S1.prop_0 + 
                #    (1 - (h^v/(h^v + (1 - h)^v*exp(u))))*S0.prop_0)];

            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(-Inf),loglik.gc:=-1000];
            setkey(MELTED.DT,loglik.gc);
            MELTED.DT[.(NaN),loglik.gc:=-1000];

            #if(return.likmeans==TRUE){
            #    setkey(MELTED.DT,index);
            #    d3[k>=(burnin+1),
			#	    loglik:=MELTED.DT[,loglik.gc]][k==(burnin+1),
            #        postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
            #        postvar.loglik.1:=0][k==(burnin+1),
            #        postmean.loglik.1:=loglik];

            #    d3[k>(burnin+1),
            #        postmean.loglik:=postmean.loglik.1 + (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
            #        postmean.exp.loglik:=postmean.exp.loglik.1 + (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
            #        postvar.loglik:=postvar.loglik.1 + (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
            #        postvar.true.loglik:=postvar.loglik/(k-(burnin+1))][k>(burnin+1),
            #        postmean.loglik.1:=postmean.loglik][k>(burnin+1),
            #        postvar.loglik.1:=postvar.loglik][k>(burnin+1),
            #        postmean.exp.loglik.1:=postmean.exp.loglik];
            #                         };
									 
									 
            if(return.likmeans==TRUE){
                setkey(MELTED.DT,index);
                d3[k>=(burnin+1),loglik:=MELTED.DT[,loglik.gc]][k==(burnin+1),
				    postmean.loglik.sqrd.1:=loglik^2][k==(burnin+1),        #****NEW 10 DEC 2021****#
                    postmean.exp.loglik.1:=exp(loglik)][k==(burnin+1),
                    postvar.loglik.1:=0][k==(burnin+1),
                    postmean.loglik.1:=loglik];

                d3[k>(burnin+1),
                    postmean.loglik:=postmean.loglik.1 + 
                    (loglik - postmean.loglik.1)/(k-burnin)][k>(burnin+1),
                    postmean.loglik.sqrd:=postmean.loglik.sqrd.1 + #****NEW 10 DEC 2021****#
                    (loglik^2 - postmean.loglik.sqrd.1)/(k-burnin)][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postmean.exp.loglik:=postmean.exp.loglik.1 + 
                    (exp(loglik) - postmean.exp.loglik.1)/(k-burnin)][k>(burnin+1),
                    postvar.loglik:=postvar.loglik.1 + 
                    (loglik - postmean.loglik.1)*(loglik - postmean.loglik)][k>(burnin+1),
                    postvar.true.loglik:=postvar.loglik/(k-(burnin+1))][k>(burnin+1),
                    postmean.loglik.1:=postmean.loglik][k>(burnin+1),
                    postmean.loglik.sqrd.1:=postmean.loglik.sqrd][k>(burnin+1),#****NEW 10 DEC 2021****#
                    postvar.loglik.1:=postvar.loglik][k>(burnin+1),
                    postmean.exp.loglik.1:=postmean.exp.loglik];
                                     };

            sum.loglik=MELTED.DT[,sum(loglik.gc, na.rm=T)];

            sum.prior=MET.gc[,sum(d.centre, na.rm=T)] + MET.gc[,sum(d.v, na.rm=T)];
            LL.all=sum.prior + sum.loglik;
            DF.all=LL.all-LF.all;
            UF.all=runif(1);
            MET.gc[,LL:=LL.all];

            setkey(MET.gc,ML,LL);
            MET.gc[ML<LL,mlcentre:=centre][ML<LL,mlv:=v][ML<LL,ML:=LL];

            if(DF.all>10){DF.all=10};
            TEMPF.all=exp(DF.all);
            DF0.all=DF.all>0;
            UFTF.all=UF.all<TEMPF.all;
            if(sum(DF0.all,UFTF.all)>0){LF.all=LL.all;accept.all=1}else{
                MET.gc[,v:=v1][,centre:=centre1];accept.all=0;
                                                                       };

            MET.gc[k==as.integer(burnin*0.2)+1,
                meanv.1:=log(v)][k==as.integer(burnin*0.2)+1,
                meancentre.1:=centre][k==as.integer(burnin*0.2)+1,
                var11.1:=0][k==as.integer(burnin*0.2)+1,
                var12.1:=0][k==as.integer(burnin*0.2)+1,
                var22.1:=0][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv:=meanv.1 + (log(v) - meanv.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre:=meancentre.1 + (centre - meancentre.1)/(k-as.integer(burnin*0.2))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11:=var11.1 + (log(v) - meanv.1)*(log(v) - meanv)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.v:=var11/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22:=var22.1 + (centre - meancentre.1)*(centre - meancentre)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var.centre:=var22/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12:=var12.1 + ((k-as.integer(burnin*0.2+1))/(k-as.integer(burnin*0.2)))*(log(v) - meanv.1)*
                    (centre - meancentre.1)][k>as.integer(burnin*0.2)+1&k<burnin+1,
                cov.vcentre:=var12/(k-as.integer(burnin*0.2+1))][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meanv.1:=meanv][k>as.integer(burnin*0.2)+1&k<burnin+1,
                meancentre.1:=meancentre][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var11.1:=var11][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var22.1:=var22][k>as.integer(burnin*0.2)+1&k<burnin+1,
                var12.1:=var12][k>as.integer(burnin*0.799) & k<burnin+2,
                indM11:=var.v*m][k>as.integer(burnin*0.799) & k<burnin+2,
                indM12:=cov.vcentre*m][k>as.integer(burnin*0.799) & k<burnin+2,
                indM22:=var.centre*m][k==(burnin+1),
                post.mean.logv.1:=log(v)][k==(burnin+1),
                post.mean.centre.1:=centre][k==(burnin+1),
                post.var.logv.1:=0][k==(burnin+1),
                post.var.centre.1:=0][k==(burnin+1),
                var12.1:=0][k>(burnin+1),
                post.mean.logv:=post.mean.logv.1 + (log(v) - post.mean.logv.1)/(k-burnin)][k>(burnin+1),
                post.mean.centre:=post.mean.centre.1 + (centre - post.mean.centre.1)/(k-burnin)][k>(burnin+1),
                post.var.logv:=post.var.logv.1 + (log(v) - post.mean.logv.1)*(log(v) - post.mean.logv)][k>(burnin+1),
                post.var.true.logv:=post.var.logv/(k-(burnin+1))][k>(burnin+1),
                post.var.centre:=post.var.centre.1 + (centre - post.mean.centre.1)*(centre - post.mean.centre)][k>(burnin+1),
                post.var.true.centre:=post.var.centre/(k-(burnin+1))][k>(burnin+1),  #Changes 4 Nov2021 start here#
                var12:=var12.1 + ((k-(burnin+1))/(k-(burnin)))*(log(v) - post.mean.logv.1)*
                    (centre - post.mean.centre.1)][k>(burnin+1),
                cov.vcentre:=var12/(k-(burnin+1))][k>(burnin+1),
                post.mean.logv.1:=post.mean.logv][k>(burnin+1),
                post.mean.centre.1:=post.mean.centre][k>(burnin+1),
                post.var.logv.1:=post.var.logv][k>(burnin+1),
                post.var.centre.1:=post.var.centre][k>(burnin+1),
                var12.1:=var12][k==burnin+2,
				meanv.1:=NULL][k==burnin+2,
				meancentre.1:=NULL][k==burnin+2,
				var11.1:=NULL][k==burnin+2,
				var22.1:=NULL][k==burnin+2,
				meanv:=NULL][k==burnin+2,
				meancentre:=NULL][k==burnin+2,
				var11:=NULL][k==burnin+2,
				var.v:=NULL][k==burnin+2,
				var22:=NULL][k==burnin+2,
				var.centre:=NULL];

            if(!is.null(plot.test.subject)){
                setkeyv(MET.gc,test.subject);

                if(k==1){
                    plot(MET.gc[plot.test.subject[1],v]~k, col=plot.col[1],xlim=c(1,nitt),ylab="Parameter value", 
                        ylim=plot.ylim,pch=plot.pch.v.centre[1],type="p");
                    abline(v=c((as.integer(burnin*0.2)+1),(as.integer(burnin*0.8)),
                        burnin,nitt),col="grey",lty=2,lwd=c(1,1,2,2));
                    abline(h=c(0,1),col=c("black","grey"),lty=2);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],pch=plot.pch.v.centre[2]);
                                                   };
                        };

                if(k>1){
                    points(MET.gc[plot.test.subject[1],v]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[1]);
                    points(MET.gc[plot.test.subject[1],centre]~k,col=plot.col[1],cex=0.5,
                        pch=plot.pch.v.centre[2]);

                    if(length(plot.test.subject)>1){
                        points(MET.gc[plot.test.subject[2],v]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[2],centre]~k,
                            col=plot.col[2],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };

                    if(length(plot.test.subject)>2){
                        points(MET.gc[plot.test.subject[3],v]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[1]);
                        points(MET.gc[plot.test.subject[3],centre]~k,
                            col=plot.col[3],cex=0.5,pch=plot.pch.v.centre[2]);
                                                   };
                       };
                                           };

            MET.gc[k==as.integer(burnin*0.8),
                v:=mlv][k==as.integer(burnin*0.8),
                centre:=mlcentre][,v1:=v][,centre1:=centre][k!=as.integer(burnin*0.8),
                c("v","centre"):=rtmvnormDT2(x=MET.gc,lenv=lenv,lenu=lenu)];

            setkeyv(MET.gc,test.subject);
            MET.gc[fix.subject.v,v:=fix.value.v][fix.subject.centre,centre:=qlogis(fix.value.centre)];

                        };
                        };


#Editing the output 04 Nov 2021#
    MET.gc[,logv.ci.lower:=qnorm(0.025,mean=post.mean.logv,sd=sqrt(post.var.true.logv))][,
        logv.ci.upper:=qnorm(0.975,mean=post.mean.logv,sd=sqrt(post.var.true.logv))];

    setkey(MET.gc,post.mean.logv);
    MET.gc[post.mean.logv<=0,
	    pval.v:=(pnorm(0,mean=post.mean.logv,sd=sqrt(post.var.true.logv),lower.tail=FALSE))*2][post.mean.logv>0,
        pval.v:=(pnorm(0,mean=post.mean.logv,sd=sqrt(post.var.true.logv),lower.tail=TRUE))*2];

    MET.gc[,
	    logit.centre.ci.lower:=qnorm(0.025,mean=post.mean.centre,sd=sqrt(post.var.true.centre))][,
        logit.centre.ci.upper:=qnorm(0.975,mean=post.mean.centre,sd=sqrt(post.var.true.centre))];

    setkey(MET.gc,post.mean.centre);
    MET.gc[post.mean.centre<=0,
        pval.centre:=(pnorm(0,mean=post.mean.centre,sd=sqrt(post.var.true.centre),lower.tail=FALSE))*2][post.mean.centre>0,
        pval.centre:=(pnorm(0,mean=post.mean.centre,sd=sqrt(post.var.true.centre),lower.tail=TRUE))*2];

    MET.gc[,
	    post.mean.v:=exp(post.mean.logv)][,
        v.ci.lower:=exp(logv.ci.lower)][,
        v.ci.upper:=exp(logv.ci.upper)];
		
	setnames(MET.gc,c("post.mean.centre"),c("post.mean.logit.centre"));
		
	MET.gc[,
	    post.mean.centre:=plogis(post.mean.logit.centre)][,
		centre.ci.lower:=plogis(logit.centre.ci.lower)][,
		centre.ci.upper:=plogis(logit.centre.ci.upper)];

    setnames(MET.gc,c("post.mean.v","post.mean.centre","v.ci.lower","v.ci.upper","pval.v","centre.ci.lower","centre.ci.upper","pval.centre",
	    "post.mean.logv","post.var.true.logv","post.mean.logit.centre","post.var.true.centre","cov.vcentre"),
        c("exp_mean_log_v","invlogit_mean_logit_centre","v_lower_95","v_upper_95","v_pvalue","centre_lower_95","centre_upper_95","centre_pvalue",
		"mean_log_v","var_log_v","mean_logit_centre","var_logit_centre","cov_log_v_logit_centre"));
    
	if(test.subject=="locus"){
	    laf=unique(data.prep.object[,.(locus,S0.prop_1,S1.prop_1)]);
		setkey(laf,locus);setkey(MET.gc,locus);
		MET.gc=laf[MET.gc];	
        MET.gc=MET.gc[,c(test.subject,"S0.prop_1","S1.prop_1","exp_mean_log_v","v_lower_95","v_upper_95","v_pvalue","invlogit_mean_logit_centre","centre_lower_95","centre_upper_95","centre_pvalue",
		"mean_log_v","var_log_v","mean_logit_centre","var_logit_centre","cov_log_v_logit_centre","npar","sum_npar"),with=F];	
	                         }else{
        MET.gc=MET.gc[,c(test.subject,"exp_mean_log_v","v_lower_95","v_upper_95","v_pvalue","invlogit_mean_logit_centre","centre_lower_95","centre_upper_95","centre_pvalue",
		"mean_log_v","var_log_v","mean_logit_centre","var_logit_centre","cov_log_v_logit_centre","npar","sum_npar"),with=F];
		                          };

    MET.gc[,
        exp_mean_log_v:=round(exp_mean_log_v,digits=3)][,
        v_lower_95:=round(v_lower_95,digits=3)][,
        v_upper_95:=round(v_upper_95,digits=3)][,
        invlogit_mean_logit_centre:=round(invlogit_mean_logit_centre,digits=3)][,
        centre_lower_95:=round(centre_lower_95,digits=3)][,
        centre_upper_95:=round(centre_upper_95,digits=3)];

    setkeyv(MET.gc,test.subject);
    MET.gc[fix.subject.v,
	    c("v_lower_95","v_upper_95","v_pvalue"):=NA][fix.subject.centre,
        c("centre_lower_95","centre_upper_95","centre_pvalue"):=NA];
    
	if(MET.gc[v_pvalue==0|centre_pvalue==0,.N]>0){
	    warning("Some p values for v and/or centre are estimated as exactly 0; please be aware that these may cause problems with downstream calculations such as false discovery rate.");
		                                         };

    output=list();
	    output$init.var.v=init.var.v;
		output$init.var.centre=init.var.centre;
		output$init.cov.vcentre=init.cov.vcentre;
        output$gc=MET.gc;
        output$test.subject=test.subject;		
		
    if(return.likmeans==TRUE){
        d3[,loglik:=NULL][,postmean.exp.loglik.1:=NULL][,postmean.loglik.sqrd.1:=NULL][,#****NEW 10 DEC 2021***#
            postmean.loglik.1:=NULL][,postvar.loglik:=NULL][,postvar.loglik.1:=NULL];
        setnames(d3,"postvar.true.loglik","postvar.loglik");
        output$likmeans=d3;
                             };

    cat(paste("Done"),fill=1);

    return(output)
               }
