#' Plot the fitted cline for one or more test subjects using base R graphics.
#'
#' @param ggcline.object A \code{gc} object from \code{ggcline}.
#' @param cline.locus Character vector. Names of the loci or \code{ggcline} test subjects for which to plot cline curves.
#' @param locus.column Character string. Name of the field containing the test subjects. Default \dQuote{locus}.
#' @param cline.col Character vector. Colours for the plotted cline curves. Can be a single value to be applied to all curves 
#'   or a vector the same length as \code{cline.locus}. Default \dQuote{black}.
#' @param null.line.locus Character vector. Names of loci for which to plot a black dashed line showing the null expectation. Default \code{NULL}.
#' @param null.line.col Character vector. Colours for the plotted null dashed lines. Default \dQuote{black} (repeated to the number of loci plotted).
#' @param cline.centre.line Character vector. Names of \code{ggcline} test subjects (loci) for which to plot a vertical dashed line 
#'   indicating the cline centre. Default \code{NULL}.
#' @param cline.centre.col Character vector. Colours for the plotted centre lines.
#' @param ... Further graphical parameters.
#' @details The function plots the cline curve(s), and genotype data (scaled from 0 to 1) can then be added as points. Samples from the parameter posteriors 
#'   can be taken using the \code{gghybrid} function \code{rtmvnormDT3}, and their resulting clines added to the plot to indicate uncertainty, 
#'   as described in the examples.
#' @return a plot. No object is returned.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' ########################################################################
#' #Plot a genomic cline for one locus and add uncertainty and data points#
#' ########################################################################
#' 
#' #Plot the curve for locus CLTA, using a previously created ggcline object (see examples for ggcline).
#' plot_clinecurve(
#' ggcline.object=gc$gc,
#' cline.locus="CLTA",
#' locus.column="locus",
#' cline.col="#E495A5",
#' cline.centre.line="CLTA",
#' cline.centre.col="black"
#' )
#' 
#' #Add a title and axis labels.
#' title(main = "CLTA",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)
#' 
#' ###################################################
#' #Add sampled cline curves to represent uncertainty#
#' ###################################################
#' 
#' #There is no option for shading credible intervals within plot_clinecurve, so here is a suggestion for how to add some: 
#' #sample parameter values from the joint posterior, and add a grey curve to the plot for each posterior sample. 
#' #This will overlay the original cline curve plotted above, so it needs to be plotted again.
#' 
#' #First generate random samples from the posterior using gghybrid's rtmvnormDT3 function.
#' 
#' gcsamp=gc$gc[locus%in%c("CLTA"),
#'  rtmvnormDT3(
#'   nsamp=1000,                                #default is 1 sample per locus#
#'   meanlogv=mean_log_v,                       #no default, column must be specified#
#'   meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
#'   varlogv=var_log_v,                         #no default, column must be specified#
#'   varlogitcentre=var_logit_centre,           #no default, column must be specified#
#'   covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
#'   lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
#'   upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
#'             ),                               #end of rtmvnormDT3 function#
#'  by="locus"];                                #indicate a grouping variable for the sampling#
#' 
#' #Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
#' locipost=gc$gc[locus=="CLTA",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
#' gcsamp=rbind(gcsamp,locipost)
#' 
#' #Add the individual curves to the plot in a 'for' loop.
#' for(i in 1:nrow(gcsamp)){
#'         v = gcsamp[i,exp(log_v)]
#'         u = gcsamp[i,logit_centre*exp(log_v)]
#'         S1.prop_1 = gc$gc[locus=="CLTA",S1.prop_1];
#'         S0.prop_1 = gc$gc[locus=="CLTA",S0.prop_1];
#' 
#'         par(new=T);
#' 
#' if(i < nrow(gcsamp)){
#'         curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
#' 	        from=0,to=1,axes=F,xlab="",ylab="", col="grey",lwd=0.5,ylim=c(0,1));
#' }else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
#'         curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
#' 	        from=0,to=1,axes=F,xlab="",ylab="", col="#E495A5",lwd=3,ylim=c(0,1));
#' };
#'                          };#end of for loop#
#' 
#' ######################################################
#' #Add data to the plot as genotypes scaled from 0 to 1#
#' ######################################################
#' 
#' ploidy=2
#' 
#' #Create a data.table containing individual reference, hybrid index, locus name and genotype. See the ggcline examples for creating prepdata and hindlabel.
#' setkey(prepdata$data.prep,INDLABEL);setkey(hindlabel$hi,INDLABEL);
#' genodat=hindlabel$hi[,.(INDLABEL,h_posterior_mode)][prepdata$data.prep[locus=="CLTA",sum(Source_allele)/ploidy,by=c("INDLABEL","locus")]]
#' 
#' #Add the resulting points to the cline plot.
#' points(genodat[,V1]~genodat[,h_posterior_mode],pch=16,cex=0.7,col="#E495A5")
#' 
#' 
#' ######################################################
#' #Plot cline curves for multiple loci on the same plot#
#' ######################################################
#' 
#' #No data points or uncertainty included this time.
#' 
#' plot_clinecurve(
#' ggcline.object=gc$gc,
#' cline.locus=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
#' locus.column="locus",
#' cline.col=c("orange","blue","green","red","magenta"),
#' null.line.locus=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
#' null.line.col=c("orange","blue","green","red","magenta"),
#' cline.centre.line=c("ACO1","EGR1","A2M","HSDL2","RPS4"),
#' cline.centre.col=c("orange","blue","green","red","magenta")
#' )
#' 
#' #Add a title and axis labels:
#' title(main = "ACO1,EGR1,A2M,HSDL2,RPS4",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)
#' }
#' @export
plot_clinecurve = function(ggcline.object,
    cline.locus,locus.column="locus",cline.col="black",
    null.line.locus=NULL,null.line.col="black",
    cline.centre.line=NULL,cline.centre.col="black",...){

    setkeyv(ggcline.object,locus.column);

    v = numeric(1);
    u = numeric(1);
    S1.prop_1 = numeric(1);
    S0.prop_1 = numeric(1);

    #plot.new();

    plot(0.5,0.5,type="n",xlim=c(0,1),ylim=c(0,1),
        ann=FALSE,...);

    if(is.null(null.line.locus)==FALSE){
        setkeyv(ggcline.object,locus.column);
	    if(length(null.line.col)==1){null.line.col=rep(null.line.col,length(null.line.locus))};
        for(i in 1:length(null.line.locus)){
            S1.prop_1 = unique(ggcline.object[null.line.locus[i],S1.prop_1]);#******************PRESENT IN $gc*********************#
            S0.prop_1 = unique(ggcline.object[null.line.locus[i],S0.prop_1]);#******************PRESENT IN $gc*********************#
            abline(a=S0.prop_1,b=(S1.prop_1 - S0.prop_1),lty=2,lwd=2,col=null.line.col[i]);
                                           };
                                       };

    if(is.null(cline.centre.line)==FALSE){
	    setkeyv(ggcline.object,locus.column);
	    if(length(cline.centre.col)==1){cline.centre.col=rep(cline.centre.col,length(cline.centre.line))};
	        for(i in 1:length(cline.centre.line)){
                centre=ggcline.object[cline.centre.line[i],invlogit_mean_logit_centre];
                abline(v=centre,col=cline.centre.col[i],lty=2,lwd=2);
                #setkey(ggcline.object,NULL);
						                         };
									     };

    setkeyv(ggcline.object,locus.column);
	if(length(cline.col)==1){cline.col=rep(cline.col,length(cline.locus))};
    for(i in 1:length(cline.locus)){
        v = ggcline.object[cline.locus[i],exp_mean_log_v];
        u = qlogis(ggcline.object[cline.locus[i],invlogit_mean_logit_centre])*ggcline.object[cline.locus[i],exp_mean_log_v];
        S1.prop_1 = unique(ggcline.object[cline.locus[i],S1.prop_1]);#******************PRESENT IN $gc*********************#
        S0.prop_1 = unique(ggcline.object[cline.locus[i],S0.prop_1]);#******************PRESENT IN $gc*********************#

        par(new=T);

        curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), 
	        from=0,to=1,axes=F,xlab="",ylab="", col=cline.col[i],lwd=3,ylim=c(0,1));

        if(i == length(cline.locus)){par(new=F)};
                                  };

                                                        }
