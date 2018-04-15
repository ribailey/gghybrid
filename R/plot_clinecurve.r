#' Plot the fitted cline for one or more test subjects.
#'
#' @param ggcline.object A \code{gc} object from \code{ggcline}.
#' @param data.prep.object A \code{data.prep} object from \code{data.prep}.
#' @param esth.object An \code{hi} object from \code{esth}.
#' @param test.subject Character string. The test subject used in \code{esth}.
#' @param include.Source Logical. If plotting individual data points, whether to include parental reference data. Default \code{FALSE}. 
#'   Set to \code{TRUE} if parental reference sets were included in \code{ggcline}.
#' @param cline.locus Character vector. Names of the loci or \code{ggcline} test subjects for which to plot cline curves.
#' @param locus.column Character string. Name of the field containing the test subjects. Default \dQuote{locus}.
#' @param cline.col Character vector. Colours for the plotted cline curves. Default \dQuote{black} (repeated to the number of loci plotted).
#' @param null.line.locus Character vector. Names of loci for which to plot a black dashed line showing the null expectation. Default \code{NULL}.
#' @param null.line.col Character vector. Colours for the plotted null dashed lines. Default \dQuote{black} (repeated to the number of loci plotted).
#' @param plot.data Logical default, otherwise a character vector. Names of loci for which to plot data points on the figure.
#' @param data.col Character vector. Colours for the plotted data (one for each locus/\code{ggcline} test subject).
#' @param plot.genotype Logical. Whether to plot genotypes of individuals (or \code{esth} test subjects) scaled from 0 to 1 rather than individual 
#'   allele copies (which will all take values 0 or 1).
#' @param PLOIDY Numeric. The number of allele copies for each \code{esth} individual test subject (typically \code{INDLABEL}).
#' @param cline.centre.line Character vector. Names of \code{ggcline} test subjects (loci) for which to plot a vertical dashed line 
#'   indicating the cline centre. Default \code{NULL}.
#' @param cline.centre.col Character vector. Colours for the plotted centre lines.
#' @param ... Further graphical parameters.
#' @return a plot.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' plot_clinecurve(
#' ggcline.object=gc1$gc,
#' data.prep.object=prepdata$data.prep,
#' esth.object=hindlabel$hi,
#' cline.locus="locusgroup1",		#Plotting the cline for only one locus#
#' locus.column="locusgroup",
#' #cline.col=c("orange","blue","green","red","magenta"),		#If you had chosen to plot 5 loci#
#' cline.col="green",
#' #null.line.locus=c("locusgroup1","locusgroup2","locusgroup3","locusgroup4","locusgroup5"),	#If you had chosen to plot 5 loci#
#' null.line.locus="locusgroup1",
#' #null.line.col=c("orange","blue","green","red","magenta"),	#If you had chosen to plot 5 loci#
#' null.line.col="green",
#' plot.data="locusgroup1",
#' data.col="black",
#' plot.genotype=TRUE,
#' PLOIDY=2,
#' cline.centre.line="locusgroup1",
#' cline.centre.col="green"
#' #cline.centre.line=c("locusgroup1","locusgroup2","locusgroup3","locusgroup4","locusgroup5"),	#If you had chosen to plot 5 loci#
#' #cline.centre.col=c("orange","blue","green","red","magenta")	#If you had chosen to plot 5 loci#
#' )
#' #Add a title:
#' title(main = "locusgroup1",xlab="Hybrid index",
#' ylab="Genotype frequency",cex.main=1.5,cex.lab=1.5)
#' }
#' @export
plot_clinecurve = function(ggcline.object,data.prep.object,esth.object,
    test.subject = "INDLABEL",
    include.Source = FALSE,
    cline.locus,locus.column="locus",cline.col="black",
    null.line.locus=NULL,null.line.col="black",
    plot.data=NULL,data.col=NULL,
    plot.genotype=FALSE,PLOIDY,
    cline.centre.line=NULL,cline.centre.col="black",...){

    if(!include.Source){
        setkey(data.prep.object,Source);
        data.prep.object = data.prep.object[Source == "TEST"];
        esth.object = esth.object[Source == "TEST"];
                       };

    setkeyv(data.prep.object,test.subject);setkeyv(esth.object,test.subject);
    prep.hi = data.prep.object[esth.object];

    setkeyv(prep.hi,locus.column);
    setkeyv(ggcline.object,locus.column);

    v = numeric(1);
    u = numeric(1);
    S1.prop_1 = numeric(1);
    S0.prop_1 = numeric(1);

    plot(prep.hi[cline.locus[1],
        Source_allele]~prep.hi[cline.locus[1],
        qlogis(h_posterior_mode)],type="n",xlim=c(0,1),
        ann=FALSE,...);

    if(is.null(null.line.locus)==FALSE){
		if(null.line.col=="black"){null.line.col=rep("black",length(null.line.locus))};
        for(i in 1:length(null.line.locus)){
            S1.prop_1 = unique(prep.hi[null.line.locus[i],S1.prop_1]);
            S0.prop_1 = unique(prep.hi[null.line.locus[i],S0.prop_1]);
            abline(a=S0.prop_1,b=(S1.prop_1 - S0.prop_1),lty=2,lwd=2,col=null.line.col[i]);
                                           };
                                       };

    if(is.null(cline.centre.line)==FALSE){
        centre=ggcline.object[cline.centre.line,centre_mean];
        abline(v=centre,col=cline.centre.col,lty=2,lwd=2);
        setkey(ggcline.object,NULL);
						                 };

    setkeyv(ggcline.object,locus.column);
	if(cline.col=="black"){cline.col=rep("black",length(cline.locus))};
    for(i in 1:length(cline.locus)){
        v = ggcline.object[cline.locus[i],v_mean];
        u = qlogis(ggcline.object[cline.locus[i],centre_mean])*ggcline.object[cline.locus[i],v_mean];
        S1.prop_1 = unique(prep.hi[cline.locus[i],S1.prop_1]);
        S0.prop_1 = unique(prep.hi[cline.locus[i],S0.prop_1]);

        par(new=T);

        curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), 
	        from=0,to=1,axes=F,xlab="",ylab="", col=cline.col[i],lwd=3,ylim=c(0,1));

        if(i == length(cline.locus)){par(new=F)};
                                  };

    if(is.null(plot.data)==FALSE){

        setnames(prep.hi,locus.column,"loc");

        if(plot.genotype==TRUE){
            setnames(prep.hi,test.subject,"INDLABEL");
            prep.hi.geno= melt(dcast(prep.hi, INDLABEL + h_posterior_mode~loc,
                value.var="Source_allele",fun=sum), 
                id = c(test.subject,"h_posterior_mode"),variable.name = locus.column, 
                value.name = "genotype");
            prep.hi.geno[,genotype:=genotype/PLOIDY];

            setkey(prep.hi.geno,loc);

            for(i in 1:length(plot.data)){
                points(prep.hi.geno[cline.locus[i],
                    genotype]~prep.hi.geno[cline.locus[i],
                    h_posterior_mode],col=data.col[i],...);
                                         };
										 
            setnames(prep.hi,"INDLABEL",test.subject);									 
                               }else{
            setkey(prep.hi,loc);
            for(i in 1:length(plot.data)){
                points(prep.hi[cline.locus[i],
                    Source_allele]~prep.hi[cline.locus[i],
                    h_posterior_mode],col=data.col[i],...);
                                         };
                                    };
                                 };
                                                        }