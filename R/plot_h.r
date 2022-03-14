#' Plot function for estimated hybrid indices and credible intervals.
#'
#' @param data Name of an \code{hi} object produced by the \code{esth} function.
#' @param test.subject \code{test.subject} object from \code{esth}. Default \dQuote{INDLABEL}. 
#' @param POPID.name Character string. Name of the designated \code{POPID} field. Default \dQuote{POPID}, but this is ignored if there is no \code{POPID} field.
#' @param mean.h.by Character string. Calculate mean hybrid index by group for the group designations in this field. Creates a column called \sQuote{mean_h} 
#'   for ordering data points on the figure using \option{sort.by}.
#' @param sort.by Character vector. Fields by which to sort the test subjects by value along the x axis in the figure. Can include \sQuote{mean_h}. 
#'   Default is \sQuote{h_posterior_mode}.
#' @param col.group Character string. Variable by which to group colouring of data points.
#' @param group.sep Character string. Variable by which to separate groups of data points with vertical grey lines in the plot.
#' @param fill.source Logical. Whether to fill the area of the plot representing parental reference samples with a grey background. 
#'   Default \code{FALSE}.
#' @param basic.lines Logical. If \code{TRUE}, adds horizontal grey dashed lines at hybrid index 0, 0.5 and 1. Default \code{TRUE}.
#' @param source.col Character vector. Colours to be used for the two parental reference set data points.
#' @param source.limits Character vector. Colours of horizontal dashed lines marking the innermost credible interval for each 
#'   parental reference set. Used to indicate whether test subjects are clearly outside the hybrid index ranges of parentals.
#' @param custom.abline An \code{abline} of the user's choosing.
#' @param ... Further graphical parameters.
#' @return Alongside the plot, returns a \code{data.table} and \code{data.frame} useful for creating a legend, with 
#'   the specified \code{test.subject} field, the \code{mean.h.by} field if specified, along with the \code{mean_h} value,
#'   the \code{Source} field, the colours in a field named \code{col.Dark2}, and the row number for ordering in field \code{rn}.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' abc = plot_h(data=hindlabel$hi,
#' test.subject=hindlabel$test.subject,
#' mean.h.by="POPID",			                    #Calculate the mean hybrid index for each value of the "POPID" column#
#' sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
#'							                         #individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
#' col.group="POPID",
#' group.sep="POPID",
#' fill.source=TRUE,
#' basic.lines=FALSE,
#' source.col=c("blue","red"),
#' source.limits=c("blue","red"),
#' custom.abline=abline(h=0.75,col="magenta",lwd=2),
#' cex=1,pch=16,
#' cex.lab=1.5,cex.main=1.5,ylim=c(0,1))
#'
#' #Adding a legend with base R's "legend" function, using the plot_h object created above, "abc":
#' setkey(abc,rn);		#Order data by row number#
#' legend("topleft",	#Place the legend in the top left of the figure#
#' abc[,POPID], 		#Name of the field by which data point colours are grouped#
#' bg="white",			#Background colour#
#' text.col=c("black"), #Text colour#
#' pch=22, 				#Text size#
#' col=abc[,col.Dark2], #Name of the field containing colour information#
#' pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
#' ncol=4,				#Number of columns for splitting the group names#
#' cex=0.6, pt.cex=0.7)
#' }
#' @export
plot_h = function(data,test.subject="INDLABEL",POPID.name="POPID",
    mean.h.by=NULL,sort.by="h_posterior_mode",
    col.group=NULL,group.sep=NULL,fill.source=FALSE,
    basic.lines=TRUE,source.col=NULL,source.limits=NULL,
    custom.abline=NULL,...){

    setkey(data,NULL);

    if(is.null(col.group)==FALSE){
        rf = colorRampPalette(rev(RColorBrewer::brewer.pal(8,'Dark2')));
        col.Dark2 = rf(nrow(unique(data[,col.group,with=F])));
        col.Dark2 = data.table(col.Dark2);

        setkeyv(data,col.group);
        pop = data.table(unique(data[,col.group,with=F]));
        pop.col = cbind(pop,col.Dark2);

        setkeyv(pop.col,col.group);
        data=data[pop.col];

        setkey(data,NULL);setkey(data,col.Dark2);
        col.summary = data[,head(.SD,1),by=col.Dark2]
                                 };

    setkey(data,NULL);

    if(is.na(chmatch("Source",sort.by))==FALSE){
        setkey(data,Source);
        data["S0",Source2:="A"];
        data["TEST",Source2:="B"];
        data["S1",Source2:="C"];
        sort.by[chmatch("Source",sort.by)]="Source2";
                                               };

    setkey(data,NULL);
    setkeyv(data,test.subject);

    if(is.null(mean.h.by)){
        setkeyv(data,sort.by);
        data[,rn:=row(data)[,1]];
		                  }else{
        setkeyv(data,mean.h.by);
        data[,mean_h:=mean(h_posterior_mode),by=mean.h.by];
        setkeyv(data,sort.by);
        data[,rn:=row(data)[,1]];
                               };

    setkey(data,NULL);

    plot(data[,h_posterior_mode]~data[,rn],type="n",
        #main="Hybrid index and 95% C.I.",               #Probably better without this - user can add a title if they want#
        xlab="Individual",ylab="Hybrid index",...);

    if(is.null(group.sep)==FALSE){
        setkeyv(data,group.sep);
        abline(v=0.5,col="grey");
        abline(v=data[,(max(rn)+0.5),by=group.sep]$V1,col="grey");
                                 };

    setkey(data,NULL);
			
	setnames(data,POPID.name,"POPID",skip_absent=TRUE); #skip_absent=TRUE added 4 March 2022#

    if(fill.source==TRUE){
        setkey(data,Source);

    for(i in 1:(data[Source=="S0",length(unique(POPID))])){
        rect(data[Source=="S0",(min(rn)),by=POPID]$V1[i] - 0.5, 
            par("usr")[3], 
            data[Source=="S0",(max(rn)),by=POPID]$V1[i] + 0.5, 
            par("usr")[4], col = "grey95");
                                                          };

    for(i in 1:(data[Source=="S1",length(unique(POPID))])){
        rect(data[Source=="S1",(min(rn)),by=POPID]$V1[i] - 0.5, 
            par("usr")[3], 
            data[Source=="S1",(max(rn)),by=POPID]$V1[i] + 0.5, 
            par("usr")[4], col = "grey95");
                                                          };
                         };

    setkey(data,NULL);

    if(basic.lines==TRUE){
        abline(h=c(0,0.5,1),col="grey",lty=2,lwd=2);
                         };

    if(is.null(source.col)==FALSE){
        setkey(data,Source);

        if(is.null(col.group)){
            data[,col.Dark2:="black"];
                              };
        data["S0",col.Dark2:=source.col[1]];
        data["S1",col.Dark2:=source.col[2]];

        setkey(data,NULL);setkey(data,col.Dark2);
        col.summary = data[,head(.SD,1),by=col.Dark2];
        setkey(col.summary,Source);
        col.summary["S0",POPID:="S0"];
        col.summary["S1",POPID:="S1"];
                                  };

    setkey(data,NULL);

    if(is.null(col.group) & is.null(source.col)){
        points(data[,h_posterior_mode]~data[,rn],...);
        arrows(data[,
            rn],data[,h_cred_int_lower],data[,
            rn],data[,h_cred_int_upper],angle=90,code=3,length=0);
                                                }else{
        points(data[,h_posterior_mode]~data[,rn],
            col=data[,col.Dark2],...);
        arrows(data[,
            rn],data[,h_cred_int_lower],data[,
            rn],data[,h_cred_int_upper],angle=90,code=3,length=0,
            col=data[,col.Dark2]);
                                                     };

    if(is.null(source.limits)==FALSE){
        abline(h=data[Source=="S0",max(h_cred_int_upper)],
            col=source.limits[1],lty=2,lwd=2);
        abline(h=data[Source=="S1",min(h_cred_int_lower)],
            col=source.limits[2],lty=2,lwd=2);
                                     };

    if(is.null(custom.abline)==FALSE){
        custom.abline;
                                     };
									 
	setnames(col.summary,"POPID",POPID.name);
									 
    if(is.null(mean.h.by)){
        col.summary=col.summary[,c(test.subject,"Source","col.Dark2","rn"),with=F];
                          }else{
        col.summary=col.summary[,c(test.subject,mean.h.by,"Source","mean_h","col.Dark2","rn"),with=F];
                               };
							   
	return(col.summary)
}