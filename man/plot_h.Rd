% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_h.r
\name{plot_h}
\alias{plot_h}
\title{Plot function for estimated hybrid indices and credible intervals.}
\usage{
plot_h(
  data,
  test.subject = "INDLABEL",
  POPID.name = "POPID",
  mean.h.by = NULL,
  sort.by = "h_posterior_mode",
  col.group = NULL,
  group.sep = NULL,
  fill.source = FALSE,
  basic.lines = TRUE,
  source.col = NULL,
  source.limits = NULL,
  custom.abline = NULL,
  ...
)
}
\arguments{
\item{data}{Name of an \code{hi} object produced by the \code{esth} function.}

\item{test.subject}{\code{test.subject} object from \code{esth}. Default \dQuote{INDLABEL}.}

\item{POPID.name}{Character string. Name of the designated \code{POPID} field. Default \dQuote{POPID}, but this is ignored if there is no \code{POPID} field.}

\item{mean.h.by}{Character string. Calculate mean hybrid index by group for the group designations in this field. Creates a column called \sQuote{mean_h} 
for ordering data points on the figure using \option{sort.by}.}

\item{sort.by}{Character vector. Fields by which to sort the test subjects by value along the x axis in the figure. Can include \sQuote{mean_h}. 
Default is \sQuote{h_posterior_mode}.}

\item{col.group}{Character string. Variable by which to group colouring of data points.}

\item{group.sep}{Character string. Variable by which to separate groups of data points with vertical grey lines in the plot.}

\item{fill.source}{Logical. Whether to fill the area of the plot representing parental reference samples with a grey background. 
Default \code{FALSE}.}

\item{basic.lines}{Logical. If \code{TRUE}, adds horizontal grey dashed lines at hybrid index 0, 0.5 and 1. Default \code{TRUE}.}

\item{source.col}{Character vector. Colours to be used for the two parental reference set data points.}

\item{source.limits}{Character vector. Colours of horizontal dashed lines marking the innermost credible interval for each 
parental reference set. Used to indicate whether test subjects are clearly outside the hybrid index ranges of parentals.}

\item{custom.abline}{An \code{abline} of the user's choosing.}

\item{...}{Further graphical parameters.}
}
\value{
Alongside the plot, returns a \code{data.table} and \code{data.frame} useful for creating a legend, with 
  the specified \code{test.subject} field, the \code{mean.h.by} field if specified, along with the \code{mean_h} value,
  the \code{Source} field, the colours in a field named \code{col.Dark2}, and the row number for ordering in field \code{rn}.
}
\description{
Plot function for estimated hybrid indices and credible intervals.
}
\examples{

\dontrun{
abc = plot_h(data=hindlabel$hi,
test.subject=hindlabel$test.subject,
mean.h.by="POPID",			                    #Calculate the mean hybrid index for each value of the "POPID" column#
sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid index calculated above and also by 
						                         #individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
col.group="POPID",
group.sep="POPID",
fill.source=TRUE,
basic.lines=FALSE,
source.col=c("blue","red"),
source.limits=c("blue","red"),
custom.abline=abline(h=0.75,col="magenta",lwd=2),
cex=1,pch=16,
cex.lab=1.5,cex.main=1.5,ylim=c(0,1))

#Adding a legend with base R's "legend" function, using the plot_h object created above, "abc":
setkey(abc,rn);		#Order data by row number#
legend("topleft",	#Place the legend in the top left of the figure#
abc[,POPID], 		#Name of the field by which data point colours are grouped#
bg="white",			#Background colour#
text.col=c("black"), #Text colour#
pch=22, 				#Text size#
col=abc[,col.Dark2], #Name of the field containing colour information#
pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
ncol=4,				#Number of columns for splitting the group names#
cex=0.6, pt.cex=0.7)
}
}
\author{
Richard Ian Bailey \email{richardianbailey@gmail.com}
}
