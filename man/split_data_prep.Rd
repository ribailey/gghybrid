% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_data_prep.r
\name{split_data_prep}
\alias{split_data_prep}
\title{Function to split a data.prep analysis object into multiple smaller objects in preparation for parallelized analysis.}
\usage{
split_data_prep(data.prep.object, splitBy, keepN)
}
\arguments{
\item{data.prep.object}{Character string. Name of the \code{data.prep} object created by the \code{data.prep} function.}

\item{splitBy}{Character string. Name of the column by which to split the \code{data.prep} object. Typically 
\sQuote{INDLABEL} for hybrid index estimation, or \sQuote{locus} for genomic cline estimation.}

\item{keepN}{Numeric scalar. The number of test subjects desired in each sub file.}
}
\value{
No output is returned in the workspace.
}
\description{
Function to split a data.prep analysis object into multiple smaller objects in preparation for parallelized analysis.
}
\details{
\command{split_data_prep} splits the \code{data.prep} analysis object into multiple files, each with equal 
  numbers of test subjects (except for the final file unless the total number of test subjects is a multiple of 
  \code{splitBy}). A 'test subject' is e.g. \sQuote{INDLABEL} for hybrid index estimation on each individual, or 
  \sQuote{locus} for genomic cline estimation on individual loci. The resulting objects are written to the working 
  directory as csv files and removed from the workspace. These can then be loaded into R and used as the 
  \code{data.prep.object} in downstream functions. csv files are quite memory-heavy but load quickly.
}
\examples{

\dontrun{
#Make a set of data.prep files, one for each individual.
split_data_prep(
 data.prep.object=prepdata$data.prep,   #The data analysis table#
 splitBy="INDLABEL",                    #The planned test subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
 keepN=1                                #The number of test subjects you want in each resulting file#
)

#Make another set of files, one for each set of 10 loci.

split_data_prep(
 data.prep.object=prepdata$data.prep,   #The data analysis table#
 splitBy="locus",                       #The planned test subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
 keepN=10                               #The number of test subjects you want in each resulting file 
                                        #(the final file will contain fewer if the total is not a multiple of this number)#
)

#The entry for the 'splitBy' option will be included in the filename for each resulting csv file. 
#Therefore, to create a list of files for analysis e.g. of the 'locus' files above:
files=list.files(pattern="_locus_")
}
}
\author{
Richard Ian Bailey \email{richardianbailey@gmail.com}
}
