#' Function to split a data.prep analysis object into multiple smaller objects in preparation for parallelized analysis.
#' @param data.prep.object Character string. Name of the \code{data.prep} object created by the \code{data.prep} function.
#' @param splitBy Character string. Name of the column by which to split the \code{data.prep} object. Typically 
#'   \sQuote{INDLABEL} for hybrid index estimation, or \sQuote{locus} for genomic cline estimation.
#' @param keepN Numeric scalar. The number of test subjects desired in each sub file.
#' @details \command{split_data_prep} splits the \code{data.prep} analysis object into multiple files, each with equal 
#'   numbers of test subjects (except for the final file unless the total number of test subjects is a multiple of 
#'   \code{splitBy}). A 'test subject' is e.g. \sQuote{INDLABEL} for hybrid index estimation on each individual, or 
#'   \sQuote{locus} for genomic cline estimation on individual loci. The resulting objects are written to the working 
#'   directory as csv files and removed from the workspace. These can then be loaded into R and used as the 
#'   \code{data.prep.object} in downstream functions. csv files are quite memory-heavy but load quickly.
#' @return No output is returned in the workspace.
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' #Make a set of data.prep files, one for each individual.
#' split_data_prep(
#'  data.prep.object=prepdata$data.prep,   #The data analysis table#
#'  splitBy="INDLABEL",                    #The planned test subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
#'  keepN=1                                #The number of test subjects you want in each resulting file#
#' )
#'
#' #Make another set of files, one for each set of 10 loci.
#' 
#' split_data_prep(
#'  data.prep.object=prepdata$data.prep,   #The data analysis table#
#'  splitBy="locus",                       #The planned test subject (usually "INDLABEL" for hybrid index estimation, "locus" for genomic cline estimation)#
#'  keepN=10                               #The number of test subjects you want in each resulting file 
#'                                         #(the final file will contain fewer if the total is not a multiple of this number)#
#' )
#' 
#' #The entry for the 'splitBy' option will be included in the filename for each resulting csv file. 
#' #Therefore, to create a list of files for analysis e.g. of the 'locus' files above:
#' files=list.files(pattern="_locus_")
#' }
#' @export
split_data_prep=function(data.prep.object,splitBy,keepN){
    setkeyv(data.prep.object,splitBy);
    split.subject=unique(data.prep.object[,..splitBy]);
    Numi=ceiling(split.subject[,.N]/keepN);

    for(i in 1:Numi){
        if(i==Numi){
            prepdata_splitBy = data.prep.object[split.subject[(keepN*(i - 1) + 1):split.subject[,.N]]];
            fwrite(prepdata_splitBy,paste("prepdata_",splitBy,"_",i,".csv",sep=""));
                   }else{
            prepdata_splitBy = data.prep.object[split.subject[(keepN*(i - 1) + 1):(keepN*i)]];
            fwrite(prepdata_splitBy,paste("prepdata_",splitBy,"_",i,".csv",sep=""));
                        };
                    };
    rm(prepdata_splitBy)
}