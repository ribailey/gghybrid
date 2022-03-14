#' Load data stored in the format of a \href{http://web.stanford.edu/group/pritchardlab/structure.html}{structure} file or similar data table.
#'
#' @param file Character string. The name of the data file to read in.
#' @param mainparams Character string. The name of an associated \file{mainparams} 
#'   file to read in. Optional as \code{mainparams} info can be entered into the 
#'   function call directly (see below).
#' @param extracol.names Character string or vector. Names of extra (i.e. none of those specifically named here) 
#'   non-marker columns in the data file. Optional; used 
#'   when the data file does not contain pre-marker \option{EXTRACOL} headers and you 
#'   wish to add them on input. Default is \code{NULL}.
#' @param precol.headers Numeric. Presence (\code{1}) or absence (\code{0}) of headers for 
#'   pre-marker columns in the header row of the data file. Default is \code{1}. Set to \kbd{0} for \sQuote{plink} structure files.
#' @param nprecol Numeric. Number of pre-marker columns in the data file. No 
#'   default and must always be entered. Set to \kbd{2} for standard \sQuote{plink} structure files.
#' @param markername.dup Numeric. Whether each marker name appears twice in 
#'   the data file header row. Default is \code{0} (FALSE, suitable for \sQuote{plink} structure files), alternative is \code{1} (TRUE).
#' @param NUMLOCI.autoAccept Logical. If the number of loci (\option{NUMLOCI}) is not entered, whether to require manual acceptance/rejection of 
#'   the internally calculated \option{NUMLOCI}. Default is \code{FALSE} with a warning, in order not to interrupt the data analysis 
#'   pipeline. If set to \kbd{TRUE}, the calculated NUMLOCI can be accepted with \kbd{y} or rejected with \kbd{n}, the latter 
#'   of which stops the function with a blank error message and returns no output.
#' @param EXTRACOL Numeric. Number of extra (i.e. none of those 
#'   specifically named here) non-marker columns in the data file. Not needed if this info is 
#'   uploaded in a \file{mainparams} file. For compatibility with \sQuote{structure} \code{mainparams} info. Default is \code{0} (suitable for \sQuote{plink} structure files).
#' @param INDLABEL Numeric. Presence (\code{1}) or absence (\code{0}) of a column of 
#'   individual references in the data file. Not needed if the info is uploaded in a \file{mainparams} file. Default is \code{1}, and this column 
#'   must always be present. See \sQuote{Details}. The first two columns of \sQuote{plink} structure files can be treated as \option{INDLABEL} and 
#'   \option{POPID} columns.
#' @param LOCDATA Numeric. Presence (\code{1}) or absence (\code{0}) of a \option{LOCDATA} column in 
#'   the data file. For compatibility with \sQuote{structure}. Not needed if the info is uploaded in a \file{mainparams} file. Default is \code{0}.
#' @param MAPDISTANCE Numeric. Presence (\code{1}) or absence (\code{0}) of a \code{MAPDISTANCE} 
#'   row in the data file, which will be removed if present. Not needed if the info is uploaded in a \file{mainparams} file. 
#'   Set to \kbd{1} for \sQuote{plink} structure files that include this row below the marker names. Default is \code{0}.
#' @param MARKERNAME Numeric. Presence (\code{1}) or absence (\code{0}) of a header row 
#'   containing marker names in the data file. Not needed if the info is uploaded in a \file{mainparams} file. Default is \code{1} and
#'   is suitable for \sQuote{plink} structure files.
#' @param MISSINGVAL The identifier for missing values. No default. Not needed if the info is uploaded in a \file{mainparams} file. Currently,
#'   \code{0} is not allowed and hence \sQuote{plink} structure files need to be modified before reading in, preferably replacing \code{0} with \code{NA} (no quotation marks needed).
#' @param NUMINDS Numeric. Number of individuals in the data set. Required to check that the marker data are of the appropriate dimensions. Default 
#'   set to \code{0}.
#' @param NUMLOCI Numeric. Number of loci in the data set. Not needed if the info is uploaded in a \file{mainparams} file. If the default value (\code{0}) 
#'   is used, the user has the option to manually accept/reject the internally calculated number of loci (calculated using \option{NUMLOCI}
#'   , \option{PLOIDY} and the dimensions of the marker data) using the \option{NUMLOCI.autoAccept} argument above, as a means of error-checking. If 
#'   the input file is a structure file containing phasing rows (\code{PHASED=1}), the correct \option{NUMLOCI} must be entered.
#' @param ONEROW Numeric. Whether (\code{1}) or not (\code{0}) there is a single data row 
#'   per \option{INDLABEL}. If the value is \code{1}, the number of columns per marker = 
#'   \option{PLOIDY}, which is suitable for \sQuote{plink} structure files. If \code{ONEROW = 0}, there is one column per marker and the number 
#'   of rows per \option{INDLABEL} = \option{PLOIDY}. Not needed if the info is uploaded in a \file{mainparams} file. Default is \code{0}.
#' @param PHASED Numeric. Presence (\code{1}) or absence (\code{0}) of phasing data rows 
#'   in a \sQuote{structure} input file. These rows will be removed if present. Not needed if the info is uploaded in a \file{mainparams} file. 
#'   Default is \code{0}.
#' @param PHENOTYPE Numeric. Presence (\code{1}) or absence (\code{0}) of a \code{PHENOTYPE} 
#'   column in the data file. For compatibility with \sQuote{structure}. Not needed if the info is uploaded in a \file{mainparams} file. Default 
#'   is \code{0}.
#' @param PLOIDY Numeric. Maximum ploidy among the markers in the data file. Not needed if the info is uploaded in a \file{mainparams} file. 
#'   Default is \code{2}.
#' @param POPFLAG Numeric. \code{0} or \code{1}, for compatibility with \sQuote{structure}. Default \code{0}.
#' @param POPID Numeric. Presence (\code{1}) or absence (\code{0}) of a column of population 
#'   identifiers in the data file. Not needed if the info is uploaded in a \file{mainparams} file. 
#'   Default is \code{1}. The first two columns of \sQuote{plink} structure files can be treated as \option{INDLABEL} and \option{POPID} columns respectively.
#' @param RECESSIVEALLELE Numeric. \code{0} or \code{1}, for compatibility with \sQuote{structure}. Default \code{0}.
#' @param marker.info.file Character string. The name of an optional file of additional marker information to be read in, with one row per marker and any number of columns.  
#'   This can be joined to the primary output table (and to the locus table if \code{return.locus.table} is set to \code{TRUE}) when using the \command{data.prep} function.
#'   The file should contain parental allele frequencies when these are not to be calculated from the data. In that case, it 
#'   must contain columns entitled \code{refAllele}, \code{alternateAllele}, \code{S0.prop_r} and \code{S1.prop_r}, as well as 
#'   the obligatory \code{locus} column. The first two columns identify an arbitrarily chosen reference allele and the alternate 
#'   allele, while \code{S0.prop_r} and \code{S1.prop_r} indicate the allele frequency of the chosen reference allele in each source 
#'   (parental reference) population.
#'   Default is \code{NULL}.
#' @param sourceAbsent Logical. Whether source (parental reference) samples are absent from the dataset. If \code{TRUE}, the parental 
#'   allele frequencies of a reference allele for each locus must be present in the \option{marker.info.file} in columns named 
#'   \code{S0.prop_r} and \code{S1.prop_r}. Default is \code{FALSE}.
#' @details \command{read.data} is designed to be compatible with data input files 
#'   for the software \href{http://web.stanford.edu/group/pritchardlab/structure.html}{structure}, including those produced by 
#'   \href{https://www.cog-genomics.org/plink2/}{plink}. However, \sQuote{plink} structure files use zero to denote missing data,
#'   which is not currently allowed in the \href{https://cran.r-project.org/web/packages/data.table/index.html}{data.table} 
#'   package, although this limitation is currently being actioned. For the time being, zeroes should be replaced prior to use of the \command{read.data} function, 
#'   preferably with \kbd{NA} (no quotations needed).
#'
#' The simplest file format to read in is
#'   a rectangular data table with a complete header row, one column per marker and one row per allele copy (equivalent to \code{ONEROW
#'   = 0}). With this format the only required fields other than those with default settings  
#'   are \option{file}, \option{nprecol}, \option{NUMINDS} and \option{MISSINGVAL}. Furthermore, allele identities are 
#'   read in as character strings and therefore each allele can be given any 
#'   name or number, as long as there are only two unique alleles per locus. 
#'   As with \href{http://web.stanford.edu/group/pritchardlab/structure.html}{structure}, all non-marker 
#'   columns should be to the left of the marker columns in the data file.
#'
#' The data input file is required to have at least one column prior to the 
#'   marker data columns. Typically there are two columns, the first (required) is referred to as \sQuote{INDLABEL} (but can be 
#'   given a different column header) and represents the finest resolution of 
#'   identification that may be required for the data set (usually a unique 
#'   individual reference, but could for example be a population-level 
#'   reference for poolseq data). The second is referred to as \sQuote{POPID} and 
#'   indicates the population from which the individual was sampled. The main 
#'   purpose of the \sQuote{POPID} column is to declare, downstream, which \option{POPID} values 
#'   represent the parental reference samples, and it can also be useful for 
#'   plotting hybrid index estimates. If the \sQuote{POPID} column is absent or if it is desired to use another column to 
#'   identify parental reference samples, an alternative column must be identified downstream in the \command{data.prep} function.
#'
#' Any ploidy is allowed, but the declared \option{PLOIDY} must be the same for all markers. So for example if 
#'   data are haplodiploid, haploid markers should either be present as 
#'   diploid homozygotes, or the second allele declared as missing data. 
#'
#' If an associated \file{mainparams} file is read in, it should contain two columns 
#'   and no header row, with the first column holding the field names and the 
#'   second the field values. The only mandatory fields are  
#'   (and must have one of these synonyms): \code{PLOIDY}, 
#'   \code{MISSINGVAL}/\code{MISSING}, \code{ONEROW}/\code{ONEROWPERIND}, \code{INDLABEL}, \code{POPID}, 
#'   \code{EXTRACOL}/\code{EXTRACOLS}, \code{MARKERNAME}/\code{MARKERNAMES}. If \code{RECESSIVEALLELE} or 
#'   \code{MAPDISTANCE} are present, they can be pluralized.
#'
#' \command{read.data} uses the \code{\link[data.table]{fread}} function from \pkg{\href{https://cran.r-project.org/web/packages/data.table/index.html}{data.table}} 
#'   to rapidly read in the data file. The 
#'   loaded data file is therefore of \href{http://127.0.0.1:23982/library/base/html/class.html}{class} \code{data.table} and \code{data.frame}. 
#'   A \code{data.table} can be treated the same as a standard \code{data.frame}
#'   for those not familiar with the package.
#' @return \code{read.data} returns a list with one component containing the loaded data in \code{ONEROW = 0} format, in the 
#'   form of a \code{data.table} and \code{data.frame}, and other components used in downstream functions or otherwise potentially 
#'   useful to the user.
#'
#'   The list contains the following components:
#'   \item{mainparams}{A \code{data.table} and \code{data.frame} with the inputted \code{mainparams} information (or default values if 
#'   no information provided).}
#'   \item{nprecols}{Numeric. The number of non-marker columns to the left of the marker columns in the imported data set.}
#'   \item{precols}{Character vector. The names of non-marker columns to the left of the marker columns in the imported data set.}
#'   \item{data}{A \code{data.table} and \code{data.frame}. The imported data.}
#'   \item{alleles}{A \code{data.table} and \code{data.frame}. The names of the two alleles at each locus.}
#'   \item{loci}{Character vector. The names of all loci in the imported data set.}
#'   \item{marker.info}{A \code{data.table} and \code{data.frame}. The imported marker.info.file.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' #First create example data files and save to the working directory#
#'
#' #1. Genomic (SNP) data as a regular table with full header row#
#' ex <- "INDLABEL,POPID,chr1:001,chr1:002\nind1,pop1,A,A\nind1,pop1,A,B\nind2,pop1,B,B\nind2,pop1,B,B\nind3,pop2,NA,A\nind3,pop2,B,A\n"
#'
#' #Save 'ex' to the working directory as 'ex.data', to be read in by read.data. Open in e.g. notepad or notepad++ to view in table format#
#' cat("INDLABEL POPID chr1:001 chr1:002","ind1 pop1 A A","ind1 pop1 A B","ind2 pop1 B B","ind2 pop1 B B","ind3 pop2 NA A","ind3 pop2 B A", file = "ex.data", sep = "\n")
#'
#' #2. Now create and save a file, ex2.data, in plink structure format (but with missing data recoded from 0 to NA), with a MAPDISTANCE row present below the marker names#
#' cat("chr1:001 chr1:002","-3 3 2 4","ind1 pop1 A A A B","ind2 pop1 B B B B","ind3 pop2 NA B A A", file = "ex2.data", sep = "\n")
#'
#' #Load ex.data, in this case without specifying the number of loci (useful if this is not known in advance). These are the minimum arguments required for a file in this format#
#' dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,MISSINGVAL=NA)#NUMLOCI is not specified so will be calculated, with a warning#
#'
#' #The same but specifying the number of loci. The resulting object is identical#
#' dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA)#No warning this time as NUMLOCI is specified#
#'
#' #The minimum arguments required to read in the PLINK structure format file with missing data recoded as 'NA' (gives a warning because NUMLOCI is not specified)#
#' dat2 <- read.data(file="ex2.data",MISSINGVAL=NA,NUMINDS=3,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
#'
#' #Same, but displaying all arguments including those where the default works for this file format#
#' dat2 <- read.data(
#'   file="ex2.data",
#'   mainparams = NULL,              #Default#
#'   extracol.names = NULL,          #Default#
#'   precol.headers = 0,
#'   nprecol=2,
#'   markername.dup = 0,             #Default#
#'   NUMLOCI.autoAccept = TRUE,      #Default#
#'   EXTRACOL = 0,                   #Default#
#'   INDLABEL = 1,                   #Default#
#'   LOCDATA = 0,                    #Default#
#'   MAPDISTANCE = 1,
#'   MARKERNAME = 1,                 #Default#
#'   MISSINGVAL = NA,                #Cannot be zero for the time being#
#'   NUMINDS = 3,
#'   NUMLOCI = 2,                    #Only mandatory when PHASED=1; otherwise will be calculated internally if not known#
#'   ONEROW = 1,
#'   PHASED = 0,                     #Default#
#'   PHENOTYPE = 0,                  #Default#
#'   PLOIDY = 2,                     #Default#
#'   POPFLAG = 0,                    #Default#
#'   POPID = 1,                      #Default#
#'   RECESSIVEALLELE = 0,            #Default#
#'   marker.info.file = NULL,        #Default#
#'   sourceAbsent = FALSE            #Default#
#' )
#'
#' #If NUMLOCI is not specified, the number of loci will be calculated with a warning, which will not interfere with downstream processes. 
#' #However, if 'NUMLOCI.autoAccept = FALSE' is set, the user is required to manually accept the calculated number of loci. This option is 
#' #included in case the user wishes to verify that the calculated number of loci is accurate. Stops with an error if the estimated NUMLOCI 
#' #is manually rejected. Example:
#' dat2 <- read.data(file="ex2.data",NUMLOCI.autoAccept = FALSE,MISSINGVAL=NA,NUMINDS=3,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
#'
#' #Uploading a marker info file.
#'
#' #Example: Create a marker info file indicating whether the locus is intronic or exonic, and save to the working directory#
#' cat("locus type","chr1:001 intronic","chr1:002 exonic", file = "ex_marker_info.data", sep = "\n")
#'
#' #Read it in alongside the SNP data#
#' dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info.data")
#'
#' #Loading a data file that does not include parental reference samples. For the situation where parental (S0 and S1) allele frequencies are 
#' #to be loaded, but not parental reference samples. A marker.info.file can be loaded regardless of whether S0 and S1 samples are present in 
#' #the dataset, but it is obligatory if they are absent. When 'sourceAbsent = TRUE', as a minimum the following columns with the exact headers 
#' #in the first set of quotation marks below are required in the marker.info.file (more columns are allowed).
#' cat("locus refAllele alternateAllele S0.prop_r S1.prop_r type","chr1:001 A B 0.1 0.9 intronic","chr1:002 B A 0.8 0.2 exonic", file = "ex_marker_info2.data", sep = "\n")
#'
#' #"...prop_r" means the allele frequency of the reference allele. Choice of reference and alternate alleles is arbitrary#
#' dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info2.data", sourceAbsent = TRUE)
#'
#' #The data.prep function will then determine which allele has higher frequency in S1, as it would if parental reference samples were included#
#'
#' unlink("ex_marker_info.data")#Tidy up#
#' unlink("ex_marker_info2.data")#Tidy up#
#' unlink("ex.data")#Tidy up#
#' unlink("ex2.data")#Tidy up#
#' }
#' @export
read.data=function(file,mainparams=NULL,extracol.names=NULL,
    precol.headers=1,nprecol,markername.dup=0,NUMLOCI.autoAccept=TRUE,
    EXTRACOL=0,INDLABEL=1,LOCDATA=0,MAPDISTANCE=0,MARKERNAME=1,MISSINGVAL,
    NUMINDS=0,NUMLOCI=0,ONEROW=0,PHASED=0,PHENOTYPE=0,PLOIDY=2,POPFLAG=0,
    POPID=1,RECESSIVEALLELE=0,marker.info.file=NULL,sourceAbsent=FALSE){
 
    cat(paste("Reading in data file..."),fill=1); flush.console();
	
#Reading in the marker.info.file#######################################BELOW	
    if(is.null(marker.info.file)==FALSE){

        marker.info=fread(marker.info.file);#********REMEMBER THIS CAN INCLUDE S0 AND S1 ALLELE FRQs***#
		
		if(sourceAbsent==TRUE){
		    if("S0.prop_r"%in%names(marker.info)==FALSE){
			    stop("If there are no source (parental reference) individuals in the dataset with which to determine parental allele frequencies, 
the frequency of a reference allele in each source population (S0 and S1) must be present in the marker.info.file, in columns 
named 'S0.prop_r' and 'S1.prop_r'. S0 is the source population to be assigned hybrid index = 0; S1 is the source population with 
hybrid index = 1");
			                                            };
														   
		    if("S1.prop_r"%in%names(marker.info)==FALSE){
			    stop("If there are no source (parental reference) individuals in the dataset with which to determine parental allele frequencies, 
the frequency of a reference allele in each source population (S0 and S1) must be present in the marker.info.file, in columns 
named 'S0.prop_r' and 'S1.prop_r'. S0 is the source population to be assigned hybrid index = 0; S1 is the source population with 
hybrid index = 1");
			                                            };
		    if("refAllele"%in%names(marker.info)==FALSE){
			    stop("If there are no source (parental reference) individuals in the dataset, 
the identity of a reference allele and the alternate allele must be declared in the marker.info.file, in columns 
named 'refAllele' and 'alternateAllele'. The choice of reference allele is arbitrary.");
			                                            };
														   
		    if("alternateAllele"%in%names(marker.info)==FALSE){
			    stop("If there are no source (parental reference) individuals in the dataset, 
the identity of a reference allele and the alternate allele must be declared in the marker.info.file, in columns 
named 'refAllele' and 'alternateAllele'. The choice of reference allele is arbitrary.");
			                                                  };														
            marker.info[,refAllele:=as.character(refAllele)];
            marker.info[,alternateAllele:=as.character(alternateAllele)];

            allele.id=data.table(t(marker.info[,.(refAllele,alternateAllele)]));
            setnames(allele.id,marker.info[,locus]);		
		                      };

                                        };
#Reading in the marker.info.file#######################################ABOVE	 
	 
#Creating the mainparams object##########################################BELOW
	if(is.null(mainparams)==FALSE){

    mainparams=fread(mainparams,colClasses=c("character","character"),
	    header=FALSE);
    setkey(mainparams,V1);
    mainparams["ONEROWPERIND",V1:="ONEROW"];
    setkey(mainparams,V1);
    mainparams["MISSING",V1:="MISSINGVAL"];
    setkey(mainparams,V1);
    mainparams["EXTRACOLS",V1:="EXTRACOL"];
    setkey(mainparams,V1);
    mainparams["RECESSIVEALLELES",V1:="RECESSIVEALLELE"];
    setkey(mainparams,V1);
    mainparams["MARKERNAMES",V1:="MARKERNAME"];
    setkey(mainparams,V1);
    mainparams["MAPDISTANCES",V1:="MAPDISTANCE"];
    
	mp2=data.table(V1=c("markername.dup","precol.headers"),
	    V2=as.character(c(markername.dup,precol.headers)));
	setkey(mainparams);setkey(mp2);
	mainparams=rbind(mainparams,mp2);
	mp2=NULL;
                              }else{
    mainparams=data.table(V1=c("EXTRACOL","INDLABEL","LOCDATA","MAPDISTANCE",
        "MARKERNAME","markername.dup","MISSINGVAL","NUMINDS","NUMLOCI",
		"ONEROW","PHASED","PHENOTYPE","PLOIDY","POPFLAG","POPID",
		"precol.headers","RECESSIVEALLELE"),
        V2=as.character(c(EXTRACOL,INDLABEL,LOCDATA,MAPDISTANCE,MARKERNAME,
        markername.dup,MISSINGVAL,NUMINDS,NUMLOCI,ONEROW,PHASED,PHENOTYPE,
        PLOIDY,POPFLAG,POPID,precol.headers,RECESSIVEALLELE)));
                                    };									
#Creating the mainparams object##########################################ABOVE

#NUMINDS, nprecol and NUMLOCI (PHASED=1) error check#####################BELOW
    setkey(mainparams,V1);
    if(mainparams["NUMINDS",V2=="0"]){
        stop("Please enter NUMINDS.");
                                       };
									   
    if(mainparams["PHASED",V2=="1"]&mainparams["NUMLOCI",V2=="0"]){
        stop("Please enter NUMLOCI.");
                                                                    };
																	
    nprecolumns=mainparams[c("INDLABEL","POPID",
	    "POPFLAG","LOCDATA","PHENOTYPE","EXTRACOL"),sum(as.numeric(V2))];
		
	if(nprecol!=nprecolumns){
	    stop("There is a discrepancy between nprecol and the number of pre-marker 
columns specified by arguments INDLABEL, POPID, POPFLAG, LOCDATA, PHENOTYPE and EXTRACOL.")
		                    };
#NUMINDS, nprecol and NUMLOCI (PHASED=1) error check#####################ABOVE

#Preparations for reading in data########################################BELOW
    ntop=mainparams[c("MARKERNAME","RECESSIVEALLELE","MAPDISTANCE"),
	    sum(as.numeric(V2))];
    precols=subset(mainparams[c("INDLABEL","POPID",
	    "POPFLAG","LOCDATA","PHENOTYPE")],V2=="1");
	precols=precols[,as.character(V1)];
#Preparations for reading in data########################################ABOVE

#Reading in data for a simple table format with all headers present and PLOIDY rows per locus#BELOW
    if(sum(mainparams[c("MAPDISTANCE","MARKERNAME","markername.dup",
        "ONEROW","precol.headers","RECESSIVEALLELE"),
        V2==c("0","1","0","0","1","0")]) == 6
       ){
        data=fread(file,header=TRUE,
		    na.strings=paste(mainparams["MISSINGVAL",V2]),
            colClasses="character");

        MARKER=names(data[,(nprecol+1):ncol(data),with=F]);#SEEMS LIKE with=F IS NO LONGER NECESSARY************************

        if(mainparams["PHASED",V2=="1"]){
		    data = data[rep(c(rep(TRUE,mainparams["PLOIDY",as.numeric(V2)]),
			    FALSE),length=.N)];
			                            };
										
#******************************************************************************
#Error check###################################################################
    if(mainparams["ONEROW",V2=="1"]){
	    if(mainparams["NUMINDS",as.numeric(V2)]!=data[,.N]){
		    stop("Number of data rows does not match NUMINDS.");
		                                                   };
		                          };

    if(mainparams["ONEROW",V2=="0"]){
	    if(mainparams["NUMINDS",as.numeric(V2)]*mainparams["PLOIDY",as.numeric(V2)]!=data[,.N]){
		    stop("Number of data rows does not match NUMINDS.");
		                                                                                       };
		                          };
								  
    if(mainparams["NUMLOCI",V2=="0"]){
#CALCULATE NUMLOCI FROM THE DATA DIMENSIONS AND ADD A WARNING.
        mainparams["NUMLOCI",V2:=prod(dim(data[,(nprecol+1):ncol(data)]))/prod(
		    as.numeric(mainparams[V1%in%c("NUMINDS","PLOIDY")]$V2))];
        if(NUMLOCI.autoAccept==FALSE){
            ANSWER=readline(paste0("NUMLOCI calculated as ",dat$mainparams[V1=="NUMLOCI",
			    as.numeric(V2)],", does this look right to you (y/n)? "));
            if (substr(ANSWER, 1, 1) == "n"){
			    stop();
				                            };
                                     }else{
        warning("NUMLOCI not entered, therefore calculated using NUMINDS, PLOIDY and nprecol along with data dimensions. 
		    Check these entries are correct before proceeding.");
                                          };
                                       }else{
        if(prod(dim(data[,(nprecol+1):ncol(data)]))!=
	        prod(as.numeric(mainparams[V1%in%c("NUMINDS","NUMLOCI","PLOIDY")]$V2))){
                stop("Data file dimensions do not match NUMINDS*NUMLOCI*PLOIDY. 
Check these entries and also nprecol are correct.");
                                                                                   };									   
									        };										
#Error check###################################################################
#******************************************************************************

		if(sourceAbsent==FALSE){
            allele.id=data[,lapply(.SD,unique.fun),.SDcols=MARKER[]];
                               };

        output=list();
            output$mainparams=mainparams;
            output$nprecols=nprecol;
            output$precols=names(data[,1:nprecol])
            output$data=data;
            output$alleles=allele.id;
            output$loci=MARKER;
			
        if(is.null(marker.info.file)==FALSE){
            output$marker.info=marker.info;			
                                            };
											
        cat(paste("Done"),fill=1);

        return(output)
#Reading in data for a simple table format with all headers present and PLOIDY rows per locus#ABOVE
		
         }else{
        data=fread(file,header=FALSE,skip=ntop,
		    na.strings=paste(mainparams["MISSINGVAL",V2]),
			colClasses="character");

if(mainparams["PHASED",V2=="0"]){
#******************************************************************************
#Error check###################################################################
    if(mainparams["ONEROW",V2=="1"]){
	    if(mainparams["NUMINDS",as.numeric(V2)]!=data[,.N]){
		    stop("Number of data rows does not match NUMINDS.");
		                                                   };
		                          };

    if(mainparams["ONEROW",V2=="0"]){
	    if(mainparams["NUMINDS",as.numeric(V2)]*mainparams["PLOIDY",as.numeric(V2)]!=data[,.N]){
		    stop("Number of data rows does not match NUMINDS.");
		                                                                                       };
		                          };

    if(mainparams[V1=="NUMLOCI",V2=="0"]){
#CALCULATE NUMLOCI FROM THE DATA DIMENSIONS AND ADD A WARNING.
        mainparams[V1=="NUMLOCI",V2:=prod(dim(data[,(nprecol+1):ncol(data)]))/prod(
		    as.numeric(mainparams[V1%in%c("NUMINDS","PLOIDY")]$V2))];
        if(NUMLOCI.autoAccept==FALSE){
            ANSWER=readline(paste0("NUMLOCI calculated as ",dat$mainparams[V1=="NUMLOCI",
			    as.numeric(V2)],", does this look right to you (y/n)? "));
            if (substr(ANSWER, 1, 1) == "n"){
			    stop();
				                            };
                                     }else{
        warning("NUMLOCI not entered, therefore calculated using NUMINDS, PLOIDY and nprecol along with data dimensions. 
		    Check these entries are correct before proceeding.");
                                          };
                                       }else{
        if(prod(dim(data[,(nprecol+1):ncol(data)]))!=
	        prod(as.numeric(mainparams[V1%in%c("NUMINDS","NUMLOCI","PLOIDY")]$V2))){
                stop("Data file dimensions do not match NUMINDS*NUMLOCI*PLOIDY. 
Check these entries and also nprecol are correct.");
                                                                                   };									   
									        };										
#Error check###################################################################
#******************************************************************************
                                };

#***NEED TO CALCULATE NUMLOCI BEFORE THIS POINT***#
        if(mainparams["MARKERNAME",V2=="0"]){
	        MARKER=make.unique(rep("L",mainparams["NUMLOCI",
			    as.numeric(V2)+1]))[-1];
		                                    };

        if(sum(mainparams[c("MARKERNAME","markername.dup","precol.headers"),
            V2==c("1","0","0")]) == 3){MARKER=
	        scan(file,nmax=mainparams["NUMLOCI",as.numeric(V2)],
			    what=character(),quiet=TRUE);
		                              };
								  
        if(sum(mainparams[c("MARKERNAME","markername.dup","precol.headers"),
            V2==c("1","0","1")]) == 3){
			cnames=scan(file,nmax=(mainparams["NUMLOCI",as.numeric(V2)] + 
			    nprecol),what=character(),quiet=TRUE);	
		    MARKER=cnames[(nprecol + 1):length(cnames)];
		                              };
									  
        if(sum(mainparams[c("markername.dup","precol.headers"),
            V2==c("1","0")]) == 2){
			MARKER=unique(scan(file,nmax=(mainparams["NUMLOCI",
			    as.numeric(V2)]*mainparams["PLOIDY",as.numeric(V2)]),
			    what=character(),quiet=TRUE));
				                  };
			
        if(sum(mainparams[c("markername.dup","precol.headers"),
            V2==c("1","1")]) == 2){
			cnames=unique(scan(file,nmax=(mainparams["NUMLOCI",
			    as.numeric(V2)]*mainparams["PLOIDY",as.numeric(V2)] + 
				nprecol),what=character(),quiet=TRUE));	
            MARKER=cnames[(nprecol + 1):length(cnames)];
			                      };

        if(mainparams["precol.headers",V2=="0"]){
            if(is.null(extracol.names)&mainparams["EXTRACOL",V2!="0"]){
                EXTRACOLN=make.unique(rep("V",mainparams["EXTRACOL",
				    as.numeric(V2) + 1]))[-1];
                cnames=c(precols,EXTRACOLN,MARKER);
		                                                              }else{
                cnames=c(precols,extracol.names,MARKER);
		                                    }                                        
										        };
	
        if(mainparams["ONEROW",V2=="0"]){
            setnames(data,cnames);

            if(mainparams["PHASED",V2=="1"]){
			    data = data[rep(c(rep(TRUE,mainparams["PLOIDY",
				    as.numeric(V2)]),FALSE),length=.N)];
                if(prod(dim(data[,(nprecol+1):ncol(data)]))!=
	                prod(as.numeric(mainparams[V1%in%c("NUMINDS","NUMLOCI","PLOIDY")]$V2))){
                        stop("Data file dimensions do not match NUMINDS*NUMLOCI*PLOIDY. 
Check these entries and also nprecol are correct.");
                                                                                           };
				                            };

		if(sourceAbsent==FALSE){
            allele.id=data[,lapply(.SD,unique.fun),.SDcols=MARKER[]];
                               };

            output<-list();
                output$mainparams=mainparams;
                output$nprecols=nprecol;
                output$precols=cnames[1:nprecol]
                output$data=data;
                output$alleles=allele.id;
                output$loci=MARKER;
				
            if(is.null(marker.info.file)==FALSE){
                output$marker.info=marker.info;			
                                                };				

            cat(paste("Done"),fill=1);

            return(output)									
									    }else{
            if(mainparams["PHASED",V2=="1"]){
			    data = data[rep(c(TRUE,FALSE),length=.N)];
                if(prod(dim(data[,(nprecol+1):ncol(data)]))!=
	                prod(as.numeric(mainparams[V1%in%c("NUMINDS","NUMLOCI","PLOIDY")]$V2))){
                        stop("Data file dimensions do not match NUMINDS*NUMLOCI*PLOIDY. 
Check these entries and also nprecol are correct.");
                                                                                           };
				                            };

            colNum=c((1:nprecol),
                seq(
				    (nprecol + 1),
                    (mainparams["NUMLOCI",as.numeric(V2)]*mainparams["PLOIDY",
					as.numeric(V2)] + nprecol - (mainparams["PLOIDY",
					as.numeric(V2)]) + 1),by=mainparams["PLOIDY",
					as.numeric(V2)]
				   )
			        );

            dat=data[0,colNum,with=F];setnames(dat,cnames);setkey(dat);
						
            for(i in 1:mainparams["PLOIDY",as.numeric(V2)]){
                colNum=seq(
				    (nprecol + i),(mainparams["NUMLOCI",
					as.numeric(V2)]*mainparams["PLOIDY",
					as.numeric(V2)] + nprecol - (mainparams["PLOIDY",
					as.numeric(V2)]) + i),
                    by=mainparams["PLOIDY",as.numeric(V2)]
					      );
                dat1=data[, .SD, .SDcols=colNum];
                dat1=cbind(data[,1:nprecol],dat1);
                setnames(dat1,cnames);setkey(dat1);
                dat=rbind(dat,dat1);
                                                           };
												   
            data=NULL;

		    if(sourceAbsent==FALSE){
                allele.id=dat[,lapply(.SD,unique.fun),.SDcols=MARKER[]];
                                   };            
			
            output=list();
                output$mainparams=mainparams;
                output$nprecols=nprecol;
                output$precols=cnames[1:nprecol]
                output$data=dat;
                output$alleles=allele.id;
                output$loci=MARKER;
				
            if(is.null(marker.info.file)==FALSE){
                output$marker.info=marker.info;			
                                                };				

            cat(paste("Done"),fill=1);

            return(output)		
	                                         };
               };
                     }
