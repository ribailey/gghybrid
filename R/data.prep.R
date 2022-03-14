#' Prepare data for hybrid index and genomic cline estimation.
#'
#' @param data A \code{data} object produced by \code{read.data}, or a custom \code{data.table} in the same format.
#' @param loci A \code{loci} object produced by \code{read.data}, or a custom character vector with locus names.
#' @param alleles An \code{alleles} object produced by \code{read.data}, or a custom \code{data.table} in the same format.
#' @param marker.info A \code{marker.info} object uploaded using \code{read.data}, or a custom \code{data.table} in the same format, 
#'   with one row per locus and any number of columns. Default is \code{NULL}.
#' @param S0 Character vector. Identifiers of the test subjects representing the first parental reference set (\dQuote{Source 0}). 
#'   Allele frequencies in this set represent a hybrid index of 0. By default \code{data.prep} searches for these identifiers in the 
#'   \code{POPID} column of the \code{data} object, but a different column can be designated using \option{POPID.name}. Default is 
#'   \code{NULL}, and parental reference samples are not needed if parental allele frequencies are included in a \code{marker.info} 
#'   file in the correct format.
#' @param S1 Character vector. Identifiers of the test subjects representing the second parental reference set (\dQuote{Source 1}). 
#'   Allele frequencies in this set represent a hybrid index of 1. By default \code{data.prep} searches for these identifiers in the 
#'   \code{POPID} column of the \code{data} object, but a different column can be designated using \option{POPID.name}.  Default is 
#'   \code{NULL}, and parental reference samples are not needed if parental allele frequencies are included in a \code{marker.info} 
#'   file in the correct format.
#' @param sourceAbsent Logical. Whether source (parental reference) samples are absent from the dataset. If \code{TRUE}, the parental 
#'   allele frequencies will be taken from the \code{marker.info} object uploaded using \code{read.data}. 
#'   Default is \code{FALSE}.
#' @param precols Character vector. A \code{precols} object produced by \code{read.data}, or a custom character vector with 
#'   names of all pre-marker columns in the \code{data} object. No default.
#' @param INDLABEL.name Character string. The name of the column in \code{data} designated as \code{INDLABEL}. 
#'   Default is \dQuote{INDLABEL}.
#' @param POPID.name Character string. The name of the column in \code{data} designated as \code{POPID}. If the \code{data} object 
#'   contains only one pre-marker column, such as \code{INDLABEL}, this must be specified as the \option{POPID.name} to avoid an 
#'   error. The main purpose of this option is to declare which column contains the identifiers for parental reference samples, and 
#'   any non-marker column can be used. Although it must be specified at least by the default, it is ignored if there are no parental 
#'   reference samples in the dataset. Default is \dQuote{POPID}.
#' @param max.S.MAF Numeric. Locus filtering option. Removes loci for which the smaller of the two parental minor allele frequencies 
#'   is greater than this value. Included because only one parental reference minor allele frequency needs to be close to 
#'   zero for the locus to be informative. Default is \code{0.5} (no filtering).
#' @param min.diff Numeric. Locus filtering option. Removes loci for which the difference in parental reference allele frequency 
#'   is less than this value. Default is \code{0} (no filtering).
#' @param min.allele.copies.S0 Numeric. Locus filtering option based on sample size. Removes loci for which the number of allele 
#'   copies in the \code{S0} parental reference set, excluding missing data, is less than this value. This option is only available 
#'   when parental reference samples \code{S0} and \code{S1} are included in the dataset. Default is \code{0} (no filtering). 
#' @param min.allele.copies.S1 Numeric. Locus filtering option based on sample size. Removes loci for which the number of allele 
#'   copies in the \code{S1} parental reference set, excluding missing data, is less than this value. This option is only available 
#'   when parental reference samples \code{S0} and \code{S1} are included in the dataset. Default is \code{0} (no filtering). 
#' @param AF.CIoverlap Logical. Recommended locus filtering option, based on posterior credible intervals of parental allele frequencies. 
#'   If \code{FALSE}, removes loci for which there is overlap between the allele frequency 95 percent credible intervals of the 
#'   two parental reference populations. If parental reference sets are identified, the posterior credible intervals are calculated 
#'   internally. If this option is to be used on parental allele frequencies uploaded in a \code{marker.info} file, that file must 
#'   contain columns entitled \sQuote{S0.prop_1_cred.int.upper} and \sQuote{S1.prop_1_cred.int.lower}, to indicate some kind of 
#'   confidence limits for the parental allele frequencies. Only the upper limit is needed for \code{S0}, and only the lower limit is 
#'   needed for \code{S1}. Default is \code{TRUE} (no filtering).
#' @param return.genotype.table Logical. Whether to create a genotype table with individual genotypes for each locus scored from 
#'   \code{0} to \code{PLOIDY}, based on their number of copies of the designated \code{S1} allele (the allele with relatively higher 
#'   frequency in the \code{S1} parental reference set). Default is \code{FALSE}.
#' @param return.locus.table Logical. Whether to return a table with one row per locus containing a wider array of locus-level data, 
#'   such as parental allele frequency credible intervals and number of non-missing allele copies in each parental reference set. 
#'   Default is \code{FALSE} to save memory, but this could be useful in many circumstances.
#' @param prior.prop_1 Numeric vector. Values of \code{shape1} and \code{shape2} parameters for a beta-distributed prior
#'   for parental allele frequencies (proportions) of the \code{S1_allele}. Default is \code{shape1=shape2=1} (uniform distribution).
#' @details The primary output of \code{data.prep} is a long-form \code{data.table} with 1 row per non-missing allele copy in the 
#'   data set, to be used in downstream hybrid index and genomic cline estimation. It also optionally produces (1) a locus table,
#'   with one row per locus and containing useful locus-level data (see \sQuote{Value}), and (2) a genotype table, with genotypes 
#'   scored according to the number of alleles from parental reference set \code{S1}. This may be useful for example in estimating 
#'   admixture-induced (\sQuote{ancestry}) linkage disequilibria. \code{data.prep} uses objects produced by \code{read.data} as input, 
#'   or custom objects in the same format.
#' 
#' Even when all filtering options are left at default values, loci with exactly zero allele frequency difference between S0 and S1 will 
#'   be removed. Therefore, it is often likely that the resulting analysis table will contain fewer loci than were loaded by \code{read.data}.
#'
#' The optional \code{marker.info} file should be a file with one row per marker and any number of columns (including the marker name in 
#'   a column named \sQuote{locus}), each for a variable of interest for grouping markers, such as chromosome, map location, gene 
#'   annotations, sliding window identifiers etc. These variables can then be used downstream for statistical comparison of hybrid index 
#'   or clines among groups of loci. In the absence of parental reference samples, certain column names and information must be included 
#'   in the \code{marker.info} file. See examples here and in \code{read.data}.
#'
#' The locus data object includes posterior mean, mode and upper and lower 95 percent credible intervals for frequency of the Source 1 allele
#'   for each locus, assuming beta-distributed proportions. Note that if the locus is fixed in a parental reference set, its posterior 
#'   mode will be 0 or 1 and will fall outside the credible intervals due to the influence of the non-zero prior shape parameters. 
#'   If prior shape parameters are set to zero (equivalent to no prior), the best posterior estimate is the mean rather than the mode, 
#'   but in this case credible intervals cannot be estimated when the locus is fixed for one allele.
#' @return list containing at least two components: (1) \code{data.prep}, a \code{data.table} and \code{data.frame} long-form version of 
#'   the imported data table, with one row per non-missing allele copy and the data columns required for analysis, listed below; (2) 
#'   \code{Nloci.postfilter}, a numeric scalar indicating the number of loci remaining after filtering. If parental reference samples are 
#'   identified, a third output is returned, \code{sourceInfo}, a \code{data.table} and \code{data.frame} listing all the unique entries in 
#'   the declared \sQuote{POPID.name} column (default \sQuote{POPID}) and indicating whether they represent \sQuote{S0}, \sQuote{S1} or 
#'   \sQuote{TEST}. The two optional outputs are: (1) \code{geno.data}, a \code{data.table} and \code{data.frame} containing genotypes for 
#'   each locus scored from \code{0} to \code{PLOIDY}, based on their number of copies of the designated \code{S1} allele (the allele with 
#'   relatively higher frequency in the \code{S1} parental reference set); (2) \code{locus.data}, a \code{data.table} and \code{data.frame} 
#'   with one row per locus and containing useful locus-specific information, including the contents of any uploaded \code{marker.info} 
#'   object.
#'
#'   \code{data.prep} contains the pre-marker columns and the following fields, plus any columns from a declared \code{marker.info} object:
#'   \item{Source}{whether the sample is from parental reference 0 (\sQuote{S0}), parental reference 1 (\sQuote{S1}) or is a test individual 
#'      (\sQuote{TEST}).}
#'   \item{locus}{the locus names.}
#'   \item{S0.prop_1}{frequency (proportion) of the \code{S1_allele} in \code{S0}.}
#'   \item{S1.prop_1}{frequency (proportion) of the \code{S1_allele} in \code{S1}.}
#'   \item{Source_allele}{numeric. Whether the allele is the \code{S0_allele} (\sQuote{0}) or \code{S1_allele} (\sQuote{1}).}
#'
#'   \code{locus.data} additionally contains, at a minimum:
#'   \item{S0_allele}{name of the allele with relatively higher frequency in \code{S0}.}
#'   \item{S1_allele}{name of the allele with relatively higher frequency in \code{S1}.}
#'   \item{Source.afdiff}{difference in allele frequency between \code{S0} and \code{S1}.}
#'   \item{S0.MAF}{\code{S0} minor allele frequency.}
#'   \item{S1.MAF}{\code{S1} minor allele frequency.}
#'   \item{min.MAF}{smaller of \code{S0.MAF} and \code{S1.MAF}.}
#'
#'   If parental reference samples are present, \code{locus.data} further contains:
#'   \item{S0.N_1}{number of non-missing copies of the \code{S1_allele} in \code{S0}.}
#'   \item{S1.N_1}{number of non-missing copies of the \code{S1_allele} in \code{S1}.}
#'   \item{S0.prop_0}{frequency (proportion) of the \code{S0_allele} in \code{S0}.}
#'   \item{S1.prop_0}{frequency (proportion) of the \code{S0_allele} in \code{S1}.}
#'   \item{N_allele_copies.S0}{total number of non-missing allele copies in \code{S0}.}
#'   \item{N_allele_copies.S1}{total number of non-missing allele copies in \code{S1}.}
#'   \item{S0.prop_1_shape1}{shape1 parameter estimate for the posterior beta distribution of the \code{S1_allele} frequency in \code{S0}.}
#'   \item{S0.prop_1_shape2}{shape2 parameter estimate for the posterior beta distribution of the \code{S1_allele} frequency in \code{S0}.}
#'   \item{S0.prop_1_postmode}{mode of the beta-distributed posterior for the frequency of the \code{S1_allele} in \code{S0}.}
#'   \item{S0.prop_1_postmean}{mean of the beta-distributed posterior for the frequency of the \code{S1_allele} in \code{S0}.}
#'   \item{S0.prop_1_cred.int.lower}{lower 95 percent credible interval (2.5 percent quantile of the posterior beta distribution) of the frequency of the \code{S1_allele} in \code{S0}.}
#'   \item{S0.prop_1_cred.int.upper}{upper 95 percent credible interval (97.5 percent quantile of the posterior beta distribution) of the frequency of the \code{S1_allele} in \code{S0}.}
#'   \item{S1.prop_1_shape1}{shape1 parameter estimate for the posterior beta distribution of the \code{S1_allele} frequency in \code{S1}.}
#'   \item{S1.prop_1_shape2}{shape2 parameter estimate for the posterior beta distribution of the \code{S1_allele} frequency in \code{S1}.}
#'   \item{S1.prop_1_postmode}{mode of the beta-distributed posterior for the frequency of the \code{S1_allele} in \code{S1}.}
#'   \item{S1.prop_1_postmean}{mean of the beta-distributed posterior for the frequency of the \code{S1_allele} in \code{S1}.}
#'   \item{S1.prop_1_cred.int.lower}{lower 95 percent credible interval (2.5 percent quantile of the posterior beta distribution) of the frequency of the \code{S1_allele} in \code{S1}.}
#'   \item{S1.prop_1_cred.int.upper}{upper 95 percent credible interval (97.5 percent quantile of the posterior beta distribution) of the frequency of the \code{S1_allele} in \code{S1}.}
#'   \item{S1.prop_0_shape1}{shape1 parameter estimate for the posterior beta distribution of the \code{S0_allele} frequency in \code{S1}.}
#'   \item{S1.prop_0_shape2}{shape2 parameter estimate for the posterior beta distribution of the \code{S0_allele} frequency in \code{S1}.}
#'   \item{S0.prop_0_shape1}{shape1 parameter estimate for the posterior beta distribution of the \code{S0_allele} frequency in \code{S0}.}
#'   \item{S0.prop_0_shape2}{shape2 parameter estimate for the posterior beta distribution of the \code{S0_allele} frequency in \code{S0}.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' ###############################################################
#' #Read in and prepare a file with no parental reference samples#
#' ###############################################################
#' 
#' #Create an input data file and save to the working directory.
#' cat("INDLABEL POPID chr1:001 chr1:002","ind1 pop1 A A","ind1 pop1 A B","ind2 pop1 B B","ind2 pop1 B B","ind3 pop2 NA A","ind3 pop2 B A", file = "ex.data", sep = "\n")
#' 
#' #Create a marker.info file with the necessary parental reference information and save to the working directory.
#' cat("locus refAllele alternateAllele S0.prop_r S1.prop_r type","chr1:001 A B 0.1 0.9 intronic","chr1:002 B A 0.8 0.2 exonic", file = "ex_marker_info2.data", sep = "\n")
#' 
#' #Read in the data and marker.info file.
#' dat <- read.data(file="ex.data",nprecol=2,NUMINDS=3,NUMLOCI=2,MISSINGVAL=NA,marker.info.file = "ex_marker_info2.data", sourceAbsent = TRUE)
#' 
#' #Run data.prep.
#' prepdata=data.prep(
#'  data=dat$data,               #part of the read.data output object#
#'  loci=dat$loci,               #part of the read.data output object#
#'  sourceAbsent = TRUE,         #Must be set to TRUE in the absence of parental reference samples. Default is FALSE#
#'  marker.info=dat$marker.info, #must be specified if sourceAbsent = TRUE; make sure it contains the headers and their 
#'                               #contents as specified above when 'sourceAbsent=TRUE'#
#'  alleles=dat$alleles,         #part of the read.data output object#
#'  #S0,                         #character vector identifying parental reference S0 samples. Leave blank in the absence of parental reference samples#
#'  #S1,                         #character vector identifying parental reference S1 samples. Leave blank in the absence of parental reference samples#
#'  precols=dat$precols,         #part of the read.data output object#
#'  return.genotype.table=TRUE,  #This returns an table, with loci in columns, of diploid (or other ploidy) genotypes, coded as number of copies of the 
#'                               #allele with higher frequency in S1 than S0#
#'  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including data already uploaded plus values 
#'                               #calculated within this function#
#' )
#' 
#' #####################################################################
#' #Read in and prepare a file that includes parental reference samples#
#' #####################################################################
#' 
#' #Make a new data input file with more entries, including S0 and S1 samples, in plink structure format but with missing data recoded to NA. This includes 
#' #a MAPDISTANCE row below the marker names.
#' cat("chr1:001 chr1:002","-3 3 2 4",
#'  "ind1 pop1 A A A B","ind2 pop1 B B B B","ind3 pop2 NA B A A", 
#'  "ind4 pop3 A A A A","ind5 pop3 A B A A","ind6 pop3 A A A A",
#'  "ind7 pop3 NA A A A","ind8 pop4 A NA B B","ind9 pop4 B B B B",
#'  "ind10 pop4 B B B B","ind11 pop5 A A B B",
#'  file = "ex3.data", sep = "\n")
#' 
#' #pop3 represents S0; pop4 and pop5 represent S1.
#' 
#' #Run read.data and data.prep.
#' dat3 <- read.data(file="ex3.data",MISSINGVAL=NA,NUMINDS=11,nprecol=2,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
#' 
#' prepdata= data.prep(
#'  data=dat3$data,
#'  loci=dat3$loci,
#'  sourceAbsent = FALSE,                         #default#
#'  #marker.info=dat$marker.info,                 #a marker.info file is not obligatory when parental reference samples are declared#
#'  alleles=dat3$alleles,
#'  S0="pop3",                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
#'  S1=c("pop4","pop5"),                          #second parental reference set; must be specified when sourceAbsent = FALSE#
#'  precols=dat3$precols,
#'  return.genotype.table=TRUE,                   #default is FALSE to save memory#
#'  return.locus.table=TRUE                       #default is FALSE to save memory#
#' )
#' 
#' #########################################
#' #Now completely exclude the POPID column#
#' #########################################
#' 
#' #Just to show that the functions work in the absence of this column#
#' 
#' #First with no parental reference samples. Modify the plink format 'ex3.data' above to remove the POPID column, 
#' #and use the existing marker info file, 'ex_marker_info2.data'#
#' 
#' cat("chr1:001 chr1:002", "-3 3 2 4",
#'  "ind1 A A A B","ind2 B B B B","ind3 NA B A A", 
#'  "ind4 A A A A","ind5 A B A A","ind6 A A A A",
#'  "ind7 NA A A A","ind8 A NA B B","ind9 B B B B","ind10 B B B B","ind11 A A B B",
#'  file = "ex4.data", sep = "\n")
#' 
#' #you must specify 'POPID=0' (meaning the column is absent) in read.data.
#' dat4 <- read.data(file="ex4.data",POPID=0,nprecol=1,NUMINDS=11,NUMLOCI=2,MISSINGVAL=NA,marker.info.file="ex_marker_info2.data",sourceAbsent=TRUE,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
#' 
#' #Now prepare the data for analysis, this time we're not creating geno.data or locus.data objects.
#' prepdata= data.prep(
#'  data=dat4$data,
#'  loci=dat4$loci,
#'  sourceAbsent = TRUE,         #default is FALSE#
#'  marker.info=dat4$marker.info,#must be specified if sourceAbsent = TRUE#
#'  alleles=dat4$alleles,
#'  #S0,                         #leave blank in the absence of parental reference samples#
#'  #S1,                         #leave blank in the absence of parental reference samples#
#'  POPID.name="INDLABEL",     #***must be specified if there's no POPID column, to avoid an error***#
#'  precols=dat4$precols
#' )
#' 
#' #Now including parental reference samples, again with no POPID column.
#' 
#' #Create and save a marker info file. Not obligatory in the presence of parental reference samples.
#' cat("locus type","chr1:001 intronic","chr1:002 exonic", file = "ex_marker_info.data", sep = "\n")
#' 
#' #Read data in and set sourceAbsent=FALSE (the default setting) so that the code doesn't search for the obligatory columns and throw up an error.
#' dat5 <- read.data(file="ex4.data",POPID=0,nprecol=1,NUMINDS=11,NUMLOCI=2,MISSINGVAL=NA,marker.info.file="ex_marker_info.data",
#'  sourceAbsent=FALSE,ONEROW=1,MAPDISTANCE=1,precol.headers=0)
#' 
#' #Prepare the analysis object.
#' prepdata= data.prep(
#'  data=dat5$data,
#'  loci=dat5$loci,
#'  sourceAbsent = FALSE,                               #default is FALSE#
#'  marker.info=dat5$marker.info,                       #this now contains only the locus (obligatory) and type columns#
#'  alleles=dat5$alleles,
#'  S0=c("ind4","ind5","ind6","ind7"),                  #must be specified when sourceAbsent = FALSE#
#'  S1=c("ind8","ind9","ind10","ind11"),                #must be specified when sourceAbsent = FALSE#
#'  POPID.name = "INDLABEL",                            #***must be specified if there's no POPID column; the S0 and S1 sample 
#'                                                      #references must be present in the column specified here***#
#'  precols=dat5$precols,
#'  return.genotype.table=TRUE,                         #Optional#
#'  return.locus.table=TRUE                             #Optional#
#' )
#' 
#' #################################################################################
#' #Locus filtering based on parental reference allele frequencies and sample sizes#
#' #################################################################################
#' 
#' #Upload a real dataset (Italian sparrows), included as part of the package.
#' dat=read.data("RB_Italy_ggcline_precol_headers_haploidND2.csv",nprecol=2,MISSINGVAL=NA,NUMINDS=569)
#' 
#' #The file contains data for 77 diploid markers for 569 male sparrow individuals.
#' #Run data.prep including the different filtering options. Try hashing out or changing these options to see how it affects the number of 
#' #loci in the resulting analysis table (see prepdata$Nloci.postfilter in the data.prep object).
#' prepdata=data.prep(
#'  data=dat$data,
#'  loci=dat$loci,
#'  alleles=dat$alleles,
#'  S0=c("Kralove","Oslo"), #POPID names for the first parental reference set#
#'  S1="LesinaSPANISH",     #POPID names for the second parental reference set#
#'  precols=dat$precols,
#' ###Filtering below###
#'  max.S.MAF = 0.2,               #minor allele frq must be below this value in at least one of S0 or S1#
#'  min.diff = 0.2,                #loci with parental allele frequency difference < this value will be removed#
#'  min.allele.copies.S0 = 30,     #loci with fewer non-missing allele copies in the S0 parental reference set will be removed#
#'  min.allele.copies.S1 = 12,     #in this dataset the S1 sample size is much smaller than S0, so I'm being less strict with filtering for S1 sample size#
#'  AF.CIoverlap = FALSE,          #***RECOMMENDED*** filtering option if parental reference samples are included - removes all loci for which there is 
#'                                 #overlap between S0 and S1 in the Bayesian posterior 95% credible intervals of allele frequency# 
#' )
#'
#' unlink("ex_marker_info.data")#Tidy up#
#' unlink("ex_marker_info2.data")#Tidy up#
#' unlink("ex.data")#Tidy up#
#' unlink("ex3.data")#Tidy up#
#' unlink("ex4.data")#Tidy up#
#' }
#' @export
data.prep= function(data,loci,alleles,marker.info=NULL,#Setting this as not null automatically means add it to the data object. ***CHANGED FROM BLANK TO 'NULL' - MAKE SURE EVERYTHING STILL WORKS***MAKE SURE THIS IS IN THE HELPFILE. ASSIGN THE OBJECT FROM READ.DATA IF IT'S ALREADY UPLOADED***#
    
	S0=NULL,S1=NULL,sourceAbsent=FALSE,precols,

    INDLABEL.name="INDLABEL",POPID.name="POPID",

    max.S.MAF=0.5,min.diff=0,min.allele.copies.S0=0,min.allele.copies.S1=0,AF.CIoverlap=TRUE,

    return.genotype.table=FALSE,return.locus.table=FALSE,
	
	prior.prop_1=c(1,1)){



    cat(paste("Preparing data for hybrid index and genomic cline estimation..."),

        fill=1); flush.console();

    if(is.null(marker.info)==FALSE){
	    d=copy(marker.info);
	    marker.info=NULL;
	    marker.info=d;
		                           };		

									 
	if(sourceAbsent==TRUE){#Adding source AF to data.prep, plus all locus.data calculations, needs to be done within these brackets in the absence of S0 & S1**********#

        data[,Source:="TEST"];
		
        DT.EXTRACOL=unique(data[,1:length(precols)]);setkeyv(DT.EXTRACOL,INDLABEL.name);
	
	    #20Apr2021 BELOW#FOCUSING ON S0 AND S1 AF BEING IN marker.info*****************
        if(is.null(S1)==FALSE){
	
	        warning("entry for S1 ignored when sourceAbsent==TRUE");
		
	                          };
									 
        if(is.null(S0)==FALSE){
	
	        warning("entry for S0 ignored when sourceAbsent==TRUE");
		
	                          };		
		
		
		if(is.null(marker.info)){
		    stop("S0 and S1 samples must be identified, or their allele frequencies uploaded in a marker.info file using the 'read.data' function.");
		                        };

		    if(sum(names(marker.info)%in%c("S0.prop_r","S1.prop_r")) < 2){
		        stop("columns named 'S0.prop_r' and 'S1.prop_r' must be present and contain allele frequencies of the chosen 
reference allele in source populations S0 and S1 respectively.");
		                                                                 };
		    if(sum(names(marker.info)%in%c("S0.prop_r","S1.prop_r","refAllele","alternateAllele"))!=4){
		        stop("the reference allele (the one whose allele frequency is presented in S0.prop_r and S1.prop_r) and alternate allele 
must both be identified, in columns named 'refAllele' and 'alternateAllele' respectively.");
		                                                                                              };
			
			marker.info[,refAllele:=as.character(refAllele)][,alternateAllele:=as.character(alternateAllele)];

#Assigning reference and alternate alleles to S0 and S1 BELOW#			
			marker.info[S1.prop_r > S0.prop_r,
			    S0.prop_1:=S0.prop_r][S1.prop_r > S0.prop_r,
			    S1.prop_1:=S1.prop_r][S1.prop_r < S0.prop_r,
			    S0.prop_1:=(1 - S0.prop_r)][S1.prop_r < S0.prop_r,
			    S1.prop_1:=(1 - S1.prop_r)];
				
            marker.info[S1.prop_r > S0.prop_r,
                S0_allele:=alternateAllele][S1.prop_r > S0.prop_r,
                S1_allele:=refAllele][S1.prop_r < S0.prop_r,
                S1_allele:=alternateAllele][S1.prop_r < S0.prop_r,
                S0_allele:=refAllele];
										 
		    marker.info[,S1.prop_r:=NULL][,S0.prop_r:=NULL][,refAllele:=NULL][,alternateAllele:=NULL]#Removing columns no longer needed#					 
#Assigning reference and alternate alleles to S0 and S1 ABOVE#											 

			#1. Add details needed for filtering.#
			
            marker.info[,Source.afdiff:= S1.prop_1 - S0.prop_1];

            setkey(marker.info,S0.prop_1,S1.prop_1);

            marker.info[,S0.MAF:= 0.5 - abs(0.5 - S0.prop_1)];

            marker.info[,S1.MAF:= 0.5 - abs(0.5 - S1.prop_1)];

            setkey(marker.info,S0.MAF,S1.MAF);

            marker.info[S0.MAF<S1.MAF,min.MAF:=S0.MAF];

            marker.info[S0.MAF>=S1.MAF,min.MAF:=S1.MAF];
			
			#2. Filter. For AF.CIoverlap, the user only has to supply 'S0.prop_1_cred.int.upper' and 'S1.prop_1_cred.int.lower'.#
			
			if(AF.CIoverlap==FALSE){
			
		        if(sum(names(marker.info)%in%c("S0.prop_1_cred.int.upper","S1.prop_1_cred.int.lower")) < 2){
		            stop("if 'AF.CIoverlap==FALSE' columns named 'S0.prop_1_cred.int.upper' and 'S1.prop_1_cred.int.lower' must be present 
in the marker.info object and contain user-supplied credible/confidence limits for the frequency of one allele (the one with higher frequency 
in source S1). Only the upper limit for source S0 and the lower limit for source S1 are required for locus 
filtering. Alternatively, include parental reference samples (S0 and S1) in the dataset. In the presence of parental 
reference samples, the 95% credible intervals of parental allele frequencies are calculated automatically and these 
are used for filtering.");						
						                                                                                   };
	
	            marker.info=marker.info[S0.prop_1_cred.int.upper<S1.prop_1_cred.int.lower];#These column names must be present if this is to be used#
	
	                               };
			
            marker.info=marker.info[min.MAF <= max.S.MAF];

            setkey(marker.info,S0.prop_1,S1.prop_1);

            marker.info=marker.info[S0.prop_1 != "NA"][(S1.prop_1 - S0.prop_1) >= min.diff];

#min.allele.copies filtering#
            if(sum(min.allele.copies.S0,min.allele.copies.S1) > 0){
		        if(sum(names(marker.info)%in%c("N_allele_copies.S0","N_allele_copies.S1")) < 2){
		            stop("if 'sourceAbsent=TRUE' and you wish to filter by the number of allele copies in the parental reference 
samples, these numbers must be user-supplied in the marker.info file in columns named 'N_allele_copies.S0' 
and 'N_allele_copies.S1'.");						
						                                                                       };
																							   
                marker.info = marker.info[N_allele_copies.S0 >= min.allele.copies.S0][N_allele_copies.S1 >= 
                    min.allele.copies.S1];																																											   
			                                                       };#End of if(sum(min.allele.copies.S0,min.allele.copies.S1) > 0)

#############3. Join to the data.#First need to create the data.prep object###########################################################
			
        keepNloc=as.integer(2.5e+07/nrow(data));

        Numi=ceiling(length(loci)/keepNloc);

        for(i in 1:Numi){

            if(i==Numi){

                DT.2.0 = melt(data, id = c("Source",INDLABEL.name), 
                    measure = (length(precols)+keepNloc*(i-1)+1):(ncol(data)-1),
                    variable.name = "locus", 
                    value.name = "allele");

                       }else{

                DT.2.0 = melt(data, id = c("Source",INDLABEL.name), 
                    measure = (length(precols)+keepNloc*(i-1)+1):(length(precols)+keepNloc*i),
                    variable.name = "locus", 
                    value.name = "allele");

                            };

            if(i==1){DT.2=DT.2.0}else{DT.2=rbind(DT.2,DT.2.0)};

                         };
			
		
#Joining any extra marker.info columns, beyond the bare minimum for sourceAbsent=TRUE, to the data analysis table##########
		
        #namesi=c("S0.prop_1","S1.prop_1","Source.afdiff","S0.MAF","S1.MAF","min.MAF");
		namesi=c("Source.afdiff","S0.MAF","S1.MAF","min.MAF");

        #if(ncol(marker.info[,-..namesi]) > 1){

            setkey(DT.2,locus);setkey(marker.info,locus);
            DT.2=DT.2[marker.info[,-..namesi]];
                                             #};
		
		DT.2=DT.EXTRACOL[DT.2,on=INDLABEL.name];
		
		
#Remove missing data#########

        DT.2 = DT.2[is.na(allele)==F];
		
#Add Source_allele column####

        setkey(DT.2,allele,S0_allele,S1_allele);
        DT.2[allele==S0_allele,Source_allele := 0];
        DT.2[allele==S1_allele,Source_allele := 1];
		
		DT.2=DT.2[,-c("S0_allele","S1_allele")];
		
		
		
#Add genotype table##########

        if(return.genotype.table==TRUE){#***USE THIS SAME CODE TO GET GENOTYPES IN THE PLOT CLINECURVE FUNCTION***#
            geno=dcast(DT.2[,c(precols,"locus","Source","Source_allele"),with=F],
                ...~locus,value.var="Source_allele",fun=sum);

            geno[is.na(geno)]=NA;
                                       };


#Add locus.data object#######

        if(return.locus.table==TRUE){

            locus.data=marker.info;
		                            };
									
	    DT.2[,index:=seq(1,.N)];#New 02 March 2022#

			
        output=list();

            output$data.prep = DT.2[,-c("allele")]#[,-c("Source","allele")];#Keeping Source column to avoid errors in downstream functions#
		
		    output$Nloci.postfilter = marker.info[,.N];	
			
        if(return.genotype.table==TRUE){

            output$geno.data = geno;

                                       };


        if(return.locus.table==TRUE){

            output$locus.data = locus.data;

                                    };
			
			
		
        cat(paste("Done."),fill=1);

        return(output);
			

			
			

													
#################################*************************************
#End of 'if(sourceAbsent==TRUE)'#*************************************
#################################*************************************			
	                                 }else{#Adding source AF to data.prep, plus all locus.data calculations, needs to be done within the previous brackets in the absence of S0 & S1**********#

    setkeyv(data,POPID.name);

    data[,Source:="TEST"];

    data[S0,Source:="S0"];

    data[S1,Source:="S1"];
	
    if(data[S0,.N]==0&data[S1,.N]!=0){
	
	    stop("either identify both S0 and S1, or neither. If neither, set 'sourceAbsent==TRUE'.");
		
	                                 };
									 
    if(data[S0,.N]!=0&data[S1,.N]==0){
	
	    stop("either identify both S0 and S1, or neither If neither, set 'sourceAbsent==TRUE'.");
		
	                                 };
									
									
    data[,loci] = data[,lapply(.SD,replace.fun.A),.SDcols=loci[]];

    data[,loci] = data[,lapply(.SD,replace.fun),.SDcols=loci[]];



    DT.EXTRACOL=unique(data[,1:length(precols)]);setkeyv(DT.EXTRACOL,INDLABEL.name);

	

    setkey(data,Source);

#******ALL THE PREPARATIONS OF THE PARENTAL REFERENCE TABLES ARE BELOW*****************************#

    keepNloc=as.integer(2.5e+07/nrow(data["S0"]));#OKAY SO THIS IS DIVIDING UP THE TABLE SO AS TO REDUCE THE CHANCE OF RUNNING OUT OF MEMORY WITH dcast(melt())#

    Numi=ceiling(length(loci)/keepNloc);



    for(i in 1:Numi){



        if(i==Numi){

            DT.S0.0 = dcast(

                melt(data["S0"], id = INDLABEL.name, measure = (length(precols)+keepNloc*(i-1)+1):(ncol(data)-1),

                variable.name = "locus", na.rm = T)[, 

                .N, by = .(locus, value)][, "prop" := N/sum(N), by = .(locus)],

                locus ~ value, value.var = c("N","prop"),fill = 0);

                   }else{

            DT.S0.0 = dcast(

                melt(data["S0"], id = INDLABEL.name, measure = (length(precols)+keepNloc*(i-1)+1):(length(precols)+keepNloc*i),

                variable.name = "locus", na.rm = T)[, 

                .N, by = .(locus, value)][, "prop" := N/sum(N), by = .(locus)],

                locus ~ value, value.var = c("N","prop"),fill = 0);

                        };



        if(i==1){DT.S0=DT.S0.0}else{DT.S0=rbind(DT.S0,DT.S0.0)};



                    };



    keepNloc=as.integer(2.5e+07/nrow(data["S1"]));

    Numi=ceiling(length(loci)/keepNloc);



    for(i in 1:Numi){



        if(i==Numi){

            DT.S1.0 = dcast(

                melt(data["S1"], id = INDLABEL.name, measure = (length(precols)+keepNloc*(i-1)+1):(ncol(data)-1),

                variable.name = "locus", na.rm = T)[, 

                .N, by = .(locus, value)][, "prop" := N/sum(N), by = .(locus)],

                locus ~ value, value.var = c("N","prop"),fill = 0);

                   }else{

            DT.S1.0 = dcast(

                melt(data["S1"], id = INDLABEL.name, measure = (length(precols)+keepNloc*(i-1)+1):(length(precols)+keepNloc*i),

                variable.name = "locus", na.rm = T)[, 

                .N, by = .(locus, value)][, "prop" := N/sum(N), by = .(locus)],

                locus ~ value, value.var = c("N","prop"),fill = 0);

                        };



        if(i==1){DT.S1=DT.S1.0}else{DT.S1=rbind(DT.S1,DT.S1.0)};



                    }



    setkey(DT.S0,locus);setkey(DT.S1,locus);

    DT.S = DT.S0[DT.S1];



    DT.S[,propdiff_A := i.prop_A - prop_A][,propdiff_B := i.prop_B - prop_B];



    setkey(DT.S,propdiff_A);

    DT.S[propdiff_A < 0, S0_allele := "A"];#****THIS NEEDS TO BE DONE ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#

    setkey(DT.S,propdiff_B);

    DT.S[propdiff_B < 0, S0_allele := "B"];#****THIS NEEDS TO BE DONE ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#



    setkey(DT.S,propdiff_A);

    DT.S[propdiff_A > 0, S1_allele := "A"][propdiff_A > 0, S0.prop_1 := prop_A][propdiff_A > 0, #****THIS NEEDS TO BE DONE ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#

        S1.prop_1 := i.prop_A][propdiff_A > 0, S0.N_1 := N_A][propdiff_A > 0, 

        S1.N_1 := i.N_A];



    setkey(DT.S,propdiff_B);

    DT.S[propdiff_B > 0, S1_allele := "B"][propdiff_B > 0, S0.prop_1 := prop_B][propdiff_B > 0, #****THIS NEEDS TO BE DONE ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#

        S1.prop_1 := i.prop_B][propdiff_B > 0, S0.N_1 := N_B][propdiff_B > 0, 

        S1.N_1 := i.N_B];



    DT.S[,S0.prop_0:= (1-S0.prop_1)][,S1.prop_0:= (1-S1.prop_1)];#***SHOULD I REMOVE THIS?***#



    DT.S[,N_allele_copies.S0 := N_A + N_B];

    DT.S[,N_allele_copies.S1 := i.N_A + i.N_B];

    DT.S[,Source.afdiff:= S1.prop_1 - S0.prop_1];



    setkey(DT.S,S0.prop_1,S1.prop_1);

    DT.S[,S0.MAF:= 0.5 - abs(0.5 - S0.prop_1)];

    DT.S[,S1.MAF:= 0.5 - abs(0.5 - S1.prop_1)];



    setkey(DT.S,S0.MAF,S1.MAF);

    DT.S[S0.MAF<S1.MAF,min.MAF:=S0.MAF];

    DT.S[S0.MAF>=S1.MAF,min.MAF:=S1.MAF];



    DT.S = DT.S[min.MAF <= max.S.MAF];#***FILTERING. I MAY NEED TO INCLUDE THIS ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE INDIVIDUALS***#



    setkey(DT.S,S0.prop_1,S1.prop_1);

    DT.S = DT.S[S0.prop_1 != "NA"][(S1.prop_1 - S0.prop_1) >= min.diff];#***FILTERING. I MAY NEED TO INCLUDE THIS ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE INDIVIDUALS***#

    setkey(DT.S,N_allele_copies.S0,N_allele_copies.S1);

    DT.S = DT.S[N_allele_copies.S0 >= min.allele.copies.S0][N_allele_copies.S1 >= 

        min.allele.copies.S1];#***FILTERING. I MAY NEED TO INCLUDE THIS ELSEWHERE IN THE ABSENCE OF PARENTAL REFERENCE INDIVIDUALS***#


#******ALL THE PREPARATIONS OF THE PARENTAL REFERENCE TABLES ARE ABOVE*****************************#

    keepNloc=as.integer(2.5e+07/nrow(data));



    Numi=ceiling(length(loci)/keepNloc);



    for(i in 1:Numi){



        if(i==Numi){

            DT.2.0 = melt(data, id = c("Source",INDLABEL.name), 

                measure = (length(precols)+keepNloc*(i-1)+1):(ncol(data)-1),

                variable.name = "locus", 

                value.name = "allele");

                   }else{

            DT.2.0 = melt(data, id = c("Source",INDLABEL.name), 

                measure = (length(precols)+keepNloc*(i-1)+1):(length(precols)+keepNloc*i),

                variable.name = "locus", 

                value.name = "allele");

                        };



        if(i==1){DT.2=DT.2.0}else{DT.2=rbind(DT.2,DT.2.0)};

                     };



    setkey(DT.2,locus); setkey(DT.S,locus);#****DON'T DO THIS IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#

    DT.join = DT.2[DT.S];#****DON'T DO THIS IN THE ABSENCE OF PARENTAL REFERENCE SAMPLES****#



    setkey(DT.join,allele,S0_allele,S1_allele);

    DT.join[allele==S0_allele,Source_allele := 0];

    DT.join[allele==S1_allele,Source_allele := 1];



    alleles = data.table(t(alleles),keep.rownames=T);

    setnames(alleles,1:3,c("locus","allele_1","allele_2"));



    setkey(DT.join,locus);setkey(alleles,locus);

    DT.join = alleles[DT.join];

    setkey(DT.join,allele,S0_allele,S1_allele);

    DT.join[S0_allele=="A",S0_allele := allele_1];

    DT.join[S0_allele=="B",S0_allele := allele_2];

    DT.join[S1_allele=="A",S1_allele := allele_1];

    DT.join[S1_allele=="B",S1_allele := allele_2];

    DT.join[allele=="A",allele := allele_1];

    DT.join[allele=="B",allele := allele_2];



    setnames(DT.join,old = c("N_A","N_B","prop_A","prop_B","i.N_A","i.N_B",

        "i.prop_A","i.prop_B"),

        new = c("N_A.S0","N_B.S0","prop_A.S0","prop_B.S0","N_A.S1","N_B.S1",

        "prop_A.S1","prop_B.S1"));



    DT.join[,allele_1:=NULL][,allele_2:=NULL][,N_A.S0:=NULL][,N_B.S0:=NULL][,

        prop_A.S0:=NULL][,prop_B.S0:=NULL][,N_A.S1:=NULL][,N_B.S1:=NULL][,

        prop_A.S1:=NULL][,prop_B.S1:=NULL][,propdiff_A:=NULL][,propdiff_B:=NULL];



    DT.join=DT.EXTRACOL[DT.join,on=INDLABEL.name];
	
	
	
	
#****	
	locus.dat=unique(DT.join[,.(locus,S0.N_1,S1.N_1,N_allele_copies.S0,N_allele_copies.S1)]);
			
	
	locus.dat[,S0.prop_1_shape1 := prior.prop_1[1] + S0.N_1][,
	
        S0.prop_1_shape2 := prior.prop_1[2] + (N_allele_copies.S0 - S0.N_1)][,
		
        S0.prop_1_postmode := (S0.prop_1_shape1 - 1)/(S0.prop_1_shape1 + S0.prop_1_shape2 - 2)][,
		
        S0.prop_1_postmean := 1/(1 + S0.prop_1_shape2/S0.prop_1_shape1)];
		
		

    setkey(locus.dat,S0.prop_1_shape1,S0.prop_1_shape2);
	
    locus.dat[S0.prop_1_shape1<1&S0.prop_1_shape2>=1,
	
        S0.prop_1_postmode := 0][S0.prop_1_shape1>=1&S0.prop_1_shape2<1,S0.prop_1_postmode := 1][,
		
        S0.prop_1_cred.int.lower := qbeta(0.025,shape1 = S0.prop_1_shape1,shape2 = S0.prop_1_shape2)][,
		
        S0.prop_1_cred.int.upper := qbeta(0.975,shape1 = S0.prop_1_shape1,shape2 = S0.prop_1_shape2)];
		
		

    locus.dat[,S1.prop_1_shape1 := prior.prop_1[1] + S1.N_1][,
	
        S1.prop_1_shape2 := prior.prop_1[2] + (N_allele_copies.S1 - S1.N_1)][,
		
        S1.prop_1_postmode := (S1.prop_1_shape1 - 1)/(S1.prop_1_shape1 + S1.prop_1_shape2 - 2)][,
		
        S1.prop_1_postmean := 1/(1 + S1.prop_1_shape2/S1.prop_1_shape1)];
		
		

    setkey(locus.dat,S1.prop_1_shape1,S1.prop_1_shape2);
	
    locus.dat[S1.prop_1_shape1<1&S1.prop_1_shape2>=1,
	
        S1.prop_1_postmode := 0][S1.prop_1_shape1>=1&S1.prop_1_shape2<1,S1.prop_1_postmode := 1][,
		
        S1.prop_1_cred.int.lower := qbeta(0.025,shape1 = S1.prop_1_shape1,shape2 = S1.prop_1_shape2)][,
		
        S1.prop_1_cred.int.upper := qbeta(0.975,shape1 = S1.prop_1_shape1,shape2 = S1.prop_1_shape2)];
		
	
	
	locus.dat[,S1.prop_0_shape1 := S1.prop_1_shape2][,
	
	    S1.prop_0_shape2 := S1.prop_1_shape1][,
		
		S0.prop_0_shape1 := S0.prop_1_shape2][,
		
		S0.prop_0_shape2 := S0.prop_1_shape1];

		
		
	setkey(locus.dat,locus);setkey(DT.join,locus);
	
	DT.join=DT.join[locus.dat[,!c("S0.N_1","S1.N_1","N_allele_copies.S0","N_allele_copies.S1")]];	 
#****

    if(AF.CIoverlap==FALSE){
	
	    DT.join=DT.join[S0.prop_1_cred.int.upper<S1.prop_1_cred.int.lower];
	
	                       };

#****



    if(return.genotype.table==TRUE){#***USE THIS SAME CODE TO GET GENOTYPES IN THE PLOT CLINECURVE FUNCTION***#

        geno=dcast(DT.join[,c(precols,"locus","Source","Source_allele"),with=F],

            ...~locus,value.var="Source_allele",fun=sum);

        geno[is.na(geno)]=NA;

                                   };



    if(return.locus.table==TRUE){

        locus.data=unique(DT.join[,.(locus,S0_allele,S1_allele,S0.prop_1,S1.prop_1,S0.N_1,S1.N_1,S0.prop_0,
		
            S1.prop_0,N_allele_copies.S0,N_allele_copies.S1,Source.afdiff,S0.MAF,S1.MAF,min.MAF,S0.prop_1_shape1,
			
			S0.prop_1_shape2,S0.prop_1_postmode,S0.prop_1_postmean,S0.prop_1_cred.int.lower,S0.prop_1_cred.int.upper,
			
			S1.prop_1_shape1,S1.prop_1_shape2,S1.prop_1_postmode,S1.prop_1_postmean,S1.prop_1_cred.int.lower,
			
			S1.prop_1_cred.int.upper,S1.prop_0_shape1,S1.prop_0_shape2,S0.prop_0_shape1,S0.prop_0_shape2)]);



        if(is.null(marker.info)==F){

            setkey(locus.data,locus);setkey(marker.info,locus);

            locus.data=locus.data[marker.info];

                                        };



                                };
								
								
#Making a sourceInfo object#

    namesj=c(POPID.name,"Source");
    sourceInfo=unique(DT.join[,..namesj]);
    setkey(sourceInfo,Source);


    namesi=c(precols,"locus","Source","S0.prop_1","S1.prop_1","Source_allele");
    DT.join=DT.join[,..namesi];
			

    if(is.null(marker.info)==F){

        setkey(DT.join,locus);setkey(marker.info,locus);

        DT.join=DT.join[marker.info];

                               };

    DT.join=DT.join[is.na(Source_allele)==F];
	
	DT.join[,index:=seq(1,.N)];#New 16 Nov 2021#
	

    output=list();

        output$data.prep = DT.join;#[,-"Source"];#I need this column to avoid errors in downstream functions#
		
		output$sourceInfo = sourceInfo;
		
		output$Nloci.postfilter = length(DT.join[,unique(locus)]);



    if(return.genotype.table==TRUE){

        output$geno.data = geno;

                                   };



    if(return.locus.table==TRUE){

        output$locus.data = locus.data;

                            };



    cat(paste("Done."),fill=1);



    return(output)

                                                          };
														  
														  }
