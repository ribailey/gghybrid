#' Load data stored in the format of a \sQuote{structure} file or similar data table.
#'
#' @param file Character string. The name of the data file to read in.
#' @param mainparams Character string. The name of an associated \code{mainparams} 
#'   file to read in. Optional as \code{mainparams} info can be entered into the 
#'   function call directly (see below), or may be unnecessary (see \sQuote{Details}).
#' @param extracol.names Character string or vector. Names of extra (i.e. none of those specifically named here) 
#'   non-marker columns in the data file. Optional; used 
#'   when the data file does not contain pre-marker \sQuote{extracol} headers and you 
#'   wish to add them on input. Default is \code{NULL}.
#' @param precol.headers Numeric. Presence (\code{1}) or absence (\code{0}) of headers for 
#'   pre-marker columns in the header row of the data file. Default is \code{1}.
#' @param nprecol Numeric. Number of pre-marker columns in the data file. No 
#'   default and must always be entered.
#' @param markername.dup Numeric. Whether each marker name appears twice in 
#'   the data file header row. Default is \code{0} (FALSE), alternative is \code{1} (TRUE).
#' @param EXTRACOL Numeric. Presence (\code{1}) or absence (\code{0}) of extra (i.e. none of those 
#'   specifically named here) non-marker columns in the data file. Not needed if this info is 
#'   uploaded in a \code{mainparams} file. For compatibility with structure \code{mainparams} info. Default is \code{0}.
#' @param INDLABEL Numeric. Presence (\code{1}) or absence (\code{0}) of a column of 
#'   individual references in the data file. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{1}, and this column 
#'   must always be present. See \sQuote{Details}.
#' @param LOCDATA Numeric. Presence (\code{1}) or absence (\code{0}) of a \code{LOCDATA} column in 
#'   the data file. For compatibility with sQuote{structure}. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{0}.
#' @param MAPDISTANCE Numeric. Presence (\code{1}) or absence (\code{0}) of a \code{MAPDISTANCE} 
#'   row in the data file, which will be removed if present. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{0}.
#' @param MARKERNAME Numeric. Presence (\code{1}) or absence (\code{0}) of a header row 
#'   containing marker names in the data file. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{1}.
#' @param MISSINGVAL The identifier for missing values. No default. Not needed if the info is uploaded in a \code{mainparams} file.
#' @param NUMINDS Numeric. Number of individuals in the data set. Not used 
#'   in downstream functions and included for user's convenience. Default 
#'   set to \code{0}.
#' @param NUMLOCI Numeric. Number of loci in the data set. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{0}, but 
#'   the correct number must be entered for most input file formats (see 
#'   \sQuote{Details}).
#' @param ONEROW Numeric. Whether (\code{1}) or not (\code{0}) there is a single data row 
#'   per \code{INDLABEL}. If the value is \code{1}, the number of columns per marker = 
#'   \code{PLOIDY}. If \code{ONEROW = 0}, there is one column per marker and the number 
#'   of rows per \code{INDLABEL} = \code{PLOIDY}. Not needed if the info is uploaded in a \code{mainparams} file. Default is \code{0}.
#' @param PHASED Numeric. Presence (\code{1}) or absence (\code{0}) of phasing data rows 
#'   in a \sQuote{structure} input file. These rows will be removed if present. Not needed if the info is uploaded in a \code{mainparams} file. 
#'   Default is \code{0}.
#' @param PHENOTYPE Numeric. Presence (\code{1}) or absence (\code{0}) of a \code{PHENOTYPE} 
#'   column in the data file. For compatibility with \sQuote{structure}. Not needed if the info is uploaded in a \code{mainparams} file. Default 
#'   is \code{0}.
#' @param PLOIDY Numeric. Maximum ploidy among the markers in the data file. Not needed if the info is uploaded in a \code{mainparams} file. 
#'   Default is \code{2}.
#' @param POPFLAG Numeric. \code{0} or \code{1}, for compatibility with structure. Not 
#'   needed and will be removed in future. Default \code{0}.
#' @param POPID Numeric. Presence (\code{1}) or absence (\code{0}) of a column of population 
#'   identifiers in the data file (which should always be present, see \sQuote{Details}). Not needed if the info is uploaded in a \code{mainparams} file. 
#'   Default is \code{1}.
#' @param RECESSIVEALLELE Numeric. \code{0} or \code{1}, for compatibility with \sQuote{structure}. Not 
#'   needed and will be removed in future. Default \code{0}.
#' @details \code{read.data} is designed to be compatible with data input files 
#'   for the software \href{http://web.stanford.edu/group/pritchardlab/structure.html}{structure}, including those produced by 
#'   \href{https://www.cog-genomics.org/plink2/}{PLINK}. However, the simplest file format to read in is
#'   a data table with a complete header row, one column per marker and one row per allele copy (equivalent to \code{ONEROW
#'   = 0}). With this format the only required fields
#'   are \code{file}, \code{nprecol} and \code{MISSINGVAL}. Furthermore, allele identities are 
#'   read in as character strings and therefore each allele can be given any 
#'   name or number, as long as there are only two unique alleles per locus. 
#'   As with \href{http://web.stanford.edu/group/pritchardlab/structure.html}{structure}, all non-marker 
#'   columns should be to the left of the marker columns in the data file.
#'
#' The data input file is required to have at least two columns prior to the 
#'   marker data columns. The first is referred to as \sQuote{INDLABEL} (but can be 
#'   given a different column header) and represents the finest resolution of 
#'   identification that may be required for the data set (usually a unique 
#'   individual reference, but could for example be a population-level 
#'   reference for poolseq data). The second is referred to as \sQuote{POPID} and 
#'   indicates the population from which the individual was sampled. The main 
#'   purpose of this column is to declare, downstream, which \code{POPID} values 
#'   represent the parental reference samples, and it can also be useful for 
#'   plotting hybrid index estimates.
#'
#' Any ploidy is allowed, but the declared \code{PLOIDY} must be the same for all markers. So for example if 
#'   data are haplodiploid, haploid markers should either be present as 
#'   diploid homozygotes, or (preferably from a statistical viewpoint) the 
#'   second allele for each \code{INDLABEL} and marker should be declared as missing data.
#'
#' If an associated \code{mainparams} file is read in, it should contain two columns 
#'   and no header row, with the first column holding the field names and the 
#'   second the field values. The only mandatory fields are  
#'   (and must have one of these synonyms): \code{NUMLOCI}, \code{PLOIDY}, 
#'   \code{MISSINGVAL}/\code{MISSING}, \code{ONEROW}/\code{ONEROWPERIND}, \code{INDLABEL}, \code{POPID}, 
#'   \code{EXTRACOL}/\code{EXTRACOLS}, \code{MARKERNAME}/\code{MARKERNAMES}. If \code{RECESSIVEALLELE} or 
#'   \code{MAPDISTANCE} are present, they can be pluralized.
#'
#' \code{read.data} uses the \code{\link[data.table]{fread}} function from \pkg{\href{https://cran.r-project.org/web/packages/data.table/index.html}{data.table}} 
#'   to rapidly read in the data file. The 
#'   loaded data file is therefore of \href{http://127.0.0.1:23982/library/base/html/class.html}{class} \code{data.table} and \code{data.frame}. 
#'   A \code{data.table} has slightly different syntax to a standard \code{data.frame}
#'   for data manipulation. Hence if you wish to manipulate the imported data file in any way, such as subsetting, and are not comfortable using 
#'   \code{data.table} syntax, convert the file to a standard \code{data.frame} first. However, it must be converted back to a 
#'   \code{data.table} before use in downstream functions.
#' @return \code{read.data} returns an object of \href{http://127.0.0.1:23982/library/base/html/class.html}{class} \dQuote{XXXX} 
#'   containing the loaded data in \code{ONEROW = 0} format, in the form of a 
#'   \code{data.table} and \code{data.frame}, and other components used in downstream functions or otherwise potentially useful to the user.
#'
#'   An object of \href{http://127.0.0.1:23982/library/base/html/class.html}{class} \dQuote{XXXX} is a list containing the following components:
#'   \item{mainparams}{A \code{data.table} and \code{data.frame} with the inputted \code{mainparams} information (or default values if 
#'   no information provided).}
#'   \item{nprecols}{Numeric. The number of non-marker columns to the left of the marker columns in the imported data set.}
#'   \item{precols}{Character vector. The names of non-marker columns to the left of the marker columns in the imported data set.}
#'   \item{data}{A \code{data.table} and \code{data.frame}. The imported data.}
#'   \item{alleles}{A \code{data.table} and \code{data.frame}. The names of the two alleles at each locus.}
#'   \item{loci}{Character vector. The names of all loci in the imported data set.}
#' @author
#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}
#' @examples
#'
#' \dontrun{
#' #For a table with complete header row, \code{PLOIDY=2} (default) and one column per locus:
#' dat=read.data("filename.ext",
#' nprecol=2,
#' MISSINGVAL=NA)
#'
#' #For a structure-format file with only locus names on the header row and \code{PLOIDY} columns per locus:
#' dat2=read.data(file="filename.ext",
#' precol.headers=0,
#' nprecol=2,
#' EXTRACOL=0,		#Same as default; included for clarity#
#' INDLABEL=1,		#Same as default; included for clarity#
#' LOCDATA=0,		#Same as default; included for clarity#
#' MAPDISTANCE=0,	#Same as default; included for clarity#
#' MARKERNAME=1,	#Same as default; included for clarity#
#' MISSINGVAL=-9,
#' NUMINDS=354,		#Not mandatory#
#' NUMLOCI=77,
#' ONEROW=1,
#' PHASED=0,		#Same as default; included for clarity#
#' PHENOTYPE=0,		#Same as default; included for clarity#
#' PLOIDY=2,		#Same as default; included for clarity#
#' POPFLAG=0,		#Same as default; included for clarity#
#' POPID=1,			#Same as default; included for clarity#
#' RECESSIVEALLELE=0#Same as default; included for clarity#
#' )
#' }
#' @export
read.data=function(file,mainparams=NULL,extracol.names=NULL,
    precol.headers=1,nprecol,markername.dup=0,
    EXTRACOL=0,INDLABEL=1,LOCDATA=0,MAPDISTANCE=0,MARKERNAME=1,MISSINGVAL,
    NUMINDS=0,NUMLOCI=0,ONEROW=0,PHASED=0,PHENOTYPE=0,PLOIDY=2,POPFLAG=0,
    POPID=1,RECESSIVEALLELE=0){
 
    cat(paste("Reading in data file..."),fill=1); flush.console();
	 
	if(is.null(mainparams)==F){

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

    setkey(mainparams,V1);
    ntop=mainparams[c("MARKERNAME","RECESSIVEALLELE","MAPDISTANCE"),
	    sum(as.numeric(V2))];
    nprecolumns=mainparams[c("INDLABEL","POPID",
	    "POPFLAG","LOCDATA","PHENOTYPE","EXTRACOL"),sum(as.numeric(V2))];
    precols=subset(mainparams[c("INDLABEL","POPID",
	    "POPFLAG","LOCDATA","PHENOTYPE")],V2=="1");
	precols=precols[,as.character(V1)];

    if(sum(mainparams[c("MAPDISTANCE","MARKERNAME","markername.dup",
        "ONEROW","precol.headers","RECESSIVEALLELE"),
        V2==c("0","1","0","0","1","0")]) == 6
       ){
        data=fread(file,header=TRUE,
		    na.strings=paste(mainparams["MISSINGVAL",V2]),
            colClasses="character");

        MARKER=names(data[,(nprecol+1):ncol(data),with=F]);

        if(mainparams["PHASED",V2=="1"]){
		    data = data[rep(c(rep(TRUE,mainparams["PLOIDY",as.numeric(V2)]),
			    FALSE),length=.N)];
			                            };

        allele.id=data[,lapply(.SD,unique.fun),.SDcols=MARKER[]];

        output=list();
            output$mainparams=mainparams;
            output$nprecols=nprecol;
            output$precols=names(data[,1:nprecol])
            output$data=data;
            output$alleles=allele.id;
            output$loci=MARKER;

        cat(paste("Done"),fill=1);

        return(output)
		
         }else{
        data=fread(file,header=FALSE,skip=ntop,
		    na.strings=paste(mainparams["MISSINGVAL",V2]),
			colClasses="character");

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
			    nprecolumns),what=character(),quiet=TRUE);	
		    MARKER=cnames[(nprecolumns + 1):length(cnames)];
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
				nprecolumns),what=character(),quiet=TRUE));	
            MARKER=cnames[(nprecolumns + 1):length(cnames)];
			                      };

        if(mainparams["precol.headers",V2=="0"]){
            if(is.null(extracol.names)){
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
				                            };

            allele.id=data[,lapply(.SD,unique.fun),.SDcols=MARKER[]];

            output<-list();
                output$mainparams=mainparams;
                output$nprecols=nprecolumns;
                output$precols=cnames[1:nprecolumns]
                output$data=data;
                output$alleles=allele.id;
                output$loci=MARKER;

            cat(paste("Done"),fill=1);

            return(output)									
									    }else{
            if(mainparams["PHASED",V2=="1"]){
			    data = data[rep(c(TRUE,FALSE),length=.N)];
				                            };

            colNum=c((1:nprecolumns),
                seq(
				    (nprecolumns + 1),
                    (mainparams["NUMLOCI",as.numeric(V2)]*mainparams["PLOIDY",
					as.numeric(V2)] + nprecolumns - (mainparams["PLOIDY",
					as.numeric(V2)]) + 1),by=mainparams["PLOIDY",
					as.numeric(V2)]
				   )
			        );

            dat=data[0,colNum,with=F];setnames(dat,cnames);setkey(dat);
						
            for(i in 1:mainparams["PLOIDY",as.numeric(V2)]){
                colNum=seq(
				    (nprecolumns + i),(mainparams["NUMLOCI",
					as.numeric(V2)]*mainparams["PLOIDY",
					as.numeric(V2)] + nprecolumns - (mainparams["PLOIDY",
					as.numeric(V2)]) + i),
                    by=mainparams["PLOIDY",as.numeric(V2)]
					      );
                dat1=data[, .SD, .SDcols=colNum];
                dat1=cbind(data[,1:nprecolumns],dat1);
                setnames(dat1,cnames);setkey(dat1);
                dat=rbind(dat,dat1);
                                                           };
												   
            data=NULL;
			
            allele.id=dat[,lapply(.SD,unique.fun),.SDcols=MARKER[]];
			
            output=list();
                output$mainparams=mainparams;
                output$nprecols=nprecolumns;
                output$precols=cnames[1:nprecolumns]
                output$data=dat;
                output$alleles=allele.id;
                output$loci=MARKER;

            cat(paste("Done"),fill=1);

            return(output)		
	                                         };
               };
                     }