#' Prepare data for hybrid index and genomic cline estimation.

#'

#' @param data A \code{data} object produced by \code{read.data}, or a custom \code{data.table} in the same format.

#' @param loci A \code{loci} object produced by \code{read.data}, or a custom character vector with locus names.

#' @param alleles An \code{alleles} object produced by \code{read.data}, or a custom \code{data.table} in the same format.

#' @param S0 Character vector. Names, taken from the designated \code{POPID} column in the \code{data} object, 

#'   indicating which populations represent the first parental reference set (\dQuote{Source 0}). Allele frequencies 

#'   in this set represent a hybrid index of 0. 

#' @param S1 Character vector. Names, taken from the designated \code{POPID} column in the \code{data} object, 

#'   indicating which populations represent the second parental reference set (\dQuote{Source 1}). Allele frequencies 

#'   in this set represent a hybrid index of 1. 

#' @param precols A \code{precols} object produced by \code{read.data}, or a custom character vector with 

#'   names of all pre-marker columns in the \code{data}.

#' @param INDLABEL.name Character string. The name of the column in \code{data} designated as \code{INDLABEL}. 

#'   Default is \dQuote{INDLABEL}.

#' @param POPID.name Character string. The name of the column in \code{data} designated as \code{POPID}. 

#'   Default is \dQuote{POPID}.

#' @param max.S.MAF Numeric locus filtering option. Removes loci for which the smaller of the two parental reference set minor allele 

#'   frequencies is greater than this value. Included because only one parental reference minor allele frequency needs to be close to 

#'   zero for the locus to be informative. Default is \code{0.5} (no filtering).

#' @param min.diff Numeric locus filtering option. Removes loci for which the difference in allele frequency between 

#'   parental reference sets is less than this value. Default is \code{0} (no filtering).

#' @param min.allele.copies.S0 Locus filtering option based on sample size. Removes loci for which the 

#'   number of allele copies in the \code{S0} 

#'   parental reference set, excluding missing data, is less than this value. Default is \code{0} (no filtering).

#' @param min.allele.copies.S1 Locus filtering option based on sample size. Removes loci for which the 

#'   number of allele copies in the \code{S1} 

#'   parental reference set, excluding missing data, is less than this value. Default is \code{0} (no filtering).

#' @param AF.CIoverlap Locus filtering option based on overlap of 95% credible intervals for parental allele frequencies. 

#'   If \code{FALSE} removes loci for which the allele frequency 95% credible intervals overlap for the two parental reference populations. 

#'   Default is \code{TRUE} (no filtering).

#' @param return.genotype.table Logical. Whether to create a genotype table with individual genotypes for each locus scored from 

#'   \code{0} to \code{PLOIDY}, based on their number of copies of the designated \code{S1} allele (the allele with relatively higher 

#'   frequency in the \code{S1} parental reference set). Default is \code{FALSE}.

#' @param return.locus.table Logical. Whether to create a shortened version of the output table, with one row 

#'   per locus and only containing locus-level data (such as parental allele frequencies). Default is \code{FALSE}.

#' @param marker.info.file Character string. The name of an optional file of additional marker information, to be read in 

#'   and joined to the primary output table (and to the locus table if \code{return.locus.table} set to \code{TRUE}). 

#'   May be useful for downstream statistical comparison of hybrid index or genomic cline estimates among groups of markers. 

#'   Default is \code{NULL}.

#' @param prior.prop_1 Numeric vector. Values of \code{shape1} and \code{shape2} parameters for a beta-distributed prior

#'   for parental allele frequencies (proportions) of the \code{S1_allele}. Default is \code{shape1=shape2=1} (uniform distribution).

#' @details The primary output of \code{data.prep} is a long-form \code{data.table} with 1 row for each allele copy in the 

#'   data set, to be used in downstream hybrid index and genomic cline estimation. It also optionally produces (1) a locus table,

#'   with one row per locus and containing the same locus data as in the long-form table (see \sQuote{Value}), and (2) a 

#'   genotype table, with genotypes 

#'   scored according to the number of alleles from parental reference set \code{S1}. This may be useful for example in estimating 

#'   parent-specific linkage disequilibria. \code{data.prep} uses objects produced by \code{read.data} as input, or custom 

#'   objects in the same format.

#'

#' The optional \code{marker.info.file} should be a file with one row per marker and any number of 

#'   columns (including the marker name), each for a variable of

#'   interest for grouping markers, such as chromosome, map location, gene function, SNP type (non-synonymous, synonymous, 

#'   intronic etc.). These variables 

#'   can then be used downstream for statistical comparison of hybrid index or clines among groups of loci.

#'

#' The output includes posterior mean, mode and upper and lower 95% credible intervals for frequency of the Source 1 allele

#'   for each locus, assuming beta-distributed proportions. Note that if the locus is fixed in a parental reference set,

#'   its posterior mode will be 0 or 1 and will fall outside the credible intervals due to the influence of the non-zero

#'   prior shape parameters. If prior shape parameters are set to zero (equivalent to no prior), the best posterior estimate

#'   is the mean rather than the mode, but in this case credible intervals cannot be estimated when the locus is fixed 

#'   for one allele.

#' @return list with one mandatory component, \code{data.prep}, a \code{data.table} and \code{data.frame} long-form 

#'   version of the imported data table

#'   with additional fields listed below. The two optional components are: (1) \code{geno.data}, a \code{data.table} and \code{data.frame} 

#'   containing genotypes for each locus scored from 

#'   \code{0} to \code{PLOIDY}, based on their number of copies of the designated \code{S1} allele (the allele with relatively higher 

#'   frequency in the \code{S1} parental reference set); (2) \code{locus.data}, a shortened version of the output table, with one row 

#'   per locus and containing the same locus-level data as \code{data.prep}.

#'

#'   \code{data.prep} contains the following additional fields:

#'   \item{Source}{whether the sample is from source 0 (\sQuote{S0}), source 1 (\sQuote{S1}) or is a test individual 

#'      (\sQuote{TEST}).}

#'   \item{S0_allele}{name of the allele with relatively higher frequency in \code{S0}.}

#'   \item{S1_allele}{name of the allele with relatively higher frequency in \code{S1}.}

#'   \item{S0.prop_1}{frequency (proportion) of the \code{S1_allele} in \code{S0}.}

#'   \item{S1.prop_1}{frequency (proportion) of the \code{S1_allele} in \code{S1}.}

#'   \item{S0.N_1}{number of non-missing copies of the \code{S1_allele} in \code{S0}.}

#'   \item{S1.N_1}{number of non-missing copies of the \code{S1_allele} in \code{S1}.}

#'   \item{S0.prop_0}{frequency (proportion) of the \code{S0_allele} in \code{S0}.}

#'   \item{S1.prop_0}{frequency (proportion) of the \code{S0_allele} in \code{S1}.}

#'   \item{N_allele_copies.S0}{total number of non-missing allele copies in \code{S0}.}

#'   \item{N_allele_copies.S1}{total number of non-missing allele copies in \code{S1}.}

#'   \item{Source.afdiff}{difference in allele frequency between \code{S0} and \code{S1}.}

#'   \item{S0.MAF}{\code{S0} minor allele frequency.}

#'   \item{S1.MAF}{\code{S1} minor allele frequency.}

#'   \item{min.MAF}{smaller of \code{S0.MAF} and \code{S1.MAF}.}

#'   \item{Source_allele}{numeric. Whether the allele is the \code{S0_allele} (\sQuote{0}) or \code{S1_allele} (\sQuote{1}).}

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

#' @author

#'   Richard Ian Bailey \email{richardianbailey@@gmail.com}

#' @examples

#'

#' \dontrun{

#' #Using objects produced by \code{read.data}:

#' prepdata=data.prep(data=dat$data,

#' loci=dat$loci,

#' alleles=dat$alleles,

#' S0=c("pop1","pop2"),	#Two populations in the POPID column chosen as parental reference set S0#

#' S1="pop3",			#One population in the POPID column chosen as parental reference set S1#

#' precols=dat$precols,

#' max.S.MAF = 0.1,		#Filtering by parental minor allele frequency only; other filtering options below are same as default#

#' min.diff=0,

#' min.allele.copies.S0=0,

#' min.allele.copies.S1=0,

#' AF.CIoverlap=TRUE,

#' return.genotype.table=T,

#' return.locus.table=T)

#' }

#' @export

data.prep.new2= function(data,loci,alleles,S0,S1,precols,

    INDLABEL.name="INDLABEL",POPID.name="POPID",

    max.S.MAF=0.5,min.diff=0,min.allele.copies.S0=0,min.allele.copies.S1=0,AF.CIoverlap=TRUE,

    return.genotype.table=FALSE,return.locus.table=FALSE,

    marker.info.file=NULL,
	
	prior.prop_1=c(1,1)){



    cat(paste("Preparing data for hybrid index and genomic cline estimation..."),

        fill=1); flush.console();



    if(is.null(marker.info.file)==F){

        marker.info=fread(marker.info.file);

                                    };



    setkeyv(data,POPID.name);

    data[,Source:="TEST"];

    data[S0,Source:="S0"];

    data[S1,Source:="S1"];



    data[,loci] = data[,lapply(.SD,replace.fun.A),.SDcols=loci[]];

    data[,loci] = data[,lapply(.SD,replace.fun),.SDcols=loci[]];



    DT.EXTRACOL=unique(data[,1:length(precols)]);setkeyv(DT.EXTRACOL,INDLABEL.name);

	

    setkey(data,Source);



    keepNloc=as.integer(2.5e+07/nrow(data["S0"]));

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

    DT.S[propdiff_A < 0, S0_allele := "A"];

    setkey(DT.S,propdiff_B);

    DT.S[propdiff_B < 0, S0_allele := "B"];



    setkey(DT.S,propdiff_A);

    DT.S[propdiff_A > 0, S1_allele := "A"][propdiff_A > 0, S0.prop_1 := prop_A][propdiff_A > 0, 

        S1.prop_1 := i.prop_A][propdiff_A > 0, S0.N_1 := N_A][propdiff_A > 0, 

        S1.N_1 := i.N_A];



    setkey(DT.S,propdiff_B);

    DT.S[propdiff_B > 0, S1_allele := "B"][propdiff_B > 0, S0.prop_1 := prop_B][propdiff_B > 0, 

        S1.prop_1 := i.prop_B][propdiff_B > 0, S0.N_1 := N_B][propdiff_B > 0, 

        S1.N_1 := i.N_B];



    DT.S[,S0.prop_0:= (1-S0.prop_1)][,S1.prop_0:= (1-S1.prop_1)];



    DT.S[,N_allele_copies.S0 := N_A + N_B];

    DT.S[,N_allele_copies.S1 := i.N_A + i.N_B];

    DT.S[,Source.afdiff:= S1.prop_1 - S0.prop_1];



    setkey(DT.S,S0.prop_1,S1.prop_1);

    DT.S[,S0.MAF:= 0.5 - abs(0.5 - S0.prop_1)];

    DT.S[,S1.MAF:= 0.5 - abs(0.5 - S1.prop_1)];



    setkey(DT.S,S0.MAF,S1.MAF);

    DT.S[S0.MAF<S1.MAF,min.MAF:=S0.MAF];

    DT.S[S0.MAF>=S1.MAF,min.MAF:=S1.MAF];



    DT.S = DT.S[min.MAF <= max.S.MAF];



    setkey(DT.S,S0.prop_1,S1.prop_1);

    DT.S = DT.S[S0.prop_1 != "NA"][(S1.prop_1 - S0.prop_1) >= min.diff];

    setkey(DT.S,N_allele_copies.S0,N_allele_copies.S1);

    DT.S = DT.S[N_allele_copies.S0 >= min.allele.copies.S0][N_allele_copies.S1 >= 

        min.allele.copies.S1];



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



    setkey(DT.2,locus); setkey(DT.S,locus);

    DT.join = DT.2[DT.S];



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



    if(return.genotype.table==TRUE){

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



        if(is.null(marker.info.file)==F){

            setkey(locus.data,locus);setkey(marker.info,locus);

            locus.data=locus.data[marker.info];

                                        };



                                };



    if(is.null(marker.info.file)==F){

        setkey(DT.join,locus);setkey(marker.info,locus);

        DT.join=DT.join[marker.info];

                                    };



    output=list();

        output$data.prep = DT.join[is.na(Source_allele)==F,
		
		    .(INDLABEL,POPID,locus,Source,S0.prop_1,S1.prop_1,S0.prop_0,S1.prop_0,Source_allele,Source.afdiff)];
		
		output$Nloci.postfilter = length(DT.join[,unique(locus)]);



    if(return.genotype.table==TRUE){

        output$geno.data = geno;

                                   };



    if(return.locus.table==TRUE){

        output$locus.data = locus.data;

                            };



    cat(paste("Done."),fill=1);



    return(output)

                                                          }