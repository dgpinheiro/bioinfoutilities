#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      aggregatoR.R - Aggregate row values by column
 
      Arguments:
      --x=someValue            - input file x path
      --sep=someValue          - to concatenate values by character
      --by.x=someValue         - to aggregate values by column name
      --colnames.x=someValues  - comma separated colnames for x
      --noh.x                  - file x doesn't have header
      --skip.x                 - skip n lines from x file
      --out=someValue          - output file
      --out.label=someValues   - output header labels
      --help                   - print this Help

      Example:
      ./aggregatoR.R --x=\"input1.txt\" --by.x=\"ID1\" --out=\"output.txt\" 
      
      Daniel Guariz Pinheiro
      FCAV/UNESP - Univ Estadual Paulista
      dgpinheiro@gmail.com
      \n\n")

  q(save="no")
}


suppressPackageStartupMessages( library('dplyr') )

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL[['x']])) {
        sink(stderr())
        cat("\nERROR: Missing input file x !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['by.x']])) {
        sink(stderr())
        cat("\nERROR: Missing by.x column !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['out']])) {
        sink(stderr())
        cat("\nERROR: Missing output file !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['noh.x']])) {
	argsL[['noh.x']] = TRUE
} else {
	argsL[['noh.x']] = FALSE
}

if(is.null(argsL[['sep']])) {
	argsL[['all.x']] = ";"
}

if( is.null(argsL[['skip.x']])) {
	argsL[['skip.x']] = 0
}


x.df <- read.delim(argsL[['x']], header=argsL[['noh.x']], stringsAsFactors=FALSE,skip=argsL[['skip.x']], sep="\t")



if(! is.null(argsL[['colnames.x']])) {
	colnames.x <- unlist(strsplit(argsL[['colnames.x']],","))
        if (dim(x.df)[2] > length(colnames.x)) {
                sink(stderr())
                cat("\nERROR: colnames.x must have more elements !\n\n")
                sink()
                q(save="no")
        }
	colnames(x.df) <- colnames.x
}

#print(dim(x.df))
#print(dim(y.df))

#print(head(x.df))
#print(head(y.df))


x.distinct.df <- distinct(x.df, .keep_all=TRUE)
x.agg.df <- aggregate(as.formula(paste(". ~",argsL[['by.x']])), x.distinct.df, paste, collapse=argsL[['sep']])

if(! is.null(argsL[['out.label']])) {
	out.label <- unlist(strsplit(argsL[['out.label']],","))
	if ( length(colnames.x) != length(out.label)) {
        	sink(stderr())
	        cat("\nERROR: output labels must have the same number of selected colnames.x !\n\n")
        	sink()
	        q(save="no")
	}
	colnames(x.agg.df) <- out.label
}


write.table(x.agg.df, file=argsL[['out']], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

