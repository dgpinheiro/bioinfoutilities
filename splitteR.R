#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      splitteR.R - Split strings by separator in a column
 
      Arguments:
      --x=someValue            - input file x path
      --by.x=someValue         - string delimiter
      --col.x=someValue        - column in x
      --colnames.x=someValues  - comma separated colnames for x
      --noh.x                  - file x doesn't have header
      --out=someValue          - output file
      --skip.x                 - skip n lines from x file
      --out.label=someValues   - output header labels
      --no.out.label           - remove output header labels
      --help                   - print this Help

      Example:
      ./splitteR.R --x=\"input1.txt\" --col.x=\"column name\" --by.x=\"delimiter\" --out=\"output.txt\" 
      
      Daniel Guariz Pinheiro
      FCAV/UNESP - Univ Estadual Paulista
      dgpinheiro@gmail.com
      \n\n")

  q(save="no")
}

suppressPackageStartupMessages( library('tidyr') )

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
        cat("\nERROR: Missing by.x delimiter !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['col.x']])) {
        sink(stderr())
        cat("\nERROR: Missing col.x column !\n\n")
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

if( is.null(argsL[['skip.x']])) {
	argsL[['skip.x']] = 0
}

if( is.null(argsL[['no.out.label']])) {
	argsL[['no.out.label']] = TRUE
} else {
	argsL[['no.out.label']] = FALSE
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

#print(head(x.df))



x.new.df <- separate_rows(x.df, grep(paste('^',argsL[['col.x']],'$',sep=''), colnames(x.df)), sep=argsL[['by.x']])


if(! is.null(argsL[['out.label']])) {
	out.label <- unlist(strsplit(argsL[['out.label']],","))
	if ( length(colnames.x) == length(out.label)) {
        	sink(stderr())
	        cat("\nERROR: output labels must have the same number of selected colnames.x !\n\n")
        	sink()
	        q(save="no")
	}
	colnames(x.new.df) <- out.label
}

write.table(x.new.df, file=argsL[['out']], col.names=argsL[['no.out.label']], row.names=FALSE, quote=FALSE, sep="\t")

