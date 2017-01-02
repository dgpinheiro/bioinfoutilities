#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      mergeR.R - Merge 2 files by column
 
      Arguments:
      --x=someValue            - input file x path
      --y=someValue            - input file y path
      --by.x=someValue         - by.x parameter
      --by.y=someValue         - by.y parameter
      --noh.x                  - file x doesn't have header
      --noh.y                  - file y doesn't have header
      --colname.x=someValues   - comma separated colnames for x
      --colname.y=someValues   - comma separated colnames for y              
      --out=someValue          - output file
      --help                   - print this Help

      Example:
      ./mergeR.R --x=\"input1.txt\" --y=\"input2.txt\" --by.x=\"ID1\" --by.y=\"ID2\" --out=\"output.txt\" 
      
      Daniel Guariz Pinheiro
      FCAV/UNESP - Univ Estadual Paulista
      dgpinheiro@gmail.com
      \n\n")

  q(save="no")
}


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

if(is.null(argsL[['y']])) {
        sink(stderr())
        cat("\nERROR: Missing input file y !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['by.x']])) {
        sink(stderr())
        cat("\nERROR: Missing by.x column !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['by.y']])) {
        sink(stderr())
        cat("\nERROR: Missing by.y column !\n\n")
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

if(is.null(argsL[['noh.y']])) {
	argsL[['noh.y']] = TRUE
} else {
	argsL[['noh.y']] = FALSE
}


x.df <- read.delim(argsL[['x']], header=argsL[['noh.x']], stringsAsFactors=FALSE)
y.df <- read.delim(argsL[['y']], header=argsL[['noh.y']], stringsAsFactors=FALSE)

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
if(! is.null(argsL[['colnames.y']])) {
	colnames.y <- unlist(strsplit(argsL[['colnames.y']],","))
	if (dim(y.df)[2] > length(colnames.y)) {
        	sink(stderr())
	        cat("\nERROR: colnames.y must have more elements !\n\n")
        	sink()
	        q(save="no")
	}
	colnames(y.df) <- colnames.y
}

#print(dim(x.df))
#print(dim(y.df))

#print(head(x.df))
#print(head(y.df))

xy.df <- merge(x=x.df, y=y.df, by.x=argsL[['by.x']], by.y=argsL[['by.y']])

write.table(xy.df, file=argsL[['out']], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

