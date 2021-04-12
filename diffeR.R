#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      diffeR.R - Set diff between 2 files by column
 
      Arguments:
      --x=someValue            - input file x path
      --y=someValue            - input file y path
      --by.x=someValue         - by.x parameter
      --by.y=someValue         - by.y parameter
      --noh.x                  - file x doesn't have header
      --noh.y                  - file y doesn't have header
      --out=someValue          - output file
      --skip.x                 - skip n lines from x file
      --skip.y                 - skip n files from y file
      --out.label=someValues   - output header labels
      --print.out.label        - print header labels in output file (Enabled by default if --out.label is used)
      --help                   - print this Help

      Example:
      ./diffeR.R --x=\"input1.txt\" --y=\"input2.txt\" --by.x=\"ID1\" --by.y=\"ID2\" --out=\"output.txt\" 
      
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


if( is.null(argsL[['skip.x']])) {
	argsL[['skip.x']] = 0
}
if( is.null(argsL[['skip.y']])) {
	argsL[['skip.y']] = 0
}


x.df <- read.delim(argsL[['x']], header=argsL[['noh.x']], stringsAsFactors=FALSE,skip=argsL[['skip.x']], sep="\t")
y.df <- read.delim(argsL[['y']], header=argsL[['noh.y']], stringsAsFactors=FALSE,skip=argsL[['skip.y']], sep="\t")


#print(dim(x.df))
#print(dim(y.df))

#print(head(x.df))
#print(head(y.df))

xy.diff <- as.data.frame(setdiff( x.df[[argsL[['by.x']]]], y.df[[argsL[['by.y']]]] ));

#print(xy.diff)

write.col.names <- FALSE
if(! is.null(argsL[['print.out.label']])) {
	write.col.names <- TRUE
}


if(! is.null(argsL[['out.label']])) {
	out.label <- unlist(strsplit(argsL[['out.label']],","))
	if ( length(out.label) != 1) {
        	sink(stderr())
	        cat("\nERROR: output label must have only one column identifier !\n\n")
        	sink()
	        q(save="no")
	}
	colnames(xy.diff) <- out.label
	write.col.names <- TRUE
}

write.table(xy.diff, file=argsL[['out']], col.names=write.col.names, row.names=FALSE, quote=FALSE, sep="\t")

