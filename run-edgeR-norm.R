#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      run-edgeR-norm.R - Run edgeR normalization method
 
      Arguments:
      --in=someValue   		- input file path (*)
      --groups=someValue	- group file path
      --ignore=someValue	- columns to ignore
      --out=someValue   	- output file path (*)
      --help            	- print this Help
 
      (*) mandatory

      Example:
      ./run-edgeR-norm.R [ --ignore=\"id,desc,other\" --groups=\"groups.txt\" ] --in=\"input.txt\" --out=\"output.txt\" 
      
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
	

if(is.null(argsL[['in']])) {
	sink(stderr())
	cat("\nERROR: Missing input file!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['out']])) {
	sink(stderr())
	cat("\nERROR: Missing output file!\n\n")
	sink()
	q(save="no")
}
if ( ! file.exists(argsL[['in']]) ) {
	sink(stderr())
	cat("\nERROR: Input file not found!\n\n")
	sink()
	q(save="no")
}
ignorecols=c()
if( ! is.null(argsL[['ignore']])) {
	ignorecols=c(unlist(strsplit(argsL[['ignore']])))
}


suppressPackageStartupMessages( library('methods') )
suppressPackageStartupMessages( library('edgeR') )


x <-read.delim( argsL[['in']], header=TRUE )


for (col in colnames(x)) { 
	if ( ! is.numeric(x[[col]]) ) {
		ignorecols <- unique(c(ignorecols, col))
	}
}

idcol <- colnames(x)[1]
samps <- setdiff(colnames(x), ignorecols )

if(is.null(argsL[['groups']])) {
	g <-data.frame( id=samps, name=samps, group=samps, stringsAsFactors=FALSE )
} else {
	g <-read.delim( argsL[['groups']], header=TRUE, stringsAsFactor=FALSE )
}


if ( any(! ( g$id %in% samps )) ) {
	sink(stderr())
	cat(paste("\nERROR: id(s) (", paste(c(g$id)[ ! ( g$id %in% samps )  ], collapse=","), ") not found in data!\n\n", sep=""))
	sink()
	q(save="no")
}


rownames(x) <- x[[idcol]]

df <- x[ , as.character(g$id)]

rownames(g) <- g$id


colnames(df)[ colnames(df) %in% g$id ] <- as.character(g[['name']])

group <- factor( g$group )

y <- DGEList(counts=df[, g$name], group=group)

y <- calcNormFactors(y,method="TMM")

df.cpm <- as.data.frame(cpm(y, log=FALSE, prior.count=0.25, normalized.lib.sizes=TRUE))

samps.reordered <- c()
for (l in levels(group)) { 
	#print(l) 
	samps.reordered <- c(samps.reordered, g[which(g$group==l), 'name'])
}

for (col in ignorecols) {
	df.cpm[[col]] <- x[rownames(df.cpm), col]
}

write.table(df.cpm[,c(ignorecols,samps.reordered)], file=argsL[['out']], quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
