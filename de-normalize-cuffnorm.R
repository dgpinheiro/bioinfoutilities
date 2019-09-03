#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      de-normalize-cuffnorm.R - \"De-normalization\" of cuffnorm raw counts
 
      Arguments:
      --in=someValue   		- input file path of normalized counts by cuffnorm (*)
      --st=someValue		- input file path of samples table (*)
      --out=someValue   	- output file path (*)
      --help            	- print this Help
 
      (*) mandatory

      Example:
      ./de-normalize-cuffnorm.R [ --in=\"cuffnorm/genes.count_table\" --st=\"cuffnorm/samples.table\" --out=\"output.txt\" 
      
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
	cat("\nERROR: Missing input file (genes.count_table)!\n\n")
	sink()
	q(save="no")
}
if ( ! file.exists(argsL[['in']]) ) {
	sink(stderr())
	cat("\nERROR: Input file not found (genes.count_table)!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['st']])) {
	sink(stderr())
	cat("\nERROR: Missing input file (samples.table)!\n\n")
	sink()
	q(save="no")
}
if ( ! file.exists(argsL[['st']]) ) {
	sink(stderr())
	cat("\nERROR: Input file not found (samples.table)!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['out']])) {
	sink(stderr())
	cat("\nERROR: Missing output file!\n\n")
	sink()
	q(save="no")
}


#suppressPackageStartupMessages( library('methods') )
#suppressPackageStartupMessages( library('edgeR') )


x <-read.delim( argsL[['in']], header=TRUE, stringsAsFactors=FALSE )

idcol <- colnames(x)[1]
samps <- setdiff(colnames(x), idcol )

rownames(x) <- x[[idcol]]
#x <- x[,samps]

g <-read.delim( argsL[['st']], header=TRUE, stringsAsFactor=FALSE )
st.idcol <- colnames(g)[1]

rownames(g) <- g[[st.idcol]]

denorm <- data.frame(round(as.matrix(x[rownames(x),samps]) %*% as.matrix( diag(as.numeric( g[samps,'internal_scale'])) )))
rownames(denorm) <- rownames(x)
colnames(denorm) <- samps
denorm[[idcol]] <- rownames(denorm)

write.table(denorm[,c(idcol,samps)], file=argsL[['out']], quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

