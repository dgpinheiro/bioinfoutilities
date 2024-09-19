#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      filteR.R - filter cuffdiff results
 
      Arguments:
      --in=someValue   		- input file path of normalized counts by cuffnorm (*)
      --out=someValue   	- output file path (*)
      --minsum			- minimum value for the sum of sample's FPKM
      --help            	- print this Help
 
      (*) mandatory

      Example:
      ./filter.R [ --minsum=0 ] --in=\"cuffnorm/gene_exp.diff\" --out=\"output.txt\" 
      
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

if(is.null(argsL[['out']])) {
	sink(stderr())
	cat("\nERROR: Missing output file!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['minsum']])) {
	argsL[['minsum']]=0
}


#suppressPackageStartupMessages( library('methods') )
#suppressPackageStartupMessages( library('edgeR') )


x <-read.delim( argsL[['in']], header=TRUE, stringsAsFactors=FALSE )

idcol <- 'test_id'
samps <- setdiff(colnames(x), idcol )

rownames(x) <- paste(x[[idcol]],x[['sample_1']],x[['sample_2']],sep=".")

x.filt <- subset(x, (value_1+value_2) > argsL[['minsum']] )
x.filt$q_value.2 <- p.adjust(x.filt$p_value, method="fdr")

write.table(x.filt[,c(idcol,samps)], file=argsL[['out']], quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

