#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      run-DESeq2.R - Run DESeq2 method
 
      Arguments:
      --in=someValue   		- input file path
      --groups=someValue	- group file path
      --out=someValue   	- output directory path
      --ncores=Integer		- number of processors available
      --ignore=someValue	- columns to ignore (comma separated values). Example: --ignore=\"id,desc,other\"
      --pseudocount		- number to add
      --controlgenes		- list of \"control genes\" file (the gene identifiers must be in the input file)
      --refgroup		- Reference sample group (this group must be in the groups file)
      --lfcshrink		- LFC shrink (\"apeglm\", \"ashr\", \"normal\")
							 
      --help            	- print this Help
     
	GROUP file must be formatted as tab-delimited columns with sample names grouped by group name:

	id		name		group
	<SAMPLE_ID>	<SAMPLE_NAME>	<GROUP_NAME>

	* <SAMPLE_ID> is the name of the sample column in abundances' matrix
	* <SAMPLE_NAME> is the name of the sample column after DESeq2 analysis
	* <GROUP_NAME> is the name of the group

      Example:
      ./run-DESeq2.R --groups=\"groups.txt\" --in=\"input.txt\" --out=\"./output\" [ --controlgenes=\"control_genes.txt\" --refgroup=\"CTL\" --lfcshrink=\"apeglm\" ]
      
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
	cat("\nERROR: Missing output dir!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['groups']])) {
	sink(stderr())
	cat("\nERROR: Missing groups file!\n\n")
	sink()
	q(save="no")
}
controlGenes<-c()
if( ! is.null(argsL[['controlgenes']])) {
	if ( ! file.exists(argsL[['controlgenes']]) ) {
		sink(stderr())
		cat("\nERROR: List of \"control genes\" file not found!\n\n")
		sink()
		q(save="no")
	}
	
	controlGenes <- unique( read.delim(file=argsL[['controlgenes']], header=FALSE,stringsAsFactors=FALSE)$V1 )

}

if ( ! is.null(argsL[['lfcshrink']]) ) {
	if ( ! argsL[['lfcshrink']] %in% c("apeglm", "ashr", "normal") ) {
		sink(stderr())
		cat("\nERROR: Wrong lfcShrink type !\n\n")
		sink()
		q(save="no")
	}
}

if ( ! file.exists(argsL[['in']]) ) {
	sink(stderr())
	cat("\nERROR: Input file not found!\n\n")
	sink()
	q(save="no")
}
if ( (! file.exists(argsL[['out']]))|(! file.info(argsL[['out']])[1,"isdir"] ) ) {
	sink(stderr())
	cat("\nERROR: Output directory not found!\n\n")
	sink()
	q(save="no")
}
if ( ! file.exists(argsL[['groups']]) ) {
	sink(stderr())
	cat("\nERROR: Groups file not found!\n\n")
	sink()
	q(save="no")
}

if(! is.null(argsL[['refgroup']])) {
	sink(stderr())
	cat(paste("\nWARNING: Using reference sample group:",argsL[['refgroup']],"\n\n"))
	sink()
}

if(is.null(argsL[['ncores']])) {
	argsL[['ncores']] = 1
}

if(is.null(argsL[['pseudocount']])) {
	argsL[['pseudocount']] = 1
}

ignorecols=c()
if( ! is.null(argsL[['ignore']])) {
	ignorecols=c(unlist(strsplit(argsL[['ignore']],split=",")))
}

suppressPackageStartupMessages( library("BiocParallel") )

register(MulticoreParam(argsL[['ncores']]))

suppressPackageStartupMessages( library('methods') )
suppressPackageStartupMessages( library('DESeq2') )

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

	if(! is.null(argsL[['refgroup']])) {
		if (! any(g$group==argsL[['refgroup']])) {
			sink(stderr())
			cat(paste("\nERROR: Not found reference group (",argsL[['refgroup']], ") in groups file (",argsL[['groups']],") !\n\n", sep=""))
			sink()
			q(save="no")
		}
	}
}

if (any(!colnames(g) %in% c('id', 'name', 'group'))) {
	sink(stderr())
	cat(paste("\nERROR: Not found in groups header (", paste(colnames(g), collapse=", "), ") the required header columns (id, name, group) using TAB delimiter!\n\n", sep=""))
	sink()
	q(save="no")
}

if ( any(! ( g$id %in% samps )) ) {
	sink(stderr())
	cat(paste("\nERROR: id(s) (", paste(c(g$id)[ ! ( g$id %in% samps )  ], collapse=","), ") not found in data!\n\n", sep=""))
	sink()
	q(save="no")
}


rownames(x) <- x[[idcol]]

countData <- x[ , as.character(g$id)]

rownames(g) <- g$id

colnames(countData)[ colnames(countData) %in% g$id ] <- as.character(g[['name']])

colData<- g[,setdiff(colnames(g),'name')]

rownames(colData) <- g$name
colData$id <- g$name
colData$group <- as.factor( colData$group )

countData <- round(countData)


dds <- DESeqDataSetFromMatrix(countData = countData+as.numeric(argsL[['pseudocount']]),
                              colData = colData,
                              design = ~ group)

if (length(controlGenes) > 0) {
	
	sink(stderr())
	cat(paste("\nWARNING: Estimating size factor using control genes: ",paste(controlGenes,collapse=", "),"\n\n"))
	sink()

	isControl <- rownames(dds) %in% controlGenes
	
	dds <- estimateSizeFactors(dds, controlGenes=which(isControl) )
} else {
	dds <- estimateSizeFactors(dds)
}	
sizeFactors(dds)
dds <- estimateDispersions(dds) 

#if (is.na(tryCatch( 
#		dds <- estimateDispersions(dds) 
#		  , 
#		  error=function(e) NA 
#		  ) )) {
#		dds <- estimateDispersionsGeneEst(dds)
#		dispersions(dds) <- mcols(dds)$dispGeneEst
#}
#
#cat("DESeq2...\n")
#dds <- DESeq(dds, parallel=TRUE)
dds = nbinomWaldTest(dds)

deseq_ncounts <- as.data.frame( counts(dds, normalized=TRUE) )
deseq_ncounts[[idcol]] <- rownames(deseq_ncounts)

mean.cols <- c()
for (lvl in levels(colData$group) ) {
	deseq_ncounts[[paste(lvl,'mean',sep=".")]] <- rowMeans(deseq_ncounts[, g[g$group==lvl, 'name'] ])
	mean.cols <- c(mean.cols, paste(lvl,'mean',sep="."))
}

design <- model.matrix(~0+group, data=colData)
colnames(design) <- levels(colData$group)


if(is.null(argsL[['refgroup']])) {
	all.comb <- apply(t(combn(levels(colData$group),2)), 1, function(x) { paste( x ,collapse="-") })
} else {
	all.comb <- paste(setdiff(levels(colData$group), argsL[['refgroup']]), argsL[['refgroup']], sep="-")
}

cat("Contrast results...\n")
res <- list()
for (c in 1:length(all.comb)) {
	
	print(all.comb[c])

	one <- unlist(strsplit(all.comb[c],'-'))[1]
	two <- unlist(strsplit(all.comb[c],'-'))[2]
	
	if ( is.null(argsL[['lfcshrink']]) ) {
		res[[ all.comb[c] ]] <- results(dds, contrast=c("group",one,two))
		dgeres.df <- as.data.frame(  res[[ all.comb[c] ]]  )
		dgeres.df[[idcol]] <- rownames(dgeres.df)

		rownames(dgeres.df) <- dgeres.df[[idcol]]
		colnames(dgeres.df) <- c("baseMean","logFC","lfcSE","stat","PValue", "FDR", idcol)
	} else {
		res[[ all.comb[c] ]] <- lfcShrink(dds, contrast=c("group",one,two), type=argsL[['lfcshrink']])
		dgeres.df <- as.data.frame(  res[[ all.comb[c] ]]  )
		dgeres.df[[idcol]] <- rownames(dgeres.df)

		rownames(dgeres.df) <- dgeres.df[[idcol]]
		colnames(dgeres.df) <- c("baseMean","logFC","lfcSE","PValue", "FDR", idcol)
	}

	dge.res <- (merge( dgeres.df, deseq_ncounts, by=idcol))[,c(idcol, g$name, mean.cols, 'logFC', 'lfcSE', 'PValue', 'FDR')]
	rownames(dge.res) <- dge.res[[idcol]]
		
	write.table(file=paste(argsL[['out']],"/",all.comb[c],".txt",sep=""), x=dge.res, row.names=FALSE, quote=FALSE, sep="\t")

}

save(list = ls(all=TRUE,envir=globalenv()), file = paste(argsL[['out']],"/","run-DESeq2.RData", sep=""))

