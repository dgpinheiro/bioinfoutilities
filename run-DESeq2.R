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
      --ignore=someValue	- columns to ignore
      --pseudocount		- number to add
      --help            	- print this Help

	GROUP file must be formatted as tab-delimited columns with sample names grouped by group name:

	id		name		group
	<SAMPLE_ID>	<SAMPLE_NAME>	<GROUP_NAME>

	* <SAMPLE_ID> is the name of the sample column in abundances' matrix
	* <SAMPLE_NAME> is the name of the sample column after DESeq2 analysis
	* <GROUP_NAME> is the name of the group

      Example:
      ./run-DESeq2.R --ignore=\"id,desc,other\" --in=\"input.txt\" --out=\"./output\" 
      
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

if(is.null(argsL[['ncores']])) {
	argsL[['ncores']] = 1
}

if(is.null(argsL[['pseudocount']])) {
	argsL[['pseudocount']] = 0
}

ignorecols=c()
if( ! is.null(argsL[['ignore']])) {
	ignorecols=c(unlist(strsplit(argsL[['ignore']],split=",")))
}

suppressPackageStartupMessages( library("BiocParallel") )

register(MulticoreParam(argsL[['ncores']]))

suppressPackageStartupMessages( library('methods') )
suppressPackageStartupMessages( library('DESeq2') )
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

# Low abundance filtering
#cat("Filtering low abundance reads...\n")
#keep <- rowSums(cpm(x[,samps])>1) >= 5
#x.filtered <- x[keep, ]
x.filtered <- x

countData <- x.filtered[ , as.character(g$id)]

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

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds <- estimateDispersions(dds)

normFactors <- matrix(runif(nrow(dds)*ncol(dds),0.5,1.5),
                      ncol=ncol(dds),nrow=nrow(dds),
                      dimnames=list(1:nrow(dds),1:ncol(dds)))

normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors


cat("DESeq...\n")
dds <- DESeq(dds, parallel=TRUE)

deseq_ncounts <- as.data.frame( counts(dds, normalized=TRUE) )
deseq_ncounts[[idcol]] <- rownames(deseq_ncounts)


design <- model.matrix(~0+group, data=colData)
colnames(design) <- levels(colData$group)


all.comb <- apply(t(combn(levels(colData$group),2)), 1, function(x) { paste( x ,collapse="-") })

cat("Contrast results...\n")
res <- list()
for (c in 1:length(all.comb)) {
	print(all.comb[c])
	#my.contrast <- makeContrasts(contrasts=all.comb[c], levels=design)
	one <- unlist(strsplit(all.comb[c],'-'))[1]
	two <- unlist(strsplit(all.comb[c],'-'))[2]
	#save(list = ls(all=TRUE,envir=globalenv()), file = "run-DESeq2.RData")
	
	res[[ all.comb[c] ]] <- results(dds, contrast=c("group",one,two))

	dgeres.df <- as.data.frame(  res[[ all.comb[c] ]]  );
	dgeres.df[[idcol]] <- rownames(dgeres.df)

	rownames(dgeres.df) <- dgeres.df[[idcol]]
	colnames(dgeres.df) <- c("baseMean","logFC","lfcSE","stat","PValue", "FDR", idcol)
	
	dge.res <- (merge( dgeres.df, deseq_ncounts, by=idcol))[,c(idcol, g$name, 'logFC', 'PValue', 'FDR')]
	rownames(dge.res) <- dge.res[[idcol]]
		
	samps.reordered <- c()
	for (l in levels(colData$group)) { 
		#print(l) 
		samps.reordered <- c(samps.reordered, g[which(g$group==l), 'name'])
	}

	for (col in ignorecols) {
		dge.res[[col]] <- x.filtered[rownames(dge.res), col]
	}

	dge.res$FDR <- p.adjust(dge.res$PValue, method="fdr")

	write.table(file=paste(argsL[['out']],"/",all.comb[c],".txt",sep=""), x=dge.res, row.names=FALSE, quote=FALSE, sep="\t")
}

#save(list = ls(all=TRUE,envir=globalenv()), file = "run-DESeq2.RData")

