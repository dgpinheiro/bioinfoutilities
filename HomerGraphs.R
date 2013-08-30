####################################################################################################################
### R Script to generate plots of Homer Quality Control Analysis                                                 ###
###    [http://biowhat.ucsd.edu/homer/ngs/tagDir.html]                                                           ###
###                                                                                                              ###
### Author(s):                                                                                                   ###
###           Daniel Guariz Pinheiro                                                                             ###
### 	      Simon Coetzee                                                                                      ###
###                                                                                                              ###
### Example:                                                                                                     ###
###           cat HomerGraphs.R | R --no-restore --no-save --silent --args LNCaP_AR /tmp > /dev/null             ###
###                                                                                                              ###
####################################################################################################################

# Required libraries

library(ggplot2)
require(reshape2)
require(scales)

# Get command-line arguments and validate them

args <- commandArgs(trailingOnly = TRUE)

if (is.na(args[1])) {
	sink(stderr())
	cat("ERROR: Missing input directory (Homer TagDirectory) !\n")
	sink()
	quit(save="no", status=10);
} else {
	if (!file.exists(args[1])){
		sink(stderr())
		cat(paste(c("ERROR: Input directory (",args[1],") not exists !\n"),collapse=""))
		sink()
		quit(save="no", status=10);
	}

}

if (is.na(args[2])) {
	args[2]="."
} else {
	if (!file.exists(args[2])){
		sink(stderr())
		cat(paste(c("ERROR: Output directory (",args[2],") not exists !\n"),collapse=""))
		sink()
		quit(save="no", status=10);
	}
}


# Tag Length Distribution
if ( file.exists( paste(args[1],"tagLengthDistribution.txt",sep="/") ) ) {
	rld <- read.delim(file=paste(args[1],"tagLengthDistribution.txt",sep="/"), header=T, sep="\t")
	colnames(rld) <- c("tagLength", "fractionTags")
	prld <- ggplot(rld, aes(x = factor(tagLength), y=fractionTags))
	prld + geom_bar(stat="identity") + xlab("Tag Length") + ylab("Fraction of Tags") + ggtitle("Tag Length Distribution")
	ggsave(file=paste(args[2],"rld.png",sep="/"))
}

## Tag Count Distribution
if ( file.exists( paste(args[1],"tagCountDistribution.txt",sep="/") ) ) {
	tcd <- read.delim(file=paste(args[1],"tagCountDistribution.txt",sep="/"), header=T, sep="\t")
	colnames(tcd) <- c("readsPerPosition", "fractionOfTotalReads")
	ptcd <- ggplot(tcd[1:11,], aes(x = factor(readsPerPosition), y = fractionOfTotalReads))
	ptcd + geom_bar(stat="identity") + xlab("Reads per position") + ylab("Fraction of Total Reads") + ggtitle("Tag Count Distribution")
	ggsave(file=paste(args[2],"tcd.png",sep="/"))
}
## Tag Autocorrelation
if ( file.exists( paste(args[1],"tagAutocorrelation.txt",sep="/") ) ) {
	ta <- read.delim(file=paste(args[1],"tagAutocorrelation.txt",sep="/"), header=T, sep="\t")
	colnames(ta) <- c("distanceInBp", "sameStrand", "oppositeStrand")
	tam <- melt(ta, id.var="distanceInBp")
	tam <- tam[tam$distanceInBp > -1000 & tam$distanceInBp < 1000, ]
	colnames(tam) <- c("relativeDistanceInBp", "strand", "totalReadPairs")
	ptam <- ggplot(tam, aes(x=relativeDistanceInBp, y=totalReadPairs, group=strand, fill=strand, width=1))
	ptam + geom_bar(stat="identity", position = "identity", alpha = 0.3) + geom_vline(xintercept = 0) + xlab("Relative Distance between Reads (bp)") + ylab("Total Read Pairs") + ggtitle("Tag Autocorrelation") +
	scale_fill_discrete(name="Source", breaks=c("sameStrand","oppositeStrand"), labels=c("Same strand", "Opposite strand"))
	ggsave(file=paste(args[2],"ta.png",sep="/"))
}

## Genomic Nucleotide Frequency relative to read positions
if ( file.exists( paste(args[1],"tagFreq.txt",sep="/") ) ) {
	tf <- read.delim(file=paste(args[1],"tagFreq.txt",sep="/"), header=T, sep="\t")
	tfm <- melt(tf[,1:5], id.var="Offset")
	colnames(tfm) <- c("distanceFrom5primeEndOfRead", "Base", "nucleotideFreq")
	ptfm <- ggplot(tfm, aes(x=distanceFrom5primeEndOfRead, y=nucleotideFreq, group=Base, color=Base))
	ptfm + geom_line() + ylim(0, .4) + xlab("Distance from 5' end of Reads") + ylab("Nucleotide Frequency") + ggtitle("Genomic Nucleotide Frequency relative to read positions")
	ggsave(file=paste(args[2],"tf.png",sep="/"))
}

## Genomic Nucleotide Frequency relative to read positions (Unique)
if ( file.exists( paste(args[1],"tagFreqUniq.txt",sep="/") ) ) {
	tfu <- read.delim(file=paste(args[1],"tagFreqUniq.txt",sep="/"), header=T, sep="\t")
	tfum <- melt(tfu[,1:5], id.var="Offset")
	colnames(tfum) <- c("distanceFrom5primeEndOfRead", "Base", "nucleotideFreq")
	ptfum <- ggplot(tfum, aes(x=distanceFrom5primeEndOfRead, y=nucleotideFreq, group=Base, color=Base))
	ptfum + geom_line() + ylim(0, .4) + xlab("Distance from 5' end of Reads") + ylab("Nucleotide Frequency") + ggtitle("Genomic Nucleotide Frequency relative to read positions (Unique)")
	ggsave(file=paste(args[2],"tfu.png",sep="/"))
}

## Fragment GC% Distribution
if ( file.exists( paste(args[1],"tagGCcontent.txt",sep="/") ) && file.exists( paste(args[1],"genomeGCcontent.txt",sep="/") ) ) {
	tgcc <- read.delim(file=paste(args[1],"tagGCcontent.txt",sep="/"), header=T, sep="\t")
	ggcc <- read.delim(file=paste(args[1],"genomeGCcontent.txt",sep="/"), header=T, sep="\t")
	gcc <- data.frame(tgcc[,c(1,3)], ggcc[,3])
	colnames(gcc) <- c("GC", "tagFraction", "genomeFraction")
	gccm <- melt(gcc, id.var="GC")
	colnames(gccm) <- c("GcContentOfFragments", "variable", "normalizedFraction")
	pgccm <- ggplot(gccm, aes(x=GcContentOfFragments, y=normalizedFraction, group=variable, color=variable))
	pgccm + geom_line() + scale_x_continuous(labels = percent) + xlab("GC-content of fragments") + ylab("Normalized Fraction") + ggtitle("Fragment GC% Distribution") + 
	scale_colour_hue(name="Source", breaks=c("genomeFraction","tagFraction"), labels=c("Genome", "ChIP-Seq"))
	ggsave(file=paste(args[2],"gcc.png",sep="/"))
}


