#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      mkGOplot.R - Run DESeq2 method
 
      Arguments:
      --deseq2=someValue   	- DESeq2 analyis file path
      --goea=someValue		- goatools enrichment analysis file path
      --gene2go=someValue	- gene-to-GO file path without header line (required for shrinkage) [Default: disabled shrinkage]
      --out=someValue   	- SVG and TXT output file path (without extension .svg or .txt)
      --threshold=someValue	- adjusted p-value threshold [Default: 1e-5]
      --shrinkage=someValue	- reduce terms using Hierarchical Clustering height cutoff [Default: 1 - no reduction]
      --help            	- print this Help

      Example:
      ./mkGOplot.R --deseq2=\"deseq2_analyis.txt\" --goea=\"goatools_analysis.txt\" --out=\"./output/goplot\" 
      
      ./mkGOplot.R --deseq2=\"goea.PER.PEA.both.tsv\" --goea=\"PER-PEA.both.txt\" --out=\"PER-PEA.both\" 

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
	

if(is.null(argsL[['deseq2']])) {
	sink(stderr())
	cat("\nERROR: Missing deseq2 analysis file!\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['out']])) {
	sink(stderr())
	cat("\nERROR: Missing SVG output file path!\n\n")
	sink()
	q(save="no")
} else {
	argsL[['out']] <- tools::file_path_sans_ext( argsL[['out']] )
}

if(is.null(argsL[['shrinkage']])) {
	argsL[['shrinkage']] <- 1
} else {
	if(is.null(argsL[['gene2go']])) {
		sink(stderr())
		cat("\nERROR: Missing gene2go file path!\n\n")
		sink()
		q(save="no")
	} else {
		if ( ! file.exists(argsL[['gene2go']]) ) {
			sink(stderr())
			cat("\nERROR: gene2go analysis file not found!\n\n")
			sink()
			q(save="no")
		}
	}
}

if(is.null(argsL[['goea']])) {
	sink(stderr())
	cat("\nERROR: Missing GO Enrichment Analysis file!\n\n")
	sink()
	q(save="no")
}

if ( ! file.exists(argsL[['deseq2']]) ) {
	sink(stderr())
	cat("\nERROR: DESeq2 analysis file not found!\n\n")
	sink()
	q(save="no")
}

if ( ! file.exists(argsL[['goea']]) ) {
	sink(stderr())
	cat("\nERROR: GO Enrichment Analysis file not found!\n\n")
	sink()
	q(save="no")
}

if(is.null(argsL[['threshold']])) {
	argsL[['threshold']] = 1e-5	
}

suppressPackageStartupMessages( library('GOplot') )

my.GOBubble <- function (data, display, title, colour, labels, ID = T, table.legend = T, 
                         table.col = T, bg.col = F, go.text.size = 1, go.text.threshold = 1.3) 
{
  zscore <- adj_pval <- category <- count <- id <- term <- NULL
  if (missing(display)) 
    display <- "single"
  if (missing(title)) 
    title <- ""
  if (missing(colour)) 
    cols <- c("chartreuse4", "brown2", "cornflowerblue")
  else cols <- colour
  if (missing(labels)) 
    labels <- 5
  if (bg.col == T & display == "single") 
    cat("Parameter bg.col will be ignored. To use the parameter change display to 'multiple'")
  colnames(data) <- tolower(colnames(data))
  if (!"count" %in% colnames(data)) {
    rang <- c(5, 5)
    data$count <- rep(1, dim(data)[1])
  }
  else {
    rang <- c(1, 30)
  }
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  zlim=round(max( c(abs(min(sub$zscore)), abs(max(sub$zscore)) )  ))+1
  g <- ggplot(sub, aes(zscore, adj_pval, fill = category, size = count)) + xlim(-zlim,zlim) +
    labs(title = title, x = "z-score", y = "-log (adj p-value)") + 
    geom_point(shape = 21, col = "black", alpha = 1/2) + 
    geom_hline(yintercept = go.text.threshold, col = "orange") + scale_size(range = rang, 
                                                              guide = "none")
  if (!is.character(labels)) 
    sub2 <- subset(sub, subset = sub$adj_pval >= labels)
  else sub2 <- subset(sub, sub$id %in% labels | sub$term %in% 
                        labels)
  if (display == "single") {
    g <- g + scale_fill_manual("Category", values = cols, 
                               labels = c("Biological Process", "Cellular Component", 
                                          "Molecular Function")) + theme(legend.position = "bottom") + 
      annotate("text", x = min(sub$zscore) + 0.2, y = go.text.threshold+0.1, 
               label = "Threshold", colour = "orange", size = 4)
    if (ID) 
      g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, 
                                          label = id), size = go.text.size)
    else g <- g + geom_text(data = sub2, aes(x = zscore, 
                                             y = adj_pval, label = term), size = go.text.size)
    if (table.legend) {
      if (table.col) 
        table <- GOplot:::draw_table(sub2, col = cols)
      else table <- GOplot:::draw_table(sub2)
      g <- g + theme(axis.text = element_text(size = 14), 
                     axis.line = element_line(colour = "grey80"), 
                     axis.ticks = element_line(colour = "grey80"), 
                     axis.title = element_text(size = 14, face = "bold"), 
                     panel.background = element_blank(), panel.grid.minor = element_blank(), 
                     panel.grid.major = element_line(colour = "grey80"), 
                     plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, table, ncol = 2)
    }
    else {
      g + theme(axis.text = element_text(size = 14), axis.line = element_line(colour = "grey80"), 
                axis.ticks = element_line(colour = "grey80"), 
                axis.title = element_text(size = 14, face = "bold"), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_line(colour = "grey80"), 
                plot.background = element_blank())
    }
  }
  else {
    if (bg.col) {
      dummy_col <- data.frame(category = c("BP", "CC", 
                                           "MF"), adj_pval = sub$adj_pval[1:3], zscore = sub$zscore[1:3], 
                              size = 1:3, count = 1:3)
      g <- g + geom_rect(data = dummy_col, aes(fill = category), 
                         xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                         alpha = 0.1) + facet_grid(. ~ category, space = "free_x", 
                                                   scales = "free_x") + scale_fill_manual(values = cols, 
                                                                                          guide = "none")
    }
    else {
      g <- g + facet_grid(. ~ category, space = "free_x", 
                          scales = "free_x") + scale_fill_manual(values = cols, 
                                                                 guide = "none")
    }
    if (ID) {
      g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, 
                                     label = id), size = go.text.size) + theme(axis.title = element_text(size = 14, 
                                                                                              face = "bold"), axis.text = element_text(size = 14), 
                                                                    axis.line = element_line(colour = "grey80"), 
                                                                    axis.ticks = element_line(colour = "grey80"), 
                                                                    panel.border = element_rect(fill = "transparent", 
                                                                                                colour = "grey80"), panel.background = element_blank(), 
                                                                    panel.grid = element_blank(), plot.background = element_blank())
    }
    else {
      g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, 
                                     label = term), size = go.text.size) + theme(axis.title = element_text(size = 14, 
                                                                                                face = "bold"), axis.text = element_text(size = 14), 
                                                                      axis.line = element_line(colour = "grey80"), 
                                                                      axis.ticks = element_line(colour = "grey80"), 
                                                                      panel.border = element_rect(fill = "transparent", 
                                                                                                  colour = "grey80"), panel.background = element_blank(), 
                                                                      panel.grid = element_blank(), plot.background = element_blank())
    }
  }
}


# PARÂMETROS
goatools.path <- file.path(argsL[['goea']])
deseq2.path <- file.path(argsL[['deseq2']])

## CÓDIGO

goatools.df <- read.delim(goatools.path,sep="\t")
deseq2.df <- read.delim(deseq2.path, sep="\t")

bubble.list <- list()

bubble.list$goatools <- subset( goatools.df[,c('NS', 'X..GO', 'name', 'study_items', 'p_fdr_bh')], p_fdr_bh <= 1)
colnames(bubble.list$goatools) <- c("Category", "ID", "Term", "Genes", "adj_pval")


bubble.list$genelist <- deseq2.df[,c('X', 'logFC', 'FDR')]
colnames(bubble.list$genelist) <- c("ID", "logFC", "adj.P.Val")

if ((dim(bubble.list$goatools)[1] > 0)&&(dim(bubble.list$genelist)[1]>0)) {
  dim(bubble.list$goatools)
 
  if (argsL[['shrinkage']]<1) {
	suppressPackageStartupMessages( require('GO.db') )
	suppressPackageStartupMessages( require('GOSemSim') )
	
        ontology <- c("MF", "BP", "CC")
        gss.GO <- list()
	
	gene2go.df <- read.delim(file=file.path( argsL[['gene2go']] ),
                                 sep="\t",
                                 comment.char = "#",
                                 header=FALSE
                                )
	gene2go.df$V3 <- 'IEA'
	colnames(gene2go.df) <- c('GID','GO','EVIDENCE')
	
	tbl=toTable(GOTERM)
	gene2go.go.df <- subset(gene2go.df,  GO != "")
	goAnno <- subset(gene2go.go.df,  GO %in% tbl$go_id)[,c('GID','GO')]
	goAnno$ONTOLOGY <- suppressMessages( select(GO.db, goAnno$GO, c("ONTOLOGY"))["ONTOLOGY"] )

        for (ont in ontology) {
		goAnnoOnt <- subset(goAnno, ONTOLOGY == ont)
		gss.data <- new(     "GOSemSimDATA", 
        	      	             keys = as.character(goAnnoOnt$GID), 
                	             ont = ont, 
	              	             geneAnno = goAnnoOnt
                	       )
		
		gss.data@IC <-  GOSemSim:::computeIC(subset(goAnno, ONTOLOGY == ont), ont)
		
		sel.df <- subset(bubble.list$goatools, Category == ont)
		sel <- sel.df$adj_pval
		names(sel) <- sel.df$ID
		
		d <- as.dist(1-mgoSim(names(sel), 
                                      names(sel), 
                                      semData=gss.data, 
                                      measure="Wang", 
                                      combine=NULL
                                     )
                            )
		
		save(list = ls(all=TRUE,envir=globalenv()), file = "plot.RData")

		hc <- hclust(d, method="average")
		hc$height <- round(hc$height, 10)
		cls <- cutree(hc, h=argsL[['shrinkage']] )

		for (cl in unique(cls)) {
			gsel <- sel[names(cls[cls==cl])]
			gsel.min <- names(gsel[gsel==min(gsel)])
        		for (m in gsel.min) {
				gss.GO[[m]] <- list(cluster=cl, adj_pval=min(gsel))
			}
		}
		
	}
	
	bubble.list$goatools <- subset(bubble.list$goatools, ID %in% names(gss.GO))
  }
  
  #Building the circ object
  circ <- circle_dat(bubble.list$goatools, bubble.list$genelist)
 
  #Faceting the plot
  # Reduce redundant terms with a gene overlap >= 0.75...
  #reduced_circ <- reduce_overlap(circ, overlap = 0.75)
  # ...and plot it
  my.GOBubble(circ, 
           display = 'multiple', 
           table.legend=FALSE, 
           ID=FALSE, 
           table.col=FALSE,
           title=basename(argsL[['out']]),
           go.text.size=1,
           go.text.threshold =-log10(argsL[['threshold']]),
           labels=-log10(argsL[['threshold']]), 
           bg.col=TRUE)+
    ggsave(
      filename=paste(argsL[['out']],'.svg',sep=""),
      device="svg",
      width=18,
      height=10,
      units="cm"
    )
   
    write.table(file=paste(argsL[['out']],'.txt',sep=""),
		x=subset(circ[!duplicated(circ$term), ], adj_pval <= argsL[['threshold']])[, c('category','ID','term','adj_pval','zscore')],
		row.names=FALSE, quote=FALSE, sep="\t")
}

