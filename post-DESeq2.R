#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
	  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      post-DESeq2.R - Generate tablse and plots from DESeq2 results of \"run-DESeq2.R\" script
 
      Arguments:
      --in=someValue   		 - input directory path (i.e. containing \"run-DESeq2.RData\" file from run-DESeq2.R script)
      --out=someValue   	 - output directory path
      --lfc=someValue   	 - Log2FC threshold (Default: 1)
      --qvalue=someValue	 - q-value threshold (Default: 0.01)
      --emapper=someValue	 - emapper annotations report file
      --gene_trans_map=someValue - gene_trans_map (gene <=> transcript)
      --bygene=1                 - annotation by gene ? (Default: No)
      --help            	 - print this Help
     
      Example:
      ./post-DESeq2.R --in=\"./output_from_run_DESeq2\" --out=\"./output\" [ --lfc=1 --qvalue=0.01 --emapper=./emapper.annotations --gene_trans_map=./transcripts.gene_trans_map --bygene=1 ]
      
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
	cat("\nERROR: Missing input directory !\n\n")
	sink()
	q(save="no")
}
if(is.null(argsL[['out']])) {
	sink(stderr())
	cat("\nERROR: Missing output dir!\n\n")
	sink()
	q(save="no")
}

if ( ! dir.exists(argsL[['in']]) ) {
	sink(stderr())
	cat("\nERROR: Input directory not found!\n\n")
	sink()
	q(save="no")
}

if ( ! file.exists( paste(argsL[['in']], "/run-DESeq2.RData",sep="") ) ) {
	sink(stderr())
	cat("\nERROR: File run-DESeq2.RData not found!\n\n")
	sink()
	q(save="no")
}

if ( (! file.exists(argsL[['out']]))|(! file.info(argsL[['out']])[1,"isdir"] ) ) {
	sink(stderr())
	cat("\nERROR: Output directory not found!\n\n")
	sink()
	q(save="no")
}

if ( is.null( argsL[['lfc']] ) ) {
	sink(stderr())
	argsL[['lfc']]<-1
	cat(paste("\n","WARNING: Using default Log2-Fold-Change threshold (",argsL[['lfc']],") !\n\n",sep="") )
	sink()
} else {
	argsL[['lfc']] <- as.integer( argsL[['lfc']] )
}

if ( is.null( argsL[['qvalue']] ) ) {
	sink(stderr())
	argsL[['qvalue']]<-0.01
	cat(paste("\n","WARNING: Using default q-value threshold (",argsL[['qvalue']],") !\n\n",sep="") )
	sink()
} else {
	argsL[['qvalue']] <- as.double( argsL[['qvalue']] )
}

if ( argsL[['lfc']] < 0 ) {
	sink(stderr())
	cat(paste("\n","ERROR: It's not possible to use a Log2FC lower than zero (",argsL[['lfc']],") !\n\n",sep="") )
	sink()
	q(save="no")
}

if ( ! ((as.double(argsL[['qvalue']]) >=0 ) & (as.double(argsL[['qvalue']]) <= 1)) ) {
	sink(stderr())
	cat(paste("\n","ERROR: It's not possible to use a qvalue outside the range 0..1 (",argsL[['qvalue']],") !\n\n",sep="") )
	sink()
	q(save="no")
}

suppressPackageStartupMessages( library("methods") )
suppressPackageStartupMessages( library("DESeq2") )
suppressPackageStartupMessages( library("pheatmap") )
suppressPackageStartupMessages( library("vsn") )
suppressPackageStartupMessages( library("RColorBrewer") )
suppressPackageStartupMessages( library("ggplot2") )
suppressPackageStartupMessages( library("ggrepel") )
suppressPackageStartupMessages( library("xlsx") )
suppressPackageStartupMessages( library('topGO') )
suppressPackageStartupMessages( library('scales') )


plotEnrichment<-function(ggdata,method,ont,k) {
	ont.desc <- list()
	ont.desc[["BP"]] <- "Biological Process"
	ont.desc[["CC"]] <- "Celular Component"
	ont.desc[["MF"]] <- "Molecular Function"
	k.desc <- list()
	k.desc[['de']] <- "Differentially expressed genes"
	k.desc[['up']] <- "Up-regulated genes"
	k.desc[['down']] <- "Down-regulated genes"

	colnames(ggdata) <- c("GO.ID","Term","Significant", "Statistic")
	ggdata$Statistic <- as.numeric(ggdata$Statistic)
	ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
	gg1 <- ggplot(ggdata,
			aes(x = Term, y = -log10(Statistic), size = Significant, fill = -log10(Statistic))) +

			expand_limits(y = 1) +
			geom_point(shape = 21) +
			scale_size(range = c(2.5,12.5)) +
			scale_fill_continuous(low = 'royalblue', high = 'red4') +

			xlab('') + ylab('Enrichment score') +
			labs(
			     title = paste('GO ',ont.desc[[ont]], sep=""),
			     subtitle = paste('Top ',dim(ggdata)[1],' terms in ',k.desc[[k]],' ordered by ',method, 'statistic', sep=""),
			     caption = 'Cut-off lines drawn at equivalents of q-value=0.2, q-value=0.1, q-value=0.05') +

			geom_hline(yintercept = c(-log10(0.2), -log10(0.1), -log10(0.05)),
				   linetype = c("dotted", "dotted", "dotted"),
				   colour = c("black", "black", "black"),
				   size = c(0.5, 0.5, 0.5)) +

			theme_bw(base_size = 24) +
			theme(
			      legend.position = 'right',
			      legend.background = element_rect(),
			      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
			      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
			      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

			      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
			      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
			      axis.title = element_text(size = 12, face = 'bold'),
			      axis.title.x = element_text(size = 12, face = 'bold'),
			      axis.title.y = element_text(size = 12, face = 'bold'),
			      axis.line = element_line(colour = 'black'),

			      #Legend
			      legend.key = element_blank(), # removes the border
			      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
			      legend.text = element_text(size = 14, face = "bold"), # Text size
			      title = element_text(size = 14, face = "bold")) +

			coord_flip()
	return(gg1)
}

new.in=argsL[['in']]
new.out=argsL[['out']]

threshold.qvalue <- as.double(argsL[['qvalue']])
threshold.lfc <- as.double(argsL[['lfc']])

emapper.file <- argsL[['emapper']]
gene_trans_map.file <- argsL[['gene_trans_map']]
bygene <- argsL[['bygene']]

load( paste(new.in, "/run-DESeq2.RData",sep="") )

options(Java.parameters = "-Xmx8000m")

if(is.null(emapper.file)) {
	sink(stderr())
	cat("\nWARNING: Missing emapper annotations report file !\n\n")
	sink()
	q(save="no")
} else {

	if ( ! file.exists(emapper.file) ) {
		sink(stderr())
		cat("\nERROR: Not found emapper annotations report file !\n\n")
		sink()
		q(save="no")
	}

	if(is.null(gene_trans_map.file)) {
		sink(stderr())
		cat("\nWARNING: Missing gene_trans_map file !\n\n")
		sink()
		q(save="no")
	}


	eggnog.df<-read.delim(file=emapper.file,sep="\t",header=TRUE)
	colnames(eggnog.df)

	genetransmap.df<- read.delim(file=gene_trans_map.file,sep="\t",header=FALSE)

	colnames(genetransmap.df) <- c("gene","transcript")
	rownames(genetransmap.df) <- genetransmap.df$transcript

	eggnog.df$query_id <- gsub("\\.p\\d+$",'',eggnog.df$query_id,perl=TRUE)
	if ( length(eggnog.df$query_id) != length(unique(eggnog.df$query_id)) ) {
		sink(stderr())
		cat("\nERROR: Not unique query_ids in emapper annotaions report file!\n\n")
		sink()
		q(save="no")
	}
	rownames(eggnog.df) <- eggnog.df$query_id

	eggnog.df$gene_id <- genetransmap.df[eggnog.df$query_id, 'gene']

	if ((! is.null(bygene) )&( as.integer(bygene) != 0 )) {

		annotation.df <- data.frame(row.names=unique(eggnog.df$gene_id))
		for (q in unique(eggnog.df$gene_id) ) {
			eggnog.df.temp <- subset(eggnog.df, gene_id==q)
			eggnog.df.temp.sel.score <- subset(eggnog.df.temp, score==max(eggnog.df.temp$score))
			eggnog.df.temp.sel.evalue <- subset(eggnog.df.temp.sel.score, evalue==min(eggnog.df.temp.sel.score$evalue))[1,]
			annotation.df <- rbind(annotation.df, eggnog.df.temp.sel.evalue)
		}
		annotation.df$X <- annotation.df$gene_id
	} else {
		annotation.df <- eggnog.df
		annotation.df$X <- eggnog.df$query_id
	}

	write.table(annotation.df[,c('gene_id','GOs')], 
	            file=paste(new.out,"/gene2GOs.txt",sep=""),
        	    col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

	gene2GOs <- readMappings(file=paste(new.out,"/gene2GOs.txt",sep=""))

	gene2GOs[sapply(sapply(gene2GOs, function(x){ x=="-" }), any)] <- NA

	for (er in names(res) ) {

		print(er)

		rm(dgeres.df)
		rm(dge.res)

		if ( is.null(argsL[['lfcshrink']]) ) {
			dgeres.df <- as.data.frame(  res[[ er ]]  )
			dgeres.df[[idcol]] <- rownames(dgeres.df)

			rownames(dgeres.df) <- dgeres.df[[idcol]]
			colnames(dgeres.df) <- c("baseMean","logFC","lfcSE","stat","PValue", "FDR", idcol)
		} else {
			dgeres.df <- as.data.frame(  res[[ er ]]  )
			dgeres.df[[idcol]] <- rownames(dgeres.df)

			rownames(dgeres.df) <- dgeres.df[[idcol]]
			colnames(dgeres.df) <- c("baseMean","logFC","lfcSE","PValue", "FDR", idcol)
		}

		dge.res <- (merge( dgeres.df, deseq_ncounts, by=idcol))[,c(idcol, g$name, mean.cols, 'logFC', 'lfcSE', 'PValue', 'FDR')]
		rownames(dge.res) <- dge.res[[idcol]]

		write.table(file=paste(new.out,"/",gsub("-","x",er),"_eggNOG.txt",sep=""), 
			    x=merge(x=dge.res, y=annotation.df, by.x=idcol, by.y='X', all.x=TRUE),
			    row.names=FALSE, 
			    quote=FALSE, 
			    na="",
			    sep="\t")

		
		write.xlsx(merge(x=dge.res, y=annotation.df, by.x=idcol, by.y='X', all.x=TRUE), file=paste(new.out, "/",gsub("-","x",er),"_eggNOG.xlsx", sep=""),
			   sheetName = gsub("-","x",er),
			   col.names = TRUE,
			   row.names = FALSE,
			   append = FALSE
			   )
		gc()

		
		geneList <- abs(as.numeric(as.data.frame(res[[er]])$log2FoldChange))
		names(geneList) <- rownames(as.data.frame(res[[er]]))

		topde <<- apply(as.data.frame(res[[er]]), 1, function(x,q,f) { return( ((!is.na(x['padj']))&(x['padj']<=q))&(abs(x['log2FoldChange'])>=f) ) }, threshold.qvalue, threshold.lfc)
		topde.up <<- apply(as.data.frame(res[[er]]), 1, function(x,q,f) { return( ((!is.na(x['padj']))&(x['padj']<=q))&(x['log2FoldChange']>=f) ) }, threshold.qvalue, threshold.lfc)
		topde.down <<- apply(as.data.frame(res[[er]]), 1, function(x,q,f) { return( ((!is.na(x['padj']))&(x['padj']<=q))&(x['log2FoldChange']<=-f) ) }, threshold.qvalue, threshold.lfc)

		geneSelFun <- list()
		geneSelFun[['de']] <- function(x) { return( topde ) }
		geneSelFun[['up']] <- function(x) { return( topde.up ) }
		geneSelFun[['down']] <- function(x) { return( topde.down ) }

		allRes<-list()

		for (onto in c("BP", "MF", "CC")) {

			allRes[[onto]] <- list()

			for (k in c("de", "up", "down")) {
				GOdata <- new("topGOdata",
					      description = er, ontology = onto,
					      allGenes = geneList, geneSelectionFun=geneSelFun[[k]],
					      nodeSize = 10,
					      annot = annFUN.gene2GO, gene2GO = gene2GOs)

				#sigGenes(GOdata[[onto]])
				weight_fisher_result<-runTest(GOdata, algorithm='weight01', statistic='fisher', scoreOrder="decreasing")
				weight_ks_result<-runTest(GOdata, algorithm='weight01', statistic='ks', scoreOrder="decreasing")
				allGO<-usedGO(GOdata)
				all_res.tmp <- GenTable(GOdata, weightFisher=weight_fisher_result, weightKS=weight_ks_result, orderBy='weightKS', topNodes=length(allGO))
				all_res <- subset(all_res.tmp, weightKS < 1)
				all_res$p.adj <- p.adjust(all_res$weightKS, method="BH")
				allRes[[onto]][[k]] <- all_res[order(all_res$p.adj),]

				rm(GOdata)
				rm(all_res)
				rm(weight_fisher_result)
				rm(all_res.tmp)
				rm(weight_ks_result)

				goEnrichment<-allRes[[onto]][[k]][1:10,c("GO.ID","Term","Significant","weightKS")]
				svg( paste(new.out, "/GOEnrichment_",gsub("-","x",er),"_",onto,"_",k,".svg", sep=""), width = 10, height = 10)
				plotEnrichment(goEnrichment,"weighted Kolmogorov-Smirnov",onto,k)
				graphics.off()
		
				write.table(file=paste(new.out,"/GOEnrichment_",gsub("-","x",er),"_",onto,"_",k,".txt",sep=""), 
					    x=allRes[[onto]][[k]],
					    row.names=FALSE, 
					    quote=FALSE, 
					    na="",
					    sep="\t")


				write.xlsx(allRes[[onto]][[k]],
					   file=paste(new.out,"/GOEnrichment_",gsub("-","x",er),"_",onto,"_",k,".xlsx",sep=""),
					   sheetName = paste("GOEA_",gsub("-","x",er),"_",onto,"_",k,sep=""),
					   col.names = TRUE,
					   row.names = FALSE,
					   append = FALSE
					   )
			}
		}

		rm(allRes)
		gc()

	}
	

}



# DESeq2 Library Normalization
# https://www.youtube.com/watch?v=UFB993xufUU

# Questão básica (p-value)
# Cada gene passará por um teste de hipótese e cada teste resultará em um nível descritivo de significância estatística ou p-valor (p-value)
# Avaliar se a diferença observada pode ter aparecido pelo acaso ("random chance"), ou seja, qual o nível de confiança no teste, naquele teste 
# em particular
# Ou seja, se utilizarmos um limiar máximo de 0,05 para o p-valor, significa que ao 
# realizarmos o mesmo experimento com o mesmo gene uma grande quantidade de vezes, em 5% desses experimentos
# resultariam em uma falha na decisão ( rejeitar a hipótese nula = neste caso dos genes de que não há diferença de expressão )
# ou seja, erro do tipo I, quando a hipótese nula é verdadeira e você a rejeita (falso positivo = falsa diferença de expressão).
# * MUITO IMPORTANTE: NÃO QUER DIZER QUE USAR UM LIMIAR DE P-VALOR MENOR QUE 0,05, PARA SELECIONAR GENES 
# ISSO RESULTARÁ EM 5% DE RESULTADOS FALSOS POSITIVOS DENTRE ESSES GENES SELECIONADOS (VEJA OS VÍDEOS SOBRE FDR)
# O método de FDR serve para ajustar os valores de p (p-values) para que a seleção seja feita em cima desses e neste caso
# considerar que usando o limiar de 5% para o valor de p ajustado, poderemos esperar essa porcentagem de falsos positivos
# dentre os resultados.
# * OUTRA COISA: o valor de p não indica o QUANTO há de diferença, quem lhe dá isso é o log2FC (razão entre os valores de expressão)
#
# No exemplo do vídeo a seguir efeito de medicamentos, em nosso caso o efeito do metal
# https://www.youtube.com/watch?v=vemZtEM63GY
#
# Veja sobre testes de hipótese (caso necessário assista antes de ver os vídeos sobre valores de p):
# https://www.youtube.com/watch?v=0oc49DyA3hU
# https://www.youtube.com/watch?v=5koKb5B_YWo
# 
# E mais sobre os valores de p (p-values) - probabilidades de ter obtido aquele resultado ao acaso, ou resultado igual ou raro
# https://www.youtube.com/watch?v=5Z9OIYA8He8
# 

# A função DESeq calcula, para cada gene e para cada amostra, um teste diagnóstico para amostras discrepantes ("outliers")
# chamado distância de Cook. A distância de Cook é uma medida de quanto uma única amostra está influenciando os 
# coeficientes ajustados para um gene, e um grande valor da distância de Cook destina-se a indicar uma 
# contagem de valores discrepantes. 
# As distâncias de Cook são armazenadas como uma matriz disponível em assays(dds)[["cooks"]].

svg( paste(new.out, "/Boxplot_Cooks.svg",sep=""), width = 20, height = 10)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
graphics.off()


# A estimativa de dispersão para cada gene poderia ser feita dentro de cada condição, caso tivéssemos centenas de réplicas.
# Normalmente por questões de viabilidade econômica, é comum os experimentos conterem apenas 3 réplicas. 
# Para contornar o problema, DESeq2 considera a suposição de que 
# genes com expressão similar possuem dispersão também similar
# então considera a expressão de outros genes para estimar uma dispersão comum a genes com expressão similar

# O gráfico de estimativa de dispersão é uma ferramenta de diagnóstico. 
# O gráfico exibe o encolhimento ("shrunk") da dispersão para uma estimativa final a partir do observado 
# em direção à curva ajustada. A quantidade do encolhimento 
# depende do tamanho da amostra, do número de coeficientes, das médias dos valores de expressão e 
# da variabilidade das estimativas obtidas

svg( paste(new.out, "/DispEst.svg",sep=""), width = 10, height = 10)
plotDispEsts(dds)
graphics.off()


# O DESeq2 executa uma filtragem independente do teste estatístico usando a média das contagens normalizadas.
# Um limiar desse filtro é encontrado, sendo aquele que otimiza o número de p-valores ajustados menor que um nível de significância alfa. 

# Um pouco sobre o ajuste do p-value e sobre a Filtragem Independente
# https://www.youtube.com/watch?v=Gi0JdrxRq5s&t=738s

# Mais sobre FDR
# https://www.youtube.com/watch?v=K8LQSvtjcEo

# O objetivo desta filtragem é descartar os genes que possuem pouca chance de apresentar evidência significativa
# sem mesmo olhar para o teste estatístico. Isso tem o efeito de aumentar o poder de detecção de erros do tipo I (FDR).

# Os valores de p-ajustados para os genes que não passam no filtro aparecem como NA.

# A visualização da otimização pode ser observada no gráfico a seguir.
# A função maximiza o número de rejeições (p valores ajustados menores que um nível alfa), 
# sobre os quantis da média de contagens normalizadas. 

# O limiar encontrado (linha vertical) é o menor quantil da filtragem onde o número de rejeições
# está dentro de 1 resíduo de desvio padrão do pico da curva ajustada aos números de rejeições
# sobre os quantis da filtragem.



# Para cada análise comparativa em relação à amostra controle
deg<-c()
for (er in names(res)) {

  svg( paste(new.out, "/IndependentFiltering_",gsub("-","x",er),".svg", sep=""), width = 10, height = 10)
  metadata(res[[er]])$alpha
  metadata(res[[er]])$filterThreshold
  plot(metadata(res[[er]])$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res[[er]])$lo.fit, col="red")
  abline(v=metadata(res[[er]])$filterTheta)
  graphics.off()
  
  svg(paste(new.out, "/MAPlot_",gsub("\\.","_",threshold.qvalue),"_",gsub("-","x",er),".svg", sep=""), width = 10, height = 10)
  plotMA(res[[er]],alpha=threshold.qvalue)
  graphics.off()
  
  selected.res <- subset(as.data.frame(res[[er]]), (padj<=threshold.qvalue)&(abs(log2FoldChange)>=threshold.lfc) )

  deg<-unique(c(deg,rownames(selected.res)))
  
  if (length(deg) > 0) { 
	  svg(paste(new.out, "/VolcanoPlot_",threshold.lfc,"_",gsub("\\.","_",threshold.qvalue),"_",gsub("-","x",er),".svg", sep=""), width = 10, height = 10)

	  genes <- as.data.frame(res[[er]])
	
	  genes$Significant <- ifelse( ( (genes$padj <= threshold.qvalue) & (genes$log2FoldChange >= threshold.lfc) ), 
				      "Up-regulated genes", 
				      ifelse( ( (genes$padj <= threshold.qvalue) & (genes$log2FoldChange <= -threshold.lfc) ), 
					     "Down-regulated genes", 
					     "Not DEG"
					    ) 
				     )

	  genes$padj[is.na(genes$padj)] <- 1
	  ggplot(genes, aes(x = log2FoldChange, y = -log10(padj)))+
	    labs(title=gsub("-"," x ", er))+
	    geom_point(aes(color = Significant), size=1) +
	    scale_color_manual(values = c("Up-regulated genes"="red", "Down-regulated genes"="green", "Not DEG"="grey")) +
	    theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5) )+
	    geom_vline(xintercept=c(-2,2),linetype='dotted')+
	    geom_hline(yintercept=-log10(threshold.qvalue), linetype='dotted' )+
	    annotate("text", x=as.integer(min(genes$log2FoldChange)), y = -log10(threshold.qvalue), label = paste("adj. p-value \u2264 ",threshold.qvalue,sep=""), vjust=-0.5)+
	    annotate("text", x=-2, y=as.integer(max(-log10(genes$padj))-1), label = paste("log2FC \u2264 -",threshold.lfc,sep=""), hjust=1.1)+
	    annotate("text", x=2, y=as.integer(max(-log10(genes$padj))-1), label = paste("log2FC \u2265 ",threshold.lfc,sep=""), hjust=-0.1)
	    
	   graphics.off()
  }
  # http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
}
rm(er)

rld <- rlog(dds, blind=FALSE)
# this gives log2(n + 1)
ntd <- normTransform(dds)
svg(paste(new.out, "/meanSdPlot_ntd",".svg", sep=""), width = 10, height = 10)
meanSdPlot(assay(ntd))
graphics.off()
svg(paste(new.out, "/meanSdPlot_rld",".svg", sep=""), width = 10, height = 10)
meanSdPlot(assay(rld))
graphics.off()

df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- colData(dds)[,c("id")]
colnames(df) <- 'group'


if (length(deg) > 0) { 
	svg(paste(new.out, "/HeatMap_rld_",threshold.lfc,"_",gsub("\\.","_",threshold.qvalue),".svg", sep=""), width = 10, height = 10)
	pheatmap(assay(rld)[deg,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE,
		 scale="row",
		 annotation_col=df)
	graphics.off()
}

rld.df <- as.data.frame(assay(rld))
rld.df$ID <- rownames(rld.df)
write.table(rld.df[,c('ID',samps)], 
	    file=paste(new.out, "/FullRLDnormData.txt", sep=""),
	    col.names=TRUE, 
	    row.names=FALSE, 
	    quote=FALSE, 
	    sep="\t")


write.xlsx(rld.df[,c('ID',samps)], file=paste(new.out, "/FullRLDnormData.xlsx", sep=""),
           sheetName = "RLD",
           col.names = TRUE,
           row.names = FALSE,
           append = FALSE
           )

svg(paste(new.out, "/HeatMap_Sample2Sample_rld",".svg", sep=""), width = 10, height = 10)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$group, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
graphics.off()

svg(paste(new.out, "/PCA_rld",".svg", sep=""), width = 10, height = 10)
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
geom_text_repel(aes(label = rownames(pcaData)),
                 show.legend=FALSE,
                 size = 3.5)
graphics.off()




# GLM
# Como a expressão do gene pode ser determinada por uma determina da condição 

# Linear regression
# https://www.youtube.com/watch?v=nk2CQITm_eo&t=63s
#
# Linear Regression and Linear Models
# https://www.youtube.com/watch?v=PaFPbb66DxQ
#
# General Linear Models
# https://www.youtube.com/watch?v=NF5_btOaCig
# * Os exemplos aqui são parao teste-t usando a média para estimar os interceptos
# pois no caso do teste-t a distribuição considerara é a normal, 
# no caso do DESeq2, a distribuição considerada para modelar os dados é uma distribuição
# binomial negativa, então não é exatamente o mesmo, mas serve para entender o conceito.
# Quando ele mostra os 0's e 1's isso é possível observar aqui armazenado na variável do R "design"
# veja a seguir a matriz de design:
#design
# quando ele fala sobre o "overall mean", é a nossa coluna "baseMean" no resultado final

# Este outro vídeo (o inglês com sotaque aqui é difícil) é de um outro canal sobre o DESeq2
# https://www.youtube.com/watch?v=mq6UvDneKc0
# mostra por volta do minuto 14 o encolhimento da estimativa de dispersão para a modelagem da
# distribuição binomial negativa, baseada em 2 parâmetros (média, dispersão). O teste estatístico
# note que não é o teste-t (como no vídeo anterior sobre GLM) mas o Teste de Wald, onde a hipótese
# nula é que o log2FC == 0 (ou seja, sem diferença de expressão)

# No link abaixo tem um tutorial (não é necessário fazer, apenas compreender)
# https://uclouvain-cbio.github.io/WSBIM2122/sec-rnaseq.html

