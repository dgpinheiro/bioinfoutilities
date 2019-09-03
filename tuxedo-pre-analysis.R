#!/usr/bin/env Rscript


# Get command-line arguments and validate them

args <- commandArgs(trailingOnly = TRUE)

if (is.na(args[1])) {
        sink(stderr())
        cat("ERROR: Missing input directory (cufflinks directories {cuffdiff & cuffmerge}) !\n")
        sink()
        quit(save="no", status=10);
} else {
        if (!file.exists(args[1])){
                sink(stderr())
                cat(paste(c("ERROR: Input directory (",args[1],") not exists !\n"),collapse=""))
                sink()
                quit(save="no", status=10)
        }
}
input.dir <- args[1]

if (is.na(args[2])) {
        sink(stderr())
        cat("ERROR: Missing output directory (CummeRbund) !\n")
        sink()
        quit(save="no", status=10)
} else {
        if (!file.exists(args[2])){
                sink(stderr())
                cat(paste(c("Creating output directory (",args[2],") !\n"),collapse=""))
                sink()
		dir.create(args[2], showWarnings = TRUE, recursive = FALSE, mode = "0777")
        }
}
output.dir <- args[2]

genome.desc <- args[3]

if (is.na(args[3])) {
	genome.desc <- "Not defined"
}

suppressMessages(library('cummeRbund'))
suppressMessages(library('ggplot2'))

cuff <- readCufflinks(  dir = paste(input.dir,"/cuffdiff/",sep=""),
                        dbFile = "cuffData.db",
                        gtfFile = paste(input.dir,"/cuffmerge/merged.gtf",sep=""),
                        driver = "SQLite",
                        genome = genome.desc,
                        rebuild = FALSE,verbose=FALSE)


gene.features<-annotation(genes(cuff))

write.table(file=paste(output.dir,"/replicates.txt",sep=""),replicates(cuff),quote=FALSE,col.names=TRUE, row.names=FALSE, sep="\t")


## Avaliação do Ajuste do Modelo
#
## Dispersão em relação à contagem
## Lembrando: Na distribuição de Poisson a média = variância, portanto a diagonal neste gráfico indicaria uma distribuição
## teórica baseada na distribuição de Poisson. Normalmente a variabilidade biológica é maior, a dispersão é maior do que essa
## estimativa teórica.
## http://slideplayer.com/slide/6291297/
## Gráfico de dispersão comparando a média de contagens com relação à dispersão estimada.

if (!file.exists( paste(output.dir,"/DispersionPlot_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/DispersionPlot_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	disp.genes<-dispersionPlot(genes(cuff))
	plot(disp.genes)
	graphics.off()
	rm(disp.genes)
}

## Dispersão mostrando o ajuste da dispersão com relação ao modelo
if (!file.exists( paste(output.dir,"/FullDispersionPlot.png",sep="") )) {
	png(filename = paste(output.dir,"/FullDispersionPlot.png",sep=""), bg="white", res=300, width=3000, height=3000)
	disp<-dispersionPlot(cuff)
	plot(disp)
	graphics.off()
	rm(disp)
}

# NÃO FUNCIONAL NO MOMENTO "Error in seq.default(min, max, by = by) : 'from' must be a finite number"
#if (!file.exists( paste(output.dir,"/DispersionPlot_isoforms.png",sep="") )) {
#	png(filename = paste(output.dir,"/DispersionPlot_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
#	disp.isoforms<-dispersionPlot(isoforms(cuff))
#	plot(disp.isoforms)
#	graphics.off()
#	rm(disp.isoforms)	
#}

#
## Avaliação do Coeficiente de Variação - qualidade das réplicas
#
## O Quadrado do Coeficiente da Variação é uma medida normalizada de variabilidade dentro das réplicas biológicas,
## o que pode ser útil para avaliar a qualidade dos dados de RNA-seq. Diferenças em cv2 podem resultar na identificação 
## de poucos genes diferencialmente expressos
## devido ao alto grau de variabilidade entre as estivativas de fpkm das réplicas
##
## Coeficiente de Variação = desvio-padrão / média
## https://pt.wikipedia.org/wiki/Coeficiente_de_varia%C3%A7%C3%A3o
## Caso não haja réplicas ele irá calcular o coeficiente de variação para todo o conjunto de dados

if (!file.exists( paste(output.dir,"/SCVPlot_genes.png",sep="") )) {
	genes.scv<-fpkmSCVPlot(genes(cuff))
	png(filename = paste(output.dir,"/SCVPlot_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(genes.scv)
	graphics.off()
	rm(genes.scv)
}
## Idem para isoformas
if (!file.exists( paste(output.dir,"/SCVPlot_isoforms.png",sep="") )) {
	png(filename = paste(output.dir,"/SCVPlot_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
	plot(isoforms.scv)
	graphics.off()
	rm(isoforms.scv)
}

#
## Avaliação da Distribuição - identificação de outliers
#
## Distribuição de FPKM (através) por todas as amostras
## pseudocount=1 para evitar erros na transformação logarítmica de valores zerados (FPKM = 0)
## Estimativa não paramétrica da função de densidade de probabilidade (pdf) - estimativa de densidade Kernel (EDK)
## https://pt.wikipedia.org/wiki/Estimativa_de_densidade_kernel
## 
## O estimador Kernel é um interpolador, que possibilita a estimação da
## densidade em toda a área, mesmo nas regiões onde o processo
## não tenha gerado nenhuma ocorrência real.

### GENES
if (!file.exists( paste(output.dir,"/DensityPlot_genes.png",sep="") )) {
	dens.genes<-csDensity(genes(cuff),pseudocount=1, replicates=FALSE)
	png(filename = paste(output.dir,"/DensityPlot_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(dens.genes)
	graphics.off()
	rm(dens.genes)
}
## Idem para réplicas
if (!file.exists( paste(output.dir,"/DensityPlot_reps_genes.png",sep="") )) {
	dens.genes.rep<-csDensity(genes(cuff),pseudocount=1, replicates=TRUE)
	png(filename = paste(output.dir,"/DensityPlot_reps_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(dens.genes.rep)
	graphics.off()
	rm(dens.genes.rep)
}

### ISOFORMAS
if (!file.exists( paste(output.dir,"/DensityPlot_isoforms.png",sep="") )) {
	dens.isoforms<-csDensity(isoforms(cuff),pseudocount=1, replicates=FALSE)
	png(filename = paste(output.dir,"/DensityPlot_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(dens.isoforms)
	graphics.off()
	rm(dens.isoforms)
}
## Idem para réplicas
if (!file.exists( paste(output.dir,"/DensityPlot_reps_isoforms.png",sep="") )) {
	dens.isoforms.rep<-csDensity(isoforms(cuff),pseudocount=1, replicates=TRUE)
	png(filename = paste(output.dir,"/DensityPlot_reps_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(dens.isoforms.rep)
	graphics.off()
	rm(dens.isoforms.rep)
}

## BoxPlots

### GENES
if (!file.exists( paste(output.dir,"/BoxPlot_genes.png",sep="") )) {
	bxplot.genes<-csBoxplot(genes(cuff),pseudocount=1,replicates=FALSE)
	png(filename = paste(output.dir,"/BoxPlot_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(bxplot.genes)
	graphics.off()
	rm(bxplot.genes)
}
## Idem para réplicas
if (!file.exists( paste(output.dir,"/BoxPlot_rep_genes.png",sep="")  )) {
	bxplot.genes.rep<-csBoxplot(genes(cuff),pseudocount=1,replicates=TRUE)
	png(filename = paste(output.dir,"/BoxPlot_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(bxplot.genes.rep)
	graphics.off()
	rm(bxplot.genes.rep)
}

### ISOFORMAS
if (!file.exists( paste(output.dir,"/BoxPlot_isoforms.png",sep="") )) {
	bxplot.isoforms<-csBoxplot(isoforms(cuff),pseudocount=1,replicates=FALSE)
	png(filename = paste(output.dir,"/BoxPlot_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(bxplot.isoforms)
	graphics.off()
	rm(bxplot.isoforms)
}
## Idem para réplicas
if (!file.exists( paste(output.dir,"/BoxPlot_rep_isoforms.png",sep="")  )) {
	bxplot.isoforms.rep<-csBoxplot(isoforms(cuff),pseudocount=1,replicates=TRUE)
	png(filename = paste(output.dir,"/BoxPlot_rep_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(bxplot.isoforms.rep)
	graphics.off()
	rm(bxplot.isoforms.rep)
}

#
## Scatterplots (gráficos de dispersão) podem ser úteis para identificar alterações globais e tendências na expressão gênica 
## em comparações em pares

### GENES
if (!file.exists( paste(output.dir,"/ScatterMatrix_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/ScatterMatrix_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	scplot.genes<-csScatterMatrix(genes(cuff))
	plot(scplot.genes)
	graphics.off()
	rm(scplot.genes)
}

### ISOFORMAS
if (!file.exists( paste(output.dir,"/ScatterMatrix_isoforms.png",sep="") )) {
	png(filename = paste(output.dir,"/ScatterMatrix_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	scplot.isoforms<-csScatterMatrix(isoforms(cuff))
	plot(scplot.isoforms)
	graphics.off()
	rm(scplot.isoforms)
}

#
## Scatterplot - comparação para todos os pares - gráfico individual para cada comparação em par
## MAplot - M (Fold-Change) x A (Intensidade) - viés na relação fold-change e intensidade
if (! file.exists( paste(output.dir,'/ScatterPlots',sep="") ) ) {
	dir.create(paste(output.dir,'/ScatterPlots',sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
## Todas as condições
conditions.all <- unique(replicates(cuff)$sample_name)

## Todas as combinações de condições (dois a dois)
conditions.combn <- t(combn(conditions.all,2))

## Scatterplots e MAplot
for (i in 1:dim(conditions.combn)[1]) {
	# Scatterplot
	
	if (!file.exists( paste(output.dir,"/ScatterPlots","/ScatterPlot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_genes.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/ScatterPlot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		scplot.genes<-csScatter(genes(cuff),conditions.combn[i,1],conditions.combn[i,2],smooth=TRUE)
		plot(scplot.genes)
		graphics.off()
		rm(scplot.genes)
	}
	
	if (!file.exists( paste(output.dir,"/ScatterPlots","/ScatterPlot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_isoforms.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/ScatterPlot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		scplot.isoforms<-csScatter(isoforms(cuff),conditions.combn[i,1],conditions.combn[i,2],smooth=TRUE)
		plot(scplot.isoforms)
		graphics.off()
		rm(scplot.isoforms)
	}

	# MAplot
	if (!file.exists( paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_fpkm_genes.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_fpkm_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		maplot.genes<-MAplot(genes(cuff),conditions.combn[i,1],conditions.combn[i,2])
		plot(maplot.genes)
		graphics.off()
		rm(maplot.genes)
	}

	if (!file.exists( paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_fpkm_isoforms.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_fpkm_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		maplot.isoforms<-MAplot(isoforms(cuff),conditions.combn[i,1],conditions.combn[i,2])
		plot(maplot.isoforms)
		graphics.off()
		rm(maplot.isoforms)
	}


	if (!file.exists( paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_ncounts_genes.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_ncounts_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		maplotcounts.genes<-MAplot(genes(cuff),conditions.combn[i,1],conditions.combn[i,2],useCount=TRUE)
		plot(maplotcounts.genes)
		graphics.off()
		rm(maplotcounts.genes)
	}
	
	if (!file.exists( paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_ncounts_isoforms.png",sep="") )) {
		png(filename = paste(output.dir,"/ScatterPlots","/MAplot_",conditions.combn[i,1],"x",conditions.combn[i,2],"_ncounts_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		maplotcounts.isoforms<-MAplot(isoforms(cuff),conditions.combn[i,1],conditions.combn[i,2],useCount=TRUE)
		plot(maplotcounts.isoforms)
		graphics.off()
		rm(maplotcounts.isoforms)
	}

}

# Avaliação dos relacionamentos entre condições e réplicas

# Dendrograma
# HeatMaps (Leitura recomendada: https://liorpachter.wordpress.com/tag/jensen-shannon-divergence/)
# 
# w <- c(10,0,10)
# x <- c(10,0,0)
# y <- c(10,5,5)
# z <- c(10,10,0)
#
# Euclidean Distance: 
#    EUD <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
#
#    EUD(x,y) # 7.071068   
#    EUD(x,z) # 10
# 
# * influência da magnitude da expressão dos genes
#
# Distância de Jensen-Shannon (medida de dissimilaridade entre distribuições de probabilidade)
#    (comparação entre distribuições): 
#    (http://enterotype.embl.de/enterotypes.html)
#
#    KLD <- function(x,y) sum(x * log(x/y))
#    JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
#
#    JSD(x+1e-10,y+1e-10) # 1.861649
#    JSD(x+1e-10,z+1e-10) # 1.861649
#
# * influência da magnitude da expressão dos genes (dominado pelos de maior expressão => maior variância)
#
# Mahalanobis Distance:
#    (distância definida pelas componentes principais)
#
#    (adicionado ruído para solucionar problema de desvio padrão igual a zero)    
#    dat <- t(apply(rbind(w, x, y, z), 1, function(x) { return (x+runif(1,min=0,max=1e-10)); } ))
#
#    (normalização = (valor - média)/desvio
#    (covariância entre variáveis normalizadas são correlações - depois da normalização cada variável tem variância 1)
#    
#    (tol - tolerância - alguns números usados no cálculo da matriz inversa podem ser muito pequenos e tratados como zero dando a entender 
#    que a matriz é singular, este parâmetro define esse limiar do que é muito pequeno. - Matriz singular não admite inversa - determinante é nulo)
#    
#    mahalanobis(x, y, cov=cor(dat),tol=1e-10) # 99.53442
#    mahalanobis(x, z, cov=cor(dat),tol=1e-10) # 141.4435
#    mahalanobis(x, w, cov=cor(dat),tol=1e-10) # 216.9467
#   
# * influência reduzida da variância
#
if (!file.exists( paste(output.dir,"/Dendrogram_JS_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/Dendrogram_JS_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	par(mar=c(10,0,0,0))
	dend.genes<-csDendro(genes(cuff), pseudocount=1, logMode=TRUE, replicates=FALSE)
	plot(dend.genes)
	graphics.off()
	rm(dend.genes)
}

if (!file.exists( paste(output.dir,"/Dendrogram_JS_rep_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/Dendrogram_JS_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	par(mar=c(12,0,0,0))
	dend.genes.rep<-csDendro(genes(cuff), pseudocount=1, logMode=TRUE, replicates=TRUE)
	plot(dend.genes.rep)
	graphics.off()
	rm(dend.genes.rep)
}	

# HeatMaps
#
if (!file.exists( paste(output.dir,"/DistHeatmap_JS_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/DistHeatmap_JS_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	myDistHeat.genes<-csDistHeat(genes(cuff))
	plot(myDistHeat.genes)
	graphics.off()
	rm(myDistHeat.genes)
}

if (!file.exists( paste(output.dir,"/DistHeatmap_JS_rep_genes.png",sep="") )) {
	png(filename = paste(output.dir,"/DistHeatmap_JS_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	myDistHeat.genes.rep<-csDistHeat(genes(cuff), replicates=TRUE)
	plot(myDistHeat.genes.rep)
	graphics.off()
	rm(myDistHeat.genes.rep)
}

if (!file.exists( paste(output.dir,"/DistHeatmap_JS_isoforms.png",sep="") )) {
	png(filename = paste(output.dir,"/DistHeatmap_JS_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	myDistHeat.isoforms<-csDistHeat(isoforms(cuff))
	plot(myDistHeat.isoforms)
	graphics.off()
	rm(myDistHeat.isoforms)
}

if (!file.exists( paste(output.dir,"/DistHeatmap_JS_rep_isoforms.png",sep="") )) {
	png(filename = paste(output.dir,"/DistHeatmap_JS_rep_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	myDistHeat.isoforms.rep<-csDistHeat(isoforms(cuff), replicates=TRUE)
	plot(myDistHeat.isoforms.rep)
	graphics.off()
	rm(myDistHeat.isoforms.rep)
}


# Redução de Dimensionalidade
# Agrupamento e exploração das relações entre condições

# PCA - Principal Component Analysis
if (!file.exists( paste(output.dir,"/PCA_PC1_PC2_genes.png",sep="") )) {
	genes.PCA<-PCAplot(genes(cuff),pseudocount=1,x="PC1",y="PC2",replicates=FALSE)
	png(filename = paste(output.dir,"/PCA_PC1_PC2_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(genes.PCA)
	graphics.off()
	rm(genes.PCA)
}

if (!file.exists( paste(output.dir,"/PCA_PC1_PC2_isoforms.png",sep="") )) {
	isoforms.PCA<-PCAplot(isoforms(cuff),pseudocount=1,x="PC1",y="PC2",replicates=FALSE)
	png(filename = paste(output.dir,"/PCA_PC1_PC2_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(isoforms.PCA)
	graphics.off()
	rm(isoforms.PCA)
}

if (!file.exists( paste(output.dir,"/PCA_PC1_PC2_rep_genes.png",sep="") )) {
	genes.PCA.rep<-PCAplot(genes(cuff),pseudocount=1,x="PC1",y="PC2",replicates=TRUE)
	png(filename = paste(output.dir,"/PCA_PC1_PC2_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(genes.PCA.rep)
	graphics.off()
	rm(genes.PCA.rep)
}

if (!file.exists( paste(output.dir,"/PCA_PC1_PC2_rep_isoforms.png",sep="") )) {
	isoforms.PCA.rep<-PCAplot(isoforms(cuff),pseudocount=1,x="PC1",y="PC2",replicates=TRUE)
	png(filename = paste(output.dir,"/PCA_PC1_PC2_rep_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(isoforms.PCA.rep)
	graphics.off()
	rm(isoforms.PCA.rep)
}

# MDS - Multi-Dimensional Scaling

if (length(samples(genes(cuff))) > 2) {
	if (!file.exists( paste(output.dir,"/MDS_M1_M2_genes.png",sep="") )) {
		genes.MDS<-MDSplot(genes(cuff),pseudocount=1, replicates=FALSE)
		png(filename = paste(output.dir,"/MDS_M1_M2_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(genes.MDS)
		graphics.off()
		rm(genes.MDS)
	}
	
	if (!file.exists( paste(output.dir,"/PCA_PC2_PC3_genes.png",sep="") )) {
		genes.PCA.rep<-PCAplot(genes(cuff),pseudocount=1,x="PC2",y="PC3",replicates=FALSE)
		png(filename = paste(output.dir,"/PCA_PC2_PC3_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(genes.PCA.rep)
		graphics.off()
		rm(genes.PCA.rep)
	}
	
	if (!file.exists( paste(output.dir,"/PCA_PC1_PC3_genes.png",sep="") )) {
		genes.PCA.rep<-PCAplot(genes(cuff),pseudocount=1,x="PC1",y="PC3",replicates=FALSE)
		png(filename = paste(output.dir,"/PCA_PC1_PC3_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(genes.PCA.rep)
		graphics.off()
		rm(genes.PCA.rep)
	}
}

if (length(samples(isoforms(cuff))) > 2) {
	if (!file.exists( paste(output.dir,"/MDS_M1_M2_isoforms.png",sep="") )) {
		isoforms.MDS<-MDSplot(isoforms(cuff),pseudocount=1, replicates=FALSE)
		png(filename = paste(output.dir,"/MDS_M1_M2_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(isoforms.MDS)
		graphics.off()
		rm(isoforms.MDS)
	}

	if (!file.exists( paste(output.dir,"/PCA_PC2_PC3_isoforms.png",sep="") )) {
		isoforms.PCA.rep<-PCAplot(isoforms(cuff),pseudocount=1,x="PC2",y="PC3",replicates=FALSE)
		png(filename = paste(output.dir,"/PCA_PC2_PC3_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(isoforms.PCA.rep)
		graphics.off()
		rm(isoforms.PCA.rep)
	}
	
	if (!file.exists( paste(output.dir,"/PCA_PC1_PC3_isoforms.png",sep="") )) {
		isoforms.PCA.rep<-PCAplot(isoforms(cuff),pseudocount=1,x="PC1",y="PC3",replicates=FALSE)
		png(filename = paste(output.dir,"/PCA_PC1_PC3_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(isoforms.PCA.rep)
		graphics.off()
		rm(isoforms.PCA.rep)
	}
}

if (length(replicates(genes(cuff))) > 2) {
	if (!file.exists( paste(output.dir,"/MDS_M1_M2_rep_genes.png",sep="") )) {
		genes.MDS.rep<-MDSplot(genes(cuff),pseudocount=1, replicates=TRUE)
		png(filename = paste(output.dir,"/MDS_M1_M2_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(genes.MDS.rep)
		graphics.off()
		rm(genes.MDS.rep)
	}

}

if (length(replicates(isoforms(cuff))) > 2) {
	if (!file.exists( paste(output.dir,"/MDS_M1_M2_rep_isoforms.png",sep="") )) {
		isoforms.MDS.rep<-MDSplot(isoforms(cuff),pseudocount=1, replicates=TRUE)
		png(filename = paste(output.dir,"/MDS_M1_M2_rep_isoforms.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(isoforms.MDS.rep)
		graphics.off()
		rm(isoforms.MDS.rep)
	}
}


my.genes.matrix <- repFpkmMatrix(genes(cuff))
my.isoforms.matrix <- repFpkmMatrix(isoforms(cuff))


suppressMessages(library('pvclust'))


if (!file.exists( paste(output.dir,"/pvclust_genes_cor_avg_100.png",sep="") )) {
	result.genes.cor <- pvclust(my.genes.matrix, method.dist="correlation", method.hclust="average", nboot=100)
	png(filename = paste(output.dir,"/pvclust_genes_cor_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(result.genes.cor)
	graphics.off()
	rm(result.genes.cor)
}

if (!file.exists( paste(output.dir,"/pvclust_genes_euc_avg_100.png", sep="") )) {
	result.genes.euc <- pvclust(my.genes.matrix, method.dist="euclidean", method.hclust="average", nboot=100)
	png(filename = paste(output.dir,"/pvclust_genes_euc_avg_100.png", sep=""), bg="white", res=300, width=3000, height=3000)
	plot(result.genes.euc)
	graphics.off()
	rm(result.genes.euc)
}

if (!file.exists( paste(output.dir,"/pvclust_isoforms_cor_avg_100.png",sep="") )) {
	result.isoforms.cor <- pvclust(my.isoforms.matrix, method.dist="correlation", method.hclust="average", nboot=100)
	png(filename = paste(output.dir,"/pvclust_isoforms_cor_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
	plot(result.isoforms.cor)
	graphics.off()
	rm(result.isoforms.cor)
}

if (!file.exists( paste(output.dir,"/pvclust_isoforms_euc_avg_100.png", sep="") )) {
	result.isoforms.euc <- pvclust(my.isoforms.matrix, method.dist="euclidean", method.hclust="average", nboot=100)
	png(filename = paste(output.dir,"/pvclust_isoforms_euc_avg_100.png", sep=""), bg="white", res=300, width=3000, height=3000)
	plot(result.isoforms.euc)
	graphics.off()
	rm(result.isoforms.euc)
}

#suppressMessages(library('HDMD'))
#
#if (!file.exists( paste(output.dir,"/pvclust_genes_mah_avg_100.png", sep="") )) {
#	result.genes.mah <- pvclust(my.genes.matrix, method.dist=function(x){ as.dist( pairwise.mahalanobis(x,grouping=c(1:dim(x)[1]))$distance) }, method.hclust="average", nboot=100)
#	png(filename = paste(output.dir,"/pvclust_mah_avg_100.png", sep=""), bg="white", res=300, width=3000, height=3000)
#	plot(result.genes.mah)
#	graphics.off()
#	rm(result.genes.mah)
#}
#
#if (!file.exists( paste(output.dir,"/pvclust_isoforms_mah_avg_100.png", sep="") )) {
#	result.isoforms.mah <- pvclust(my.isoforms.matrix, method.dist=function(x){ as.dist( pairwise.mahalanobis(x,grouping=c(1:dim(x)[1]))$distance) }, method.hclust="average", nboot=100)
#	png(filename = paste(output.dir,"/pvclust_mah_avg_100.png", sep=""), bg="white", res=300, width=3000, height=3000)
#	plot(result.isoforms.mah)
#	graphics.off()
#	rm(result.isoforms.mah)
#}


# http://bioinf.wehi.edu.au/RNAseqCaseStudy/
# voom normalization

# https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/

suppressMessages(library('limma'))
suppressMessages(library('edgeR'))

if ((! file.exists( paste(output.dir,"/pvclust_genes_voom_cor_avg_100.png",sep="") ) )||
    (! file.exists( paste(output.dir,"/pvclust_genes_voom_euc_avg_100.png",sep="") ) ) ) {
	
	samp.genes.counts <- reshape(repFpkm(genes(cuff))[,c('gene_id','rep_name','raw_frags')], direction="wide", 
						  idvar=c('gene_id'), timevar=c("rep_name"))
	samp.genes.counts.rownames <- samp.genes.counts$gene_id
	samp.genes.counts <- samp.genes.counts[,-c(which(colnames(samp.genes.counts)=='gene_id'))]
	rownames(samp.genes.counts) <- samp.genes.counts.rownames
	colnames(samp.genes.counts) <- gsub('raw_frags.','',colnames(samp.genes.counts))
	
	dge.genes <- DGEList(counts=samp.genes.counts)
	dge.genes <- calcNormFactors(dge.genes)
	
	design.genes<-model.matrix(~as.factor(replicates(cuff)$sample_name))
	isexpr.genes <- rowSums(cpm(dge.genes) > 10) >= 2
	dge.genes <- dge.genes[isexpr.genes,]
	
	v.genes <- voom(dge.genes,design.genes,plot=FALSE)
	
	if (! file.exists( paste(output.dir,"/pvclust_genes_voom_cor_avg_100.png",sep="") ) ) {
		result.genes.voom.cor <- pvclust(v.genes$E, method.dist="correlation", method.hclust="average", nboot=100)
		png(filename = paste(output.dir,"/pvclust_genes_voom_cor_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(result.genes.voom.cor)
		graphics.off()
		rm(result.genes.voom.cor)
	}
	
    	if (! file.exists( paste(output.dir,"/pvclust_genes_voom_euc_avg_100.png",sep="") ) ) {
		result.genes.voom.euc <- pvclust(v.genes$E, method.dist="euclidean", method.hclust="average", nboot=100)
		png(filename = paste(output.dir,"/pvclust_genes_voom_euc_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(result.genes.voom.euc)
		graphics.off()
		rm(result.genes.voom.euc)

	}
}

if ((! file.exists( paste(output.dir,"/pvclust_isoforms_voom_cor_avg_100.png",sep="") ) )||
    (! file.exists( paste(output.dir,"/pvclust_isoforms_voom_euc_avg_100.png",sep="") ) ) ) {
	
	samp.isoforms.counts <- reshape(repFpkm(isoforms(cuff))[,c('isoform_id','rep_name','raw_frags')], direction="wide", 
						  idvar=c('isoform_id'), timevar=c("rep_name"))
	samp.isoforms.counts.rownames <- samp.isoforms.counts$gene_id
	samp.isoforms.counts <- samp.isoforms.counts[,-c(which(colnames(samp.isoforms.counts)=='isoform_id'))]
	rownames(samp.isoforms.counts) <- samp.isoforms.counts.rownames
	colnames(samp.isoforms.counts) <- gsub('raw_frags.','',colnames(samp.isoforms.counts))
	
	dge.isoforms <- DGEList(counts=samp.isoforms.counts)
	dge.isoforms <- calcNormFactors(dge.isoforms)
	
	design.isoforms<-model.matrix(~as.factor(replicates(cuff)$sample_name))
	isexpr.isoforms <- rowSums(cpm(dge.isoforms) > 10) >= 2
	dge.isoforms <- dge.isoforms[isexpr.isoforms,]
	
	v.isoforms <- voom(dge.isoforms,design.isoforms,plot=FALSE)
	
	if (! file.exists( paste(output.dir,"/pvclust_isoforms_voom_cor_avg_100.png",sep="") ) ) {
		result.isoforms.voom.cor <- pvclust(v.isoforms$E, method.dist="correlation", method.hclust="average", nboot=100)
		png(filename = paste(output.dir,"/pvclust_isoforms_voom_cor_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(result.isoforms.voom.cor)
		graphics.off()
		rm(result.isoforms.voom.cor)
	}
	
    	if (! file.exists( paste(output.dir,"/pvclust_isoforms_voom_euc_avg_100.png",sep="") ) ) {
		result.isoforms.voom.euc <- pvclust(v.isoforms$E, method.dist="euclidean", method.hclust="average", nboot=100)
		png(filename = paste(output.dir,"/pvclust_isoforms_voom_euc_avg_100.png",sep=""), bg="white", res=300, width=3000, height=3000)
		plot(result.isoforms.voom.euc)
		graphics.off()
		rm(result.isoforms.voom.euc)

	}
}

quit(save="no")

#suppressMessages(library('MASS'))
#
#test.mds <- cmdscale(
#  dist(
#    t(as.matrix(v$E))
#  )
#)
#
#Dim1 <- test.mds [,1]
#Dim2 <- test.mds [,2]
#
###
### Plot the solution. 
###
#
#png(filename = paste(output.dir,"/MDS_M1_M2_voom_rep_genes.png",sep=""), bg="white", res=300, width=3000, height=3000)
#plot(Dim1, Dim2, xlab="", ylab="", main="MDS")
#segments(-1500, -0, 1500, 0, lty="dotted")
#segments(0, -1500, 0, 1500, lty="dotted")
#text(Dim1, Dim2, adj = c(1, 2), rownames(test.mds), srt=45, cex=0.6, col="red")
#graphics.off()

#temp <- my.gene.matrix 
#samps <- colnames(my.gene.matrix)
#
#temp$SD <- apply(temp[,samps],1,sd,na.rm=T) #sd across all samples
#temp$VAR <- apply(temp[,samps],1,var,na.rm=T) #var across all samples
#temp$MAD <- apply(temp[,samps],1,mad,na.rm=T) #mad across all samples
#
#var.cutoff <- quantile(temp$VAR,0.75)
#
### VAR Plot
#m <- ggplot(temp, aes(x=SD))  +
#  ggtitle("Density Plot\nVariance across all samples")
#m + geom_histogram(aes(y=..density..), binwidth=0.1, fill="gray", colour="black") +
#  geom_density(fill=NA, colour="green", size=0.5, alpha=0.2) +
#  theme(title=element_text(size=15, vjust= 0.25, colour="black", face = "bold")) +
#  scale_y_continuous("Density") +
#  geom_vline(xintercept= var.cutoff, colour="red", linetype=2, size=0.5) +
#  scale_x_log10("Variance")
#
#ggsave(file = paste(output.dir,"/VARplot_genes.png",sep=''), width=10, height=10)
#graphics.off()
#
#temp.var = subset(temp,VAR>=var.cutoff)
#
#
#myGenes<-getGenes(cuff,myGeneIds)
#th<-csHeatmap(myGenes,cluster='both',labRow=F) 


