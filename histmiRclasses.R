# TODO: Add comment
# 
# Author: dgpinheiro
###############################################################################

library('ggplot2')
library('reshape')
library('gridExtra')
require('grid')


act.passed.df <- read.delim(file="/tmp/PASSED-LEN-DIST-ACT.txt",header=TRUE)
rownames(act.passed.df) <- act.passed.df$Length
colnames(act.passed.df)
head(act.passed.df[,setdiff(colnames(act.passed.df), 'Length')])
act.passed.df.m <- melt(act.passed.df, id.var=c("Length"),variable_name=c("class") )
act.passed.df.m$Length <- as.factor(act.passed.df.m$Length)
head(act.passed.df.m)


for.passed.df <- read.delim(file="/tmp/PASSED-LEN-DIST-FOR.txt",header=TRUE)
rownames(for.passed.df) <- for.passed.df$Length
colnames(for.passed.df)
head(for.passed.df[,setdiff(colnames(for.passed.df), 'Length')])
for.passed.df.m <- melt(for.passed.df, id.var=c("Length"),variable_name=c("class") )
for.passed.df.m$Length <- as.factor(for.passed.df.m$Length)
head(for.passed.df.m)



nur.passed.df <- read.delim(file="/tmp/PASSED-LEN-DIST-NUR.txt",header=TRUE)
rownames(nur.passed.df) <- nur.passed.df$Length
colnames(nur.passed.df)
head(nur.passed.df[,setdiff(colnames(nur.passed.df), 'Length')])
nur.passed.df.m <- melt(nur.passed.df, id.var=c("Length"),variable_name=c("class") )
nur.passed.df.m$Length <- as.factor(nur.passed.df.m$Length)
head(nur.passed.df.m)

act.class.selected <- unique(as.character(act.passed.df.m[order(act.passed.df.m$value,decreasing=TRUE),'class']))[1:7]
for.class.selected <- unique(as.character(for.passed.df.m[order(for.passed.df.m$value,decreasing=TRUE),'class']))[1:7]
nur.class.selected <- unique(as.character(nur.passed.df.m[order(nur.passed.df.m$value,decreasing=TRUE),'class']))[1:7]

class.selected <- unique(c(act.class.selected, for.class.selected, nur.class.selected))

class.other<- setdiff(unique(c(
						unique(as.character(act.passed.df.m[order(act.passed.df.m$value,decreasing=TRUE),'class'])),
						unique(as.character(for.passed.df.m[order(for.passed.df.m$value,decreasing=TRUE),'class'])),
						unique(as.character(nur.passed.df.m[order(nur.passed.df.m$value,decreasing=TRUE),'class']))
				)
		), class.selected)


#cbPalette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
cbPalette <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#636363')
length(cbPalette)

names(cbPalette) <-c(class.selected, 'other')


### ACT

act.new.passed.df.m <- as.data.frame(matrix(ncol=3, nrow=length(unique(act.passed.df$Length))*length(c(class.selected,'other')) ))
colnames(act.new.passed.df.m) <- c('Length', 'class', 'value') 

head(act.new.passed.df.m)
x<-1
for (l in unique(act.passed.df$Length)) {
	#print(l)
	for (c in c(class.selected)) {
		#print(paste(l,c))
		act.new.passed.df.m[x, c('Length','class')] <- c(l, c)
		act.new.passed.df.m[x, 'value'] <- as.numeric( subset(act.passed.df.m, Length==l & class==c)[,'value'] )
		x=x+1
	}
	act.new.passed.df.m[x,c('Length','class','value')] <- c(l, 'other', as.numeric(sum( as.numeric(subset(act.passed.df.m, Length==l & class %in% class.other)[,'value']) )))
	x=x+1
}
head(act.new.passed.df.m,10)

act.new.passed.df.m$Length <- as.factor(act.new.passed.df.m$Length)

act.new.passed.df.m[,'value'] <- as.numeric(act.new.passed.df.m[,'value']) 



### FOR

for.new.passed.df.m <- as.data.frame(matrix(ncol=3, nrow=length(unique(for.passed.df$Length))*length(c(class.selected,'other')) ))
colnames(for.new.passed.df.m) <- c('Length', 'class', 'value') 

head(for.new.passed.df.m)
x<-1
for (l in unique(for.passed.df$Length)) {
	#print(l)
	for (c in c(class.selected)) {
		#print(paste(l,c))
		for.new.passed.df.m[x, c('Length','class')] <- c(l, c)
		for.new.passed.df.m[x, 'value'] <- as.numeric( subset(for.passed.df.m, Length==l & class==c)[,'value'] )
		x=x+1
	}
	for.new.passed.df.m[x,c('Length','class','value')] <- c(l, 'other', as.numeric(sum( as.numeric(subset(for.passed.df.m, Length==l & class %in% class.other)[,'value']) )))
	x=x+1
}
head(for.new.passed.df.m,10)

for.new.passed.df.m$Length <- as.factor(for.new.passed.df.m$Length)

for.new.passed.df.m[,'value'] <- as.numeric(for.new.passed.df.m[,'value']) 




### NUR

nur.new.passed.df.m <- as.data.frame(matrix(ncol=3, nrow=length(unique(nur.passed.df$Length))*length(c(class.selected,'other')) ))
colnames(nur.new.passed.df.m) <- c('Length', 'class', 'value') 

head(nur.new.passed.df.m)
x<-1
for (l in unique(nur.passed.df$Length)) {
	#print(l)
	for (c in c(class.selected)) {
		#print(paste(l,c))
		nur.new.passed.df.m[x, c('Length','class')] <- c(l, c)
		nur.new.passed.df.m[x, 'value'] <- as.numeric( subset(nur.passed.df.m, Length==l & class==c)[,'value'] )
		x=x+1
	}
	nur.new.passed.df.m[x,c('Length','class','value')] <- c(l, 'other', as.numeric(sum( as.numeric(subset(nur.passed.df.m, Length==l & class %in% class.other)[,'value']) )))
	x=x+1
}
head(nur.new.passed.df.m,10)

nur.new.passed.df.m$Length <- as.factor(nur.new.passed.df.m$Length)

nur.new.passed.df.m[,'value'] <- as.numeric(nur.new.passed.df.m[,'value']) 


ggp.mir <- list()

### ACT


act.new.passed.df.m$class <-
		gsub("mature", "miRNA.mature",
				gsub("hairpin","miRNA.hairpin",
						gsub("frna.ncrna", "fRNAdb.ncRNA",
								gsub("dbest","Amel.EST",
										gsub("rrna","rRNA",
												gsub("trna","Amel.tRNA",
														gsub("mirdeep2", "M", gsub("genome", "G", gsub("other", "...", gsub("no.mirdeep2"," ",gsub("no.genome", " ", gsub("_"," - ",
																										act.new.passed.df.m$class
																								))))))))))))

for.new.passed.df.m$class <-
		gsub("mature", "miRNA.mature",
				gsub("hairpin","miRNA.hairpin",
						gsub("frna.ncrna", "fRNAdb.ncRNA",
								gsub("dbest","Amel.EST",
										gsub("rrna","rRNA",
												gsub("trna","Amel.tRNA",
														gsub("mirdeep2", "M", gsub("genome", "G", gsub("other", "...", gsub("no.mirdeep2"," ",gsub("no.genome", " ", gsub("_"," - ",
																										for.new.passed.df.m$class
																								))))))))))))

nur.new.passed.df.m$class <-
		gsub("mature", "miRNA.mature",
				gsub("hairpin","miRNA.hairpin",
						gsub("frna.ncrna", "fRNAdb.ncRNA",
								gsub("dbest","Amel.EST",
										gsub("rrna","rRNA",
												gsub("trna","Amel.tRNA",
														gsub("mirdeep2", "M", gsub("genome", "G", gsub("other", "...", gsub("no.mirdeep2"," ",gsub("no.genome", " ", gsub("_"," - ",
																										nur.new.passed.df.m$class
																								))))))))))))

names(cbPalette) <-
		gsub("mature", "miRNA.mature",
				gsub("hairpin","miRNA.hairpin",
						gsub("frna.ncrna", "fRNAdb.ncRNA",
								gsub("dbest","Amel.EST",
										gsub("rrna","rRNA",
												gsub("trna","Amel.tRNA",
														gsub("mirdeep2", "M", gsub("genome", "G", gsub("other", "...", gsub("no.mirdeep2"," ",gsub("no.genome", " ", gsub("_"," - ",
																										c(class.selected, 'other')
																								))))))))))))




vlimits=c(0,14)
vbreas=0:15


a <- ggplot(act.new.passed.df.m, aes(x = Length, y = value/1e+06,
						fill = class)) + ggtitle("ACT") + 
		theme(plot.title = element_text(lineheight=.8, face="bold"))+
		xlab("Length")+ylab("Number of reads (in millions)")+theme(legend.text = element_text(colour = 'black', angle = 0, size = 10, hjust = 1, vjust = 1, face = 'bold'),legend.key.size = unit(1, "cm") )+
		scale_x_discrete()+scale_fill_manual(values=cbPalette)+scale_y_continuous(limits = vlimits, breaks= vbreaks)

ggp.mir[['act']]  <- a + geom_bar(stat = "identity", position = "stack", width=0.9)+theme(plot.margin=unit(x=c(1,1,1,1),units="cm") )
#ggp.mir[['act']] 
#ggsave(file="/tmp/NEW/ACT-length-dist.png",width=25, height=10)


### FOR

a <- ggplot(for.new.passed.df.m, aes(x = Length, y = value/1e+06,
						fill = class)) + ggtitle("FOR") + 
		theme(plot.title = element_text(lineheight=.8, face="bold"))+
		xlab("Length")+ylab("Number of reads (in millions)")+theme(legend.text = element_text(colour = 'black', angle = 0, size = 10, hjust = 1, vjust = 1, face = 'bold'),legend.key.size = unit(1, "cm") )+
		scale_x_discrete()+scale_fill_manual(values=cbPalette)+scale_y_continuous(limits = vlimits, breaks= vbreaks)

ggp.mir[['for']]  <- a + geom_bar(stat = "identity", position = "stack", width=0.9)+theme(plot.margin=unit(x=c(1,1,1,1),units="cm") )
#ggp.mir[['for']] 
#ggsave(file="/tmp/NEW/FOR-length-dist.png",width=25, height=10)


### NUR

a <- ggplot(nur.new.passed.df.m, aes(x = Length, y = value/1e+06,
						fill = class)) + ggtitle("NUR") + 
		theme(plot.title = element_text(lineheight=.8, face="bold"))+
		xlab("Length")+ylab("Number of reads (in millions)")+theme(legend.text = element_text(colour = 'black', angle = 0, size = 10, hjust = 1, vjust = 1, face = 'bold'),legend.key.size = unit(1, "cm") )+
		scale_x_discrete()+scale_fill_manual(values=cbPalette)+scale_y_continuous(limits = vlimits, breaks= vbreaks)

ggp.mir[['nur']]  <- a + geom_bar(stat = "identity", position = "stack", width=0.9)+theme(plot.margin=unit(x=c(1,1,1,1),units="cm") )
#ggp.mir[['nur']] 
#ggsave(file="/tmp/NEW/NUR-length-dist.png",width=25, height=10)


#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)}

mylegend<-g_legend(ggp.mir[['act']]+theme(text=element_text(family="Courier")))

dev.off()

tiff(filename = "/tmp/miRNA-classes.tiff", height = 3000, width = 3000, compression="lzw", res=300, units="px", type="cairo")

grid.arrange(arrangeGrob( ggp.mir[['act']] + theme(legend.position="none"),
				ggp.mir[['for']] + theme(legend.position="none"),
				ggp.mir[['nur']] + theme(legend.position="none"),
				main = textGrob( "Histogram of read lengths from small RNA-Seq libraries", vjust = 1, gp = gpar( fontface = "bold", cex = 1.5 )),
				nrow=3,
				heights  = unit( 1 / 3, "npc" ),
				ncol=1),
		mylegend, nrow=1,ncol=2, widths  = c(0.75,0.25))

dev.off()



