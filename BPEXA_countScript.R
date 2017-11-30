#Import libaries:
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

#Read tsv file:
table<-read.table("total_countTable.tsv", header=T, row.names=1)
specialAttr<-tail(table, 5)
countTable<-head(table, -5)

#Remove lowly expressed genes:
myCPM<- cpm(countTable)
thresh<- myCPM > 0.5
keep<- rowSums(thresh) >= 2
counts.keep<- countTable[keep,]

#Remove globin genes
globin<-c("HBB", "B2M", "HBA1", "HBA2", "ARHGDIB")
counts.keep.filtered<-counts.keep[!rownames(counts.keep) %in% globin, ]

#Clear some memory:
remove(countTable, thresh, keep, table, myCPM)

#DGEList:
ObjectDGE<-DGEList(counts.keep.filtered)
logcounts<-cpm(ObjectDGE, log=T)

#BARPLOT OF LIBRARY SIZES:
png(file="Barplot_LibrarySizes.png")
barplot(ObjectDGE$samples$lib.size,names=colnames(ObjectDGE),las=2)
title("Barplot of library sizes")
dev.off()

#MDSPLOT:
png(file="MDSPlot_DGE.png")
plotMDS(ObjectDGE)
dev.off()

#HEATMAP top 500 most variable genes:
var_genes<-apply(logcounts, 1, var)
select_var<-names(sort(var_genes, decreasing = T))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
png(file="High_var_genes.heatmap.png", width=900, height=600, res=100)
heatmap.2(highly_variable_lcpm,
          col=rev(morecols(50)),
          trace="none", 
          main="Top 500 most variable genes across samples",
          scale="row", 
          ColSideColors=c("Red","green","Yellow"))
dev.off()
remove(var_genes, select_var, mypalette, morecols, highly_variable_lcpm)

#Normalise:
NormObjectDGE<-calcNormFactors(ObjectDGE)
remove(ObjectDGE)

#Voom transform:
png(file="VoomPlot_meanVariance.png")
voomObject<-voom(NormObjectDGE, plot=T)
dev.off()
#remove(NormObjectDGE)

#BOXPLOT UNNORMALISED vs NORMALISED:
png(file="Boxplot_VoomTransformedvsUnnormalised.png")
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(voomObject$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(voomObject$E),col="blue")
dev.off()
remove(logcounts)

#Contrast matrices:
cont.P1P2<-makeContrasts(P1vsP2=PT1-PT2, levels=voomObject$E)
cont.P1P3<-makeContrasts(P1vsP3=PT1-PT3, levels=voomObject$E)
cont.P2P3<-makeContrasts(P2vsP3=PT2-PT3, levels=voomObject$E)
cont.matrix<-cbind(cont.P1P2,cont.P1P3,cont.P2P3)
remove(cont.P1P2, cont.P1P3, cont.P2P3)


#DEG matrix:
contrast.matrix <- cbind(PT1=c(1,0),PT2=c(0,1),"PT3"=c(-1,1))
design<-cbind(V1=c(1,1,0),V2=c(0,0,1))
fit <- lmFit(voomObject$E, design=design)
cf<-contrasts.fit(fit, contrast.matrix)
cf<-eBayes(cf)
