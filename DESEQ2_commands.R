library("DESeq2")
cts <- as.matrix(read.csv("counttable.txt",sep="\t",header=TRUE,check.names=FALSE,row.names=1))
head(cts)
coldata <- read.csv("coldata.txt",sep="\t",header=TRUE,check.names=FALSE,row.names=1)
head(coldata)
coldata <- coldata[,c("condition","type")]
dds<-DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
dds<-DESeq(dds)

library(ggplot2)
library(ggrepel)
rld<-rlog(dds,blind=FALSE)
data <- plotPCA(rld, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
qplot(data$PC1, data$PC2, color=data$group, data=data) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_point(size=2) + geom_label_repel(aes(label = dds$condition),box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50')



counts <- counts(dds, normalized=TRUE)
write.table(counts, file = "normalized_counts.txt", sep = "\t", quote = FALSE)

comp1 <- results(dds, contrast = c("condition", "HL", "CT"))
resultsNames(comp1)
res1sig <- comp1[which(comp1$padj<0.05),]
write.table(res1sig, file = "Comp1_CTvsHL.txt", sep = "\t", quote = FALSE)

comp2 <- results(dds, contrast = c("condition", "HS", "CT"))
resultsNames(comp2)
res2sig <- comp2[which(comp2$padj<0.05),]
write.table(res2sig, file = "Comp2_CTvsHS.txt", sep = "\t", quote = FALSE)

comp3 <- results(dds, contrast = c("condition", "HL_HS", "CT"))
resultsNames(comp3)
res3sig <- comp3[which(comp3$padj<0.05),]
write.table(res3sig, file = "Comp3_CTvsHL_HS.txt", sep = "\t", quote = FALSE)







