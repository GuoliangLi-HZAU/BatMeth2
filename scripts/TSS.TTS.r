#inputPrefix[${prefix}.TSS.1.txt.sorted].${context} outputPrefix color_cg color_chg color_chh cg-ceil chg-ceil chh-ceil
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)

#library(RColorBrewer)
library(pheatmap)
########################## CG
l<-read.table(paste(args[1],".cg",sep=""),sep="\t")
l<-l[,2:162]
y<-data.matrix(l)
	ceil <- as.double(args[6])
pdf(paste(args[2],".cg.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[3]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],".cg.png",sep=""),res=128)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[3]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CG down!")
########################## CHG
l<-read.table(paste(args[1],".chg",sep=""),sep="\t")
l<-l[,2:162]
y<-data.matrix(l)
#####CHG max for heatmap
       ceil<- as.double(args[7])
#################
pdf(paste(args[2],".chg.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[4]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],".chg.png",sep=""),res=128)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[4]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CHG down!")
######################### CHH
l<-read.table(paste(args[1],".chh",sep=""),sep="\t")
l<-l[,2:162]
y<-data.matrix(l)
######CHH max for heatmap
       ceil<- as.double(args[8])
##############
pdf(paste(args[2],".chh.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[5]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],".chh.png",sep=""),res=128)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[5]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CHH down!")
