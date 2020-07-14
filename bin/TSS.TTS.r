#inputPrefix[${prefix}.TSS.1.txt.sorted].${context} outputPrefix color_cg color_chg color_chh cg-ceil chg-ceil chh-ceil
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)

##Check if necessary libraries are installed
check_pkg <- function(pkg) {
  if(require(pkg, character.only = TRUE)){
    print(paste("Package", pkg, "is loaded correctly", sep = " "))
  } else {
    print(paste("Trying to install package", pkg, sep = " "))
    install.packages(pkg, repos="http://cran.us.r-project.org", dep = TRUE)
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package", pkg, sep = " "));
      }
    }
    if(require(pkg, character.only = TRUE)){
      print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
    } else{
      install.packages(pkg, repos="https://mirrors.shu.edu.cn", dep = TRUE)
      if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
      } else{
        stop(paste("Couldn't install package", pkg, sep = " "));
      }
    }
  }
}
check_pkg("pheatmap")

#library(RColorBrewer)
library(pheatmap)
########################## CG
l<-read.table(paste(args[1],".cg",sep=""),sep="\t")
l<-l[,2:ncol(l)]
y<-data.matrix(l)
	ceil <- as.double(args[6])
pdf(paste(args[2],"-CG.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[3]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],"-CG.png",sep=""),res=200, width=460, height=560)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[3]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CG down!")
########################## CHG
l<-read.table(paste(args[1],".chg",sep=""),sep="\t")
l<-l[,2:ncol(l)]
y<-data.matrix(l)
#####CHG max for heatmap
       ceil<- as.double(args[7])
#################
pdf(paste(args[2],"-CHG.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[4]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],"-CHG.png",sep=""),res=200, width=460, height=560)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[4]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CHG down!")
######################### CHH
l<-read.table(paste(args[1],".chh",sep=""),sep="\t")
l<-l[,2:ncol(l)]
y<-data.matrix(l)
######CHH max for heatmap
       ceil<- as.double(args[8])
##############
pdf(paste(args[2],"-CHH.pdf",sep=""))
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[5]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
png(paste(args[2],"-CHH.png",sep=""),res=200, width=460, height=560)
pheatmap(pmin(y,ceil), cluster_rows=F,cluster_cols=F, col=colorRampPalette(c("black",args[5]))(32), show_rownames=FALSE,show_colnames=FALSE) #,legend=FALSE)
dev.off()
cat("CHH down!")
