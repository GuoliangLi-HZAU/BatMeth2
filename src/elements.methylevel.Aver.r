cat("Batmeth2: elements.methylevel.Aver\n")
cat("Usage:Rscript elements.methylevel.Aver.r Input.from.Batmeth2:methyGff outfile.pdf\n")
cat("eg: Rscript elements.methylevel.Aver.r gene.meth.AverMethylevel.1.txt elements.pdf\n\n")

#install.packages("ggplot2")
library(ggplot2)
Args <- commandArgs()
Infile<-Args[6] #:gene.meth.AverMethylevel.1.txt
outFile<-Args[7]
a<-read.table(Infile,header=T,sep="\t")
pdf(outFile)
ggplot(a)+geom_boxplot(aes(x=Context, y=DNA.methylation.level,color=Regions))+ scale_fill_discrete(breaks=c("UP","BODY","DOWN"))
dev.off()
