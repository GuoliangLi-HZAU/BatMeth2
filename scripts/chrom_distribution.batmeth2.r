cat("BatMeth2: chrom_distribution\n")
cat("Usage: Rscript InputFile.from.Batmeth2:Split OutFile.pdf\n")
cat("eg:Rscript chrom_distribution.batmeth2.r test.bins.txt chrosome.methy.distri.pdf\n")

Args <- commandArgs()
Infile<-Args[6] #"test.bins.txt"
outFile<-Args[7] #"chrosome.methy.distri.pdf"
n<-read.table(Infile,header=F,sep="\t")
l<-levels(n$V1)
c<-levels(n$V4)
pdf(outFile,width=12,height=6)
par(mfrow=c(2,1), mar=c(4,4,3,2))
for(i in 1:length(l))
{
  for(j in 1:length(c))
  {
    chr <-n[n$V1==l[i] & n$V4==c[j],]
    col=rep("black",nrow(chr))
    col[chr[,3]>=0]<-"indianred2"
    col[chr[,3] < 0]<-"deepskyblue3"
    plot(chr[,2],chr[,3],type="p",col=col,axes=F,xlab="",ylab="",pch=20,cex=0.5)
    #plot(chr[,2],chr[,3],type="p",col=col,axes=F,xlab="",ylab="",ylim=c(-1.2,1.2),pch=20,cex=0.5)
    axis(1, col="black", col.axis="black",line = 1)
    axis(2, col="black", col.axis="black",line = -1)
    mtext("methy.level", side = 2, line = 1)
    text<-paste(l[i],c[j],sep=".")
    mtext(text, side = 4, line = -0.5)
  }
}
dev.off()

png(paste(gsub("pdf","",outFile),"Methylgenome.png", sep=""), width=860, height= (800 * length(l)),res=128)
par(mfrow=c(length(l) * length(c),1), mar=c(4,4,3,2))
for(i in 1:length(l)) ##chr
{
#  png(paste(gsub("pdf","",outFile),l[i], ".png", sep=""), width=860, height=480)
#  par(mfrow=c(3,1), mar=c(4,4,3,2))
  for(j in 1:length(c)) ##context
  {
    chr <-n[n$V1==l[i] & n$V4==c[j],]
    col=rep("black",nrow(chr))
    col[chr[,3]>=0]<-"indianred2"
    col[chr[,3] < 0]<-"deepskyblue3"
    plot(chr[,2],chr[,3],type="p",col=col,axes=F,xlab="",ylab="",pch=20,cex=0.5)
    #plot(chr[,2],chr[,3],type="p",col=col,axes=F,xlab="",ylab="",ylim=c(-1.2,1.2),pch=20,cex=0.5)
    axis(1, col="black", col.axis="black",line = 1)
    axis(2, col="black", col.axis="black",line = -1)
    mtext("methy.level", side = 2, line = 1)
    text<-paste(l[i],c[j],sep=".")
    mtext(text, side = 4, line = -0.5)
  }
#  dev.off()
}
dev.off()
