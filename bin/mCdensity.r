cat("BatMeth2: mC\n")
cat("Usage: Rscript thisRfile Inputfile1(mCdensity) OutFile1.pdf Inputfile2(mCcatero) OutFile2.pdf\n")

Args <- commandArgs()
Infile<-Args[6] #"IMR90.mCdensity.txt"
outFile<-Args[7] #"mCdensity.pdf"
a<-read.table(Infile,row.names=1)

###pdf
pdf(outFile, width=12,height=6)
par(mfrow=c(3,1), mar=c(4,4,3,2))
plot(x=seq(1:100), y =a[1,], type="h", lwd = 10, lend = 2, col="red", xlab = "", ylab = "number", main = "mCG")
plot(x=seq(1:100), y =a[2,], type="h", lwd = 10, lend = 2, col="green", xlab = "", ylab = "number", main = "mCHG")
plot(x=seq(1:100), y =a[3,], type="h", lwd = 10, lend = 2, col="black", xlab = "DNA methylation level", ylab = "number", main = "mCHH")
dev.off()

###png
png(paste(gsub("pdf","",outFile),"png", sep=""), width=660, height= (200 * 3),res=96)
par(mfrow=c(3,1), mar=c(4,4,3,2))
plot(x=seq(1:100), y =a[1,], type="h", lwd = 10, lend = 2, col="red", xlab = "", ylab = "number", main = "mCG")
plot(x=seq(1:100), y =a[2,], type="h", lwd = 10, lend = 2, col="green", xlab = "", ylab = "number", main = "mCHG")
plot(x=seq(1:100), y =a[3,], type="h", lwd = 10, lend = 2, col="black", xlab = "DNA methylation level", ylab = "number", main = "mCHH")
dev.off()

###file2
Infile2<-Args[8] #"IMR90.mCdensity.txt"
outFile2<-Args[9] #"mCdensity.pdf"
a<-read.table(Infile2,row.names=1)

###pdf
pdf(outFile2, width=12, height=6)
par(mfrow=c(2,1), mar=c(4,4,3,2))
plot(x=seq(1:5), y=head(a, n=5L)[,1], xaxt="n",type="h", lwd = 20, lend = 2, xlab = "", ylab = "number", main = "Total Cytosines")
axis(1,at=c(1,2,3,4,5),labels=c("M","Mh","H","hU","U"))
plot(x=seq(1:5), y=tail(a, n=5L)[,1], xaxt="n",type="h", lwd = 20, lend = 2, xlab = "DNA methylation category", ylab = "number", main = "mCG")
axis(1,at=c(1,2,3,4,5),labels=c("CpG_M","CpG_Mh","CpG_H","CpG_hU","CpG_U"))
dev.off()

###png
png(paste(gsub("pdf","",outFile2),"png", sep=""), width=660, height= (200 * 2),res=96)
par(mfrow=c(2,1), mar=c(4,4,3,2))

plot(x=seq(1:5), y=head(a, n=5L)[,1], xaxt="n",type="h", lwd = 20, lend = 2, xlab = "", ylab = "number", main = "Total Cytosines")
axis(1,at=c(1,2,3,4,5),labels=c("M","Mh","H","hU","U"))

plot(x=seq(1:5), y=tail(a, n=5L)[,1], xaxt="n",type="h", lwd = 20, lend = 2, xlab = "DNA methylation category", ylab = "number", main = "mCG")
axis(1,at=c(1,2,3,4,5),labels=c("CpG_M","CpG_Mh","CpG_H","CpG_hU","CpG_U"))

dev.off()
