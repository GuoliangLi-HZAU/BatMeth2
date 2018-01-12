##Input OutPut
options(warn = -1)
Args <- commandArgs(trailingOnly = TRUE)
Infile<-Args[1]
a<- read.table(Infile,header=F)
symbol<-as.matrix(a[a$V3=="DMC",]$V1)

outfile<-Args[2]
pdf(outfile)

totalDMC<-a[a$V3=="totalDMC",]$V2
dmc<-as.matrix(a[a$V3=="DMC",]$V2)
ydmc <- matrix(c(dmc), ncol=1, nrow=nrow(dmc), byrow=FALSE) 
rownames(ydmc) <-  symbol[,1]
bp <- barplot(t(ydmc), horiz=F, beside=FALSE,  las=3,ylab="DMC") 
#text(bp,t(dmc), signif(dmc/totalDMC,2), pos=3)

totalC<-a[a$V3=="totalC",]$V2
allC<-as.matrix(a[a$V3=="AllC",]$V2)
yallC <- matrix(c(allC), ncol=1, nrow=nrow(allC), byrow=FALSE)
rownames(yallC) <-  symbol[,1]
bp <- barplot(t(yallC), horiz=F, beside=FALSE,  las=3,ylab="AllC")
#text(bp,t(allC), signif(allC/totalC,2), pos=3)

dev.off()

outpng<-paste(gsub("pdf","",outfile),".png", sep="")
png(outpng,res=128)

totalDMC<-a[a$V3=="totalDMC",]$V2
dmc<-as.matrix(a[a$V3=="DMC",]$V2)
ydmc <- matrix(c(dmc), ncol=1, nrow=nrow(dmc), byrow=FALSE)
rownames(ydmc) <-  symbol[,1]
bp <- barplot(t(ydmc), horiz=F, beside=FALSE,  las=3,ylab="DMC")
#text(bp,t(dmc), signif(dmc/totalDMC,2), pos=3)

totalC<-a[a$V3=="totalC",]$V2
allC<-as.matrix(a[a$V3=="AllC",]$V2)
yallC <- matrix(c(allC), ncol=1, nrow=nrow(allC), byrow=FALSE)
rownames(yallC) <-  symbol[,1]
bp <- barplot(t(yallC), horiz=F, beside=FALSE,  las=3,ylab="AllC")
#text(bp,t(allC), signif(allC/totalC,2), pos=3)

dev.off()
