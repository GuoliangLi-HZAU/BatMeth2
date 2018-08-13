cat("BatMeth2: density_plot_with_methyl.r\n")
cat("Usage: Rscript density_plot_with_methyl.r inputFile1 geneDensityFile output.pdf label1\n")
cat("eg: Rscript density_plot_with_methyl.r Cr_DJ.bins.strand.aver.txt Cr_DJ.geneBody.count.C.gffDensity.1.txt density.pdf Cr_DJ \n\n")

Args <- commandArgs()
Infile1<-Args[6]
#Infile2<-Args[7]
genefile<-Args[7]
outPDF<-Args[8]
outPDFf<-paste(gsub(".pdf","",outPDF),"plot.pdf", sep="")
outpng<-paste(gsub("pdf","",outPDF),"png", sep="")
outpngf<-paste(gsub(".pdf","",outPDF),"plot.png", sep="")
label1<-Args[9]
#label2<-Args[12]
#install.packages("ggplot2")
library(ggplot2)
library(grid)
Methylp <- read.table(Infile1,sep="\t") #"Cr_DJ.bins.strand.aver.txt"
colnames(Methylp) <- c('chr', 'pos', 'Methyl', 'context')
Methyl<-Methylp[Methylp$Methyl>=0,]

#Methyl2p <- read.table(Infile2,sep="\t")  #"S706.bins.strand.aver.txt"
#colnames(Methyl2p) <- c('chr', 'pos', 'Methyl', 'context')
#Methyl2<-Methyl2p[Methyl2p$Methyl>=0,]
#Methyl2[,3] <- Methyl2[,3]+1
maxmethlevel<-max(Methyl$Methyl) + 0.1
if (maxmethlevel>1) maxmethlevel <- 1
################################## CG
chrMethyl<-Methyl[Methyl$context=="CG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
maxMeth=1:chrNum
maxMethOri=1:chrNum
POINT=TRUE
for (i in 1:chrNum){ 
  ndx <- which(chrMethyl[, 1]==chr[i] )
  lstMeth <- max(chrMethyl[ndx, 2])
  if(i > 1) maxMeth[i] <- (lstMeth - maxMethOri[i-1])
  else maxMeth[i] <- lstMeth
  maxMethOri[i] <- lstMeth
  if(maxMeth[i]>200) POINT=FALSE
  if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
  if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}
#as.factor(gsub("Chr|chr","CG",chr) )

if(POINT){
	p <- ggplot(data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CG",chr) ))) + geom_point(colour="indianred1") + theme_bw(base_size=15)+theme(legend.position='none')
}else{
	p <- ggplot() + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CG",chr) )), colour="indianred1", method="auto", linetype=1) + theme_bw(base_size=15)+theme(legend.position='none')
}
##############################CHG
chrMethyl<-Methyl[Methyl$context=="CHG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
maxMeth=1:chrNum
maxMethOri=1:chrNum
POINT=TRUE
for (i in 1:chrNum){ 
  ndx <- which(chrMethyl[, 1]==chr[i] )
  lstMeth <- max(chrMethyl[ndx, 2])
  if(i > 1) maxMeth[i] <- (lstMeth - maxMethOri[i-1])
  else maxMeth[i] <- lstMeth
  maxMethOri[i] <- lstMeth
  if(maxMeth[i]>200) POINT=FALSE
  if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
  if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

if(POINT) {
	p1<- p + geom_point(data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CHG",chr) )), colour="khaki4") + theme_bw(base_size=15)+theme(legend.position='none')
}else{
	(p1 <- p + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CHG",chr) )) , colour="khaki4" ,method="auto",linetype=2)  + theme_bw(base_size=15)+theme(legend.position='none') )
}
##############################CHH
chrMethyl<-Methyl[Methyl$context=="CHH",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
maxMeth=1:chrNum
maxMethOri=1:chrNum
POINT=TRUE
for (i in 1:chrNum){ 
  ndx <- which(chrMethyl[, 1]==chr[i] )
  lstMeth <- max(chrMethyl[ndx, 2])
  if(i > 1) maxMeth[i] <- (lstMeth - maxMethOri[i-1])
  else maxMeth[i] <- lstMeth
  maxMethOri[i] <- lstMeth
  if(maxMeth[i]>200) POINT=FALSE
  if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
  if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

if(POINT){
	p2<- p1 + geom_point(data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CHH",chr) )), colour="cyan1") + theme_bw(base_size=15)+theme(legend.position='none')
}else{
	(p2 <- p1 + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,group=as.factor(gsub("Chr|chr","CHH",chr) )) , colour="cyan1",method="auto",linetype=3)  + theme_bw(base_size=15)+theme(legend.position='none') )
}
####################gene density

density <- read.table(genefile,sep="\t")  #"Cr_DJ.geneBody.count.C.gffDensity.1.txt"
colnames(density) <- c('chr', 'pos', 'density', 'strand')
density[,3] <- density[,3]+maxmethlevel #1

chrDensity<-density[density$strand=="+-",]
chr<-levels(as.factor(chrDensity$chr))
chrNum=length(chr)
maxMeth=1:chrNum
maxMethOri=1:chrNum
POINT=TRUE
for (i in 1:chrNum){ 
  ndx <- which(chrDensity[, 1]==chr[i] )
  lstMeth <- max(chrDensity[ndx, 2])
  if(i > 1) maxMeth[i] <- (lstMeth - maxMethOri[i-1])
  else maxMeth[i] <- lstMeth
  maxMethOri[i] <- lstMeth
  if(maxMeth[i]>200) POINT=FALSE
  if (i < chrNum) ndx2 <- which(chrDensity[, 1]== chr[i+1] )
  if (i < chrNum) chrDensity[ndx2, 2] <- chrDensity[ndx2, 2] + lstMeth
}

bpMidVec <- vector(length=chrNum)
bpMidVeX <- vector(length=chrNum)

for (i in 1:chrNum){ndx <- which(chrDensity[, 1]==chr[i] )
                    posSub <- chrDensity[ndx, 2]
                    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub) #
                    if((max(posSub) - min(posSub)) < 100){chr[i]=""}
                    bpMidVeX[i] <- max(posSub)
}

if(POINT) {
	p6<- p2 + geom_point(data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","gene",chr))),colour="deepskyblue") + theme_bw(base_size=15)+theme(legend.position='none')
}else{
	(p6 <- p2 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","gene",chr)) ) ,colour="deepskyblue" ,method="auto",linetype=1)+ theme_bw(base_size=15)+theme(legend.position='none') )
}
##########################TE gff density
density <- read.table(TEfile,sep="\t")  #"Cr_DJ.TE.count.C.gffDensity.1.txt"
colnames(density) <- c('chr', 'pos', 'density', 'strand')
density[,3] <- density[,3]+maxmethlevel #1

chrDensity<-density[density$strand=="+-",]
chr<-levels(as.factor(chrDensity$chr))
chrNum=length(chr)
maxMeth=1:chrNum
maxMethOri=1:chrNum
POINT=TRUE
for (i in 1:chrNum){
  ndx <- which(chrDensity[, 1]==chr[i] )
  lstMeth <- max(chrDensity[ndx, 2])
  if(i > 1) maxMeth[i] <- (lstMeth - maxMethOri[i-1])
  else maxMeth[i] <- lstMeth
  maxMethOri[i] <- lstMeth
  if(maxMeth[i]>200) POINT=FALSE
  if (i < chrNum) ndx2 <- which(chrDensity[, 1]== chr[i+1] )
  if (i < chrNum) chrDensity[ndx2, 2] <- chrDensity[ndx2, 2] + lstMeth
}

bpMidVec <- vector(length=chrNum)
bpMidVeX <- vector(length=chrNum)

for (i in 1:chrNum){ndx <- which(chrDensity[, 1]==chr[i] )
                    posSub <- chrDensity[ndx, 2]
                    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub) #
                    if((max(posSub) - min(posSub)) < 100){chr[i]=""}
                    bpMidVeX[i] <- max(posSub)
}
pdf(outPDF,width=12,height=7)
if(POINT){
	(p7 <- p6 + geom_point(data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","TE",chr))),colour="purple") +
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,maxmethlevel,maxmethlevel+0.1, maxmethlevel + 0.5, maxmethlevel +1),labels=c("0",maxmethlevel,"0","0.5","1") )+
   annotate("text", label = c(label1,"Density"), x = c(bpMidVeX[length(bpMidVeX)]+1,bpMidVeX[length(bpMidVeX)]+1), y = c(maxmethlevel/2,maxmethlevel+ 0.8), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,maxmethlevel,maxmethlevel+1), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.spacing=unit(0,"line") )
   )

}else{
	(p7 <- p6 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","TE",chr))) ,colour="purple",method="auto",linetype=1)+
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,maxmethlevel,maxmethlevel+0.1, maxmethlevel + 0.5, maxmethlevel +1),labels=c("0",maxmethlevel,"0","0.5","1") )+
   annotate("text", label = c(label1,"Density"), x = c(bpMidVeX[length(bpMidVeX)]+1,bpMidVeX[length(bpMidVeX)]+1), y = c(maxmethlevel/2,maxmethlevel+ 0.8), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,maxmethlevel,maxmethlevel+1), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.spacing=unit(0,"line") )
   )
}
dev.off()
png(outpng,width=860, height=480,res=128)
if(POINT){
	(p7 <- p6 + geom_point(data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","TE",chr))),colour="purple") +
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,maxmethlevel,maxmethlevel+0.1, maxmethlevel + 0.5, maxmethlevel +1),labels=c("0",maxmethlevel,"0","0.5","1") )+
   annotate("text", label = c(label1,"Density"), x = c(bpMidVeX[length(bpMidVeX)]+1,bpMidVeX[length(bpMidVeX)]+1), y = c(maxmethlevel/2,maxmethlevel+ 0.8), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,maxmethlevel,maxmethlevel+1), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.spacing=unit(0,"line") )
   )
}else{
	(p7 <- p6 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,group=as.factor(gsub("Chr|chr","TE",chr))) ,colour="purple",method="auto",linetype=1)+
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,maxmethlevel,maxmethlevel+0.1, maxmethlevel + 0.5, maxmethlevel +1),labels=c("0",maxmethlevel,"0","0.5","1") )+
   annotate("text", label = c(label1,"Density"), x = c(bpMidVeX[length(bpMidVeX)]+1,bpMidVeX[length(bpMidVeX)]+1), y = c(maxmethlevel/2,maxmethlevel+ 0.8), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,maxmethlevel,maxmethlevel+1), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.spacing=unit(0,"line") )
   )
}
dev.off()
############################legend
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
       panel.background=element_blank(), 
       axis.text.x=element_blank(), axis.text.y=element_blank(),           
       axis.title.x=element_blank(), axis.title.y=element_blank())

bottom<-ggplot()+ geom_line(aes(x=c(0,0.02),y=c(2,2)),linetype=1,colour="purple" ,size=1)+ annotate("text",label="TE density",x=0.04,y=2,size=4,color="azure4")+
  geom_line(aes(x=c(0.07,0.09),y=c(2,2)),linetype=1,colour="deepskyblue",size=1 )+annotate("text",label="Gene density",x=0.12,y=2,size=4,color="azure4")+ylim(2,2)+xlim(0,0.29)+
  geom_line(aes(x=c(0.15,0.17),y=c(2,2)),linetype=1,colour="indianred1" ,size=1)+annotate("text",label="CG",x=0.18,y=2,size=4,color="indianred1")+
  geom_line(aes(x=c(0.20,0.22),y=c(2,2)),linetype=2,colour="khaki4" ,size=1)+annotate("text",label="CHG",x=0.23,y=2,size=4,color="azure4")+
  geom_line(aes(x=c(0.25,0.27),y=c(2,2)),linetype=3,colour="cyan1" ,size=1)+annotate("text",label="CHH",x=0.28,y=2,size=4,color="cyan1")+
  theme(panel.grid=element_blank())+ xlab('') + ylab('')+theme_bw()+theme(panel.grid=element_blank(),panel.border = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=unit(0,"line"))
pdf(outPDFf,width=12,height=7)
library(gridExtra)
grid.arrange(p7,empty, bottom, empty, ncol=2, nrow=2, widths=c(1, 0.06),heights=c(7.6,1) ,padding= unit(0, "line"))
dev.off()
png(outpngf,width=860, height=580,res=96)
grid.arrange(p7,empty, bottom, empty, ncol=2, nrow=2, widths=c(1, 0.06),heights=c(7.6,1) ,padding= unit(0, "line"))
dev.off()
