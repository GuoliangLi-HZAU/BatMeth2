cat("BatMeth2: density_plot_with_methyl.r\n")
cat("Usage: Rscript density_plot_with_methyl.r inputFile1 input2 genedensityFile TEdensity output.pdf label1 label2\n")
cat("eg: Rscript density_plot_with_methyl.r Cr_DJ.bins.strand.aver.txt ../S706osdrm2/S706.bins.strand.aver.txt Cr_DJ.geneBody.count.C.gffDensity.1.txt Cr_DJ.TE.count.C.gffDensity.1.txt density.pdf Cr_DJ OsDrm2\n\n")

Args <- commandArgs()
Infile1<-Args[6]
Infile2<-Args[7]
genefile<-Args[8]
TEfile<-Args[9]
outPDF<-Args[10]
outPDFf<-paste(gsub("pdf","",outPDF),"final.pdf", sep="")
outpng<-paste(gsub("pdf","",outPDF),".png", sep="")
outpngf<-paste(gsub("pdf","",outPDF),"final.png", sep="")
label1<-Args[11]
label2<-Args[12]
#install.packages("ggplot2")
library(ggplot2)
library(grid)
Methylp <- read.table(Infile1,sep="\t") #"Cr_DJ.bins.strand.aver.txt"
colnames(Methylp) <- c('chr', 'pos', 'Methyl', 'context')
Methyl<-Methylp[Methylp$Methyl>=0,]

Methyl2p <- read.table(Infile2,sep="\t")  #"S706.bins.strand.aver.txt"
colnames(Methyl2p) <- c('chr', 'pos', 'Methyl', 'context')
Methyl2<-Methyl2p[Methyl2p$Methyl>=0,]
Methyl2[,3] <- Methyl2[,3]+1
################################## CG
chrMethyl<-Methyl[Methyl$context=="CG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
         lstMeth <- max(chrMethyl[ndx, 2])
         if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
         if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

p <- ggplot() + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CG",chr) )),method="auto",linetype=1)  +
  theme_bw(base_size=15)+theme(legend.position='none')

#p <- ggplot() #chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(chr))
#(p2 <- p  + stat_smooth(se=F,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(chr)),method="auto",linetype=1) )# + geom_point()
           #(p2 <- p + geom_line(aes(x=pos, y=Methyl, size=3.5, colour=as.factor(chr)), alpha=1/3))
#(p3 <- p2 )#+ scale_color_manual(values=rep(c('black', 'dark green'), 6)) )
#(p4 <- p3 + theme_bw(base_size=15)) 
#(p5 <- p4 + theme(legend.position='none'))  
#(p6 <- p5 + scale_x_continuous(labels=as.character(chr[1:chrNum]), breaks=bpMidVec))   
#(p7 <- p6 + geom_vline(x=bpMidVeX, linetype=2, col='gray', lwd=0.5)) 
         #(p8 <- p7 + ggtitle('Chromsome Methylation Level') + xlab('') + ylab('Methylation.Level')) + theme(panel.grid=element_blank())

##############################CHG
chrMethyl<-Methyl[Methyl$context=="CHG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
                     lstMeth <- max(chrMethyl[ndx, 2])
                     if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
                     if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

(p1 <- p + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CHG",chr) )),method="auto",linetype=2)  +
  theme_bw(base_size=15)+theme(legend.position='none') )
 
##############################CHH
chrMethyl<-Methyl[Methyl$context=="CHH",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
                     lstMeth <- max(chrMethyl[ndx, 2])
                     if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
                     if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

#p <- ggplot() #chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(chr))
(p2 <- p1 + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CHH",chr) )),method="auto",linetype=3)  +
   theme_bw(base_size=15)+theme(legend.position='none') )
 
################################file2 CG
chrMethyl<-Methyl2[Methyl2$context=="CG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
                     lstMeth <- max(chrMethyl[ndx, 2])
                     if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
                     if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}


(p3 <- p2 + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CG",chr) )),method="auto",linetype=1)  +
   theme_bw(base_size=15)+theme(legend.position='none') )
 
################################file2 CHG
chrMethyl<-Methyl2[Methyl2$context=="CHG",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
                     lstMeth <- max(chrMethyl[ndx, 2])
                     if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
                     if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}


(p4 <- p3 + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CHG",chr) )),method="auto",linetype=2)  +
   theme_bw(base_size=15)+theme(legend.position='none') )
 
################################file2 CHH
chrMethyl<-Methyl2[Methyl2$context=="CHH",]
chr<-levels(as.factor(chrMethyl$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrMethyl[, 1]==chr[i] )
                     lstMeth <- max(chrMethyl[ndx, 2])
                     if (i < chrNum) ndx2 <- which(chrMethyl[, 1]== chr[i+1] )
                     if (i < chrNum) chrMethyl[ndx2, 2] <- chrMethyl[ndx2, 2] + lstMeth
}

(p5 <- p4 + stat_smooth(se=F,size=1,data=chrMethyl,aes(x=pos, y=Methyl,colour=as.factor(gsub("Chr|chr","CHH",chr)) ),method="auto",linetype=3)  +
   theme_bw(base_size=15)+theme(legend.position='none') )
 
####################gene density

density <- read.table(genefile,sep="\t")  #"Cr_DJ.geneBody.count.C.gffDensity.1.txt"
colnames(density) <- c('chr', 'pos', 'density', 'strand')
density[,3] <- density[,3]+2

chrDensity<-density[density$strand=="+-",]
chr<-levels(as.factor(chrDensity$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrDensity[, 1]==chr[i] )
                     lstMeth <- max(chrDensity[ndx, 2])
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

(p6 <- p5 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,colour=as.factor(gsub("Chr|chr","gene",chr)) ),method="auto",linetype=1)+
   theme_bw(base_size=15)+theme(legend.position='none') )

##########################TE gff density
density <- read.table(TEfile,sep="\t")  #"Cr_DJ.TE.count.C.gffDensity.1.txt"
colnames(density) <- c('chr', 'pos', 'density', 'strand')
density[,3] <- density[,3]+2

chrDensity<-density[density$strand=="+-",]
chr<-levels(as.factor(chrDensity$chr))
chrNum=length(chr)
for (i in 1:chrNum){ ndx <- which(chrDensity[, 1]==chr[i] )
                     lstMeth <- max(chrDensity[ndx, 2])
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
(p7 <- p6 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,colour=as.factor(gsub("Chr|chr","TE",chr)) ),method="auto",linetype=1)+
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3),labels=c("0","0.5","0","0.5","0","0.5","1") )+
   annotate("text", label = c(label1,label2,"Density","CHH","CHG","CG"), x = c(100,100,100,bpMidVeX[length(bpMidVeX)],bpMidVeX[length(bpMidVeX)],bpMidVeX[length(bpMidVeX)]), y = c(0.8,1.8,2.8,0.06,0.35,0.6), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,1,2), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.margin=unit(0,"line") )
   )
dev.off()
png(outpng,width=860, height=480,res=128)
(p7 <- p6 + geom_smooth(se=F,size=1,data=chrDensity,aes(x=pos, y=density,colour=as.factor(gsub("Chr|chr","TE",chr)) ),method="auto",linetype=1)+
   theme_bw(base_size=15)+theme(legend.position='none')+
   scale_x_continuous(labels=as.character(gsub("Chr|chr","",chr[1:chrNum]) ), breaks=bpMidVec)+
   scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3),labels=c("0","0.5","0","0.5","0","0.5","1") )+
   annotate("text", label = c(label1,label2,"Density","CHH","CHG","CG"), x = c(100,100,100,bpMidVeX[length(bpMidVeX)],bpMidVeX[length(bpMidVeX)],bpMidVeX[length(bpMidVeX)]), y = c(0.8,1.8,2.8,0.06,0.35,0.6), size = 4, colour = "azure4")+
   geom_vline(xintercept=bpMidVeX, linetype=2, col='gray', lwd=0.5)+geom_hline(yintercept=c(0,1,2), linetype=5, col='gray', lwd=0.2,alpha=0.6)+
   ggtitle('Chromsome') + xlab('') + ylab('Methylation.Level') + theme(panel.grid=element_blank(),panel.margin=unit(0,"line") )
   )
dev.off()
############################legend
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
       panel.background=element_blank(), 
       axis.text.x=element_blank(), axis.text.y=element_blank(),           
       axis.title.x=element_blank(), axis.title.y=element_blank())

bottom<-ggplot()+geom_line(aes(x=c(0,0.02),y=c(2,2)),linetype=1,colour="purple" ,size=1)+annotate("text",label="TE density",x=0.04,y=2,size=4,color="azure4")+
  geom_line(aes(x=c(0.07,0.09),y=c(2,2)),linetype=1,colour="deepskyblue",size=1 )+annotate("text",label="Gene density",x=0.12,y=2,size=4,color="azure4")+ylim(2,2)+xlim(0,0.29)+
  geom_line(aes(x=c(0.15,0.17),y=c(2,2)),linetype=1,colour="indianred1" ,size=1)+annotate("text",label="CG",x=0.18,y=2,size=4,color="azure4")+
  geom_line(aes(x=c(0.20,0.22),y=c(2,2)),linetype=2,colour="khaki4" ,size=1)+annotate("text",label="CHG",x=0.23,y=2,size=4,color="azure4")+
  geom_line(aes(x=c(0.25,0.27),y=c(2,2)),linetype=3,colour="cyan1" ,size=1)+annotate("text",label="CHH",x=0.28,y=2,size=4,color="azure4")+
  theme(panel.grid=element_blank())+ xlab('') + ylab('')+theme_bw()+theme(panel.grid=element_blank(),panel.border = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.margin=unit(0,"line"))
pdf(outPDFf,width=12,height=7)
library(gridExtra)
grid.arrange(p7,empty, bottom, empty, ncol=2, nrow=2, widths=c(1, 0.06),heights=c(7.6,1) ,padding= unit(0, "line"))
dev.off()
png(outpngf,width=860, height=580,res=96)
grid.arrange(p7,empty, bottom, empty, ncol=2, nrow=2, widths=c(1, 0.06),heights=c(7.6,1) ,padding= unit(0, "line"))
dev.off()

