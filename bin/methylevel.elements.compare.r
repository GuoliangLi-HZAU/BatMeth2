#install.packages("ggplot2")
cat("BatMeth2: methylevel.elements\n")
cat("Usage: Rscript methylevel.elements.r step(default:0.025) Input.Sample1.from.Batmeth2:methyGff Input.Sample2 outfilePrefix xLab1 xLab2 Sample1Prefix Sample2Prefix\n")
cat("eg: Rscript methylevel.elements.r 0.025 gene.meth.Methylevel.1.txt sample2.gene.meth.Methylevel.1.txt methlevel TSS TTS mutant WT\n\n")

library(ggplot2)
library(grid)
options(warn = -1)
Args <- commandArgs(trailingOnly = TRUE)

step<-as.double(Args[1]) #0.025
Infile1<-Args[2] #"gene.meth.Methylevel.1.txt"
Infile2<-Args[3]
OutFile<-Args[4]
start<-Args[5]
end<-Args[6]
sample1<-Args[7]
sample2<-Args[8]

length<-(ceiling(1/step)-1)*3 #117
num<-c(0:(length-1))
lengthp<-length+1
b1<-read.table(Infile1,header=F,sep="\t")
b2<-read.table(Infile2,header=F,sep="\t")
context<-c("CG","CHG","CHH")

#for(line in 1:3)
#{
line=1
pdf(paste(OutFile,".",context[line],".pdf",sep=""),width=10,height=6)
	data<-rbind(data.frame(num,meth=as.double(b1[line,2:lengthp]),Sample=rep(sample1,length)),data.frame(num,meth=as.double(b2[line,2:lengthp]),Sample=rep(sample2,length)) )
	p <- ggplot(data,aes(x=num,y=meth,color=Sample)) 
	max<-max(data$meth)
	max<-max*1.1
	p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+
  theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()
png(paste(OutFile,".",context[line],".png",sep=""),width=860, height=480,res=126)
p <- ggplot(data,aes(x=num,y=meth,color=Sample))
        max<-max(data$meth)
        max<-max*1.1
        p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+
  theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()

#}
line=2
pdf(paste(OutFile,".",context[line],".pdf",sep=""),width=10,height=6)
data<-rbind(data.frame(num,meth=as.double(b1[line,2:lengthp]),Sample=rep(sample1,length)),data.frame(num,meth=as.double(b2[line,2:lengthp]),Sample=rep(sample2,length)) )
p <- ggplot(data,aes(x=num,y=meth,color=Sample)) 
max<-max(data$meth)
max<-max*1.1
p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()

png(paste(OutFile,".",context[line],".png",sep=""),width=860, height=480,res=126)
data<-rbind(data.frame(num,meth=as.double(b1[line,2:lengthp]),Sample=rep(sample1,length)),data.frame(num,meth=as.double(b2[line,2:lengthp]),Sample=rep(sample2,length)) )
p <- ggplot(data,aes(x=num,y=meth,color=Sample))
max<-max(data$meth)
max<-max*1.1
p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()

line=3
pdf(paste(OutFile,".",context[line],".pdf",sep=""),width=10,height=6)
data<-rbind(data.frame(num,meth=as.double(b1[line,2:lengthp]),Sample=rep(sample1,length)),data.frame(num,meth=as.double(b2[line,2:lengthp]),Sample=rep(sample2,length)) )
p <- ggplot(data,aes(x=num,y=meth,color=Sample)) 
max<-max(data$meth)
max<-max*1.1
p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()

png(paste(OutFile,".",context[line],".png",sep=""),width=860, height=480,res=126)
data<-rbind(data.frame(num,meth=as.double(b1[line,2:lengthp]),Sample=rep(sample1,length)),data.frame(num,meth=as.double(b2[line,2:lengthp]),Sample=rep(sample2,length)) )
p <- ggplot(data,aes(x=num,y=meth,color=Sample))
max<-max(data$meth)
max<-max*1.1
p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(xintercept=c(length/3,length/3*2),linetype="dotted",color="gray",size=1)+theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()
