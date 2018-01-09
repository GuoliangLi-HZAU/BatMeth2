#install.packages("ggplot2")
cat("BatMeth2: methylevel.elements\n")
cat("Usage: Rscript methylevel.elements.r step(default:0.025) Input.from.Batmeth2:methyGff outfile.pdf xLab1 xLab2\n")
cat("eg: Rscript methylevel.elements.r 0.025 gene.meth.Methylevel.1.txt methlevel.pdf TSS TTS\n\n")

library(ggplot2)
library(grid)
options(warn = -1)
Args <- commandArgs(trailingOnly = TRUE)

step<-as.double(Args[1]) #0.025
Infile<-Args[2] #"gene.meth.Methylevel.1.txt"
OutFile<-Args[3]
start<-Args[4]
end<-Args[5]

length<-(ceiling(1/step)-1)*3 #117
num<-c(0:(length-1))
lengthp<-length+1
b<-read.table(Infile,header=F,sep="\t")
data<-rbind(data.frame(num,meth=as.double(b[1,2:lengthp]),Context=rep(" CG",length)),data.frame(num,meth=as.double(b[2,2:lengthp]),Context=rep(" CHG",length)),data.frame(num,meth=as.double(b[3,2:lengthp]),Context=rep(" CHH",length)) )
pdf(OutFile,width=10,height=6)
max<-max(data$meth)
max<-max*1.1
p <- ggplot(data,aes(x=num,y=meth,color=Context)) #rbind(data[data$Context==" CG",],data[data$Context==" CHG",],data[data$Context==" CHH",])
p + geom_line(lwd=0.70)+xlab("")+theme_bw()+theme(panel.grid=element_blank())+scale_colour_manual(values = c(" CG" = "black"," CHG" = "purple"," CHH"="blue"))+scale_x_continuous("",limits=c(0, length-1),breaks=c(0,length/3,length/3*2,length-1),labels=c("up",start,end,"down"))+geom_vline(x=c(length/3,length/3*2),linetype="dotted",color="blue",size=1)+
  theme(legend.key=element_rect(linetype='dashed',color="white"),axis.text.y = element_text(size=13),axis.text.x = element_text(size=18),legend.title = element_text(size=16),legend.text = element_text(size=15),legend.key.height=unit(1.2,'cm')) +ylim(0,max)
dev.off()

