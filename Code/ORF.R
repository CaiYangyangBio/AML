library(ggplot2)
library(ggpubr)
library(forestplot)
#======MAP_gene
immune<-read.table("~/Desktop/22白血病/Runtime/1预处理/All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

#======MAP
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",check.names = F,fill = T)
MAP<-MAP[MAP$Type!="Canonical",]
MAP<-merge(MAP,immune)

MAPRegion<-MAP[,c("Leading razor protein","PlotType")]
MAPRegion<-MAPRegion[!duplicated(MAPRegion),]

MAPRegionNum<-data.frame(table(MAPRegion$PlotType))
MAPRegionNum<-MAPRegionNum[order(MAPRegionNum$Freq),]
Region<-as.character(MAPRegionNum$Var1)
MAPRegionNum$Var1<-factor(MAPRegionNum$Var1,levels = Region)

color<-c("lncRNA"="#81BC7C","3' Overlap dORF"="#FFB9B4",
         "5' uORF"="#D8E289","3' dORF"="#FCAE3F","5' Overlap uORF"="#C3B5E6",
         "Pseudogene"="#FADD5B","Out-of-Frame"="#FFDDB8","Other"="#7466A3")

#==================ORF数量========================#
g1<-ggplot(MAPRegionNum,aes(x = log2(Freq+1),y = Var1,fill = Var1))+
  geom_bar(stat="identity")+
  xlab("log2( MS detected proteins + 1)") +
  ylab("")+
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 12,color="black"),
        axis.title = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "none",
        legend.title = element_blank())+
  scale_fill_manual(values = color)+
  geom_text(aes(label=Freq),size=4,vjust=0.5,position = position_stack(1))+
  scale_x_continuous(limits = c(0,10),expand = c(0,0),breaks = seq(0,10,2))


#==================参考库========================#
non_canonical<-read.table("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/1-SeqName/5nuORFRename.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,fill = T,quote = "")
ref<-data.frame(table(non_canonical$plotType))
ref$Var1<-factor(ref$Var1,levels = Region)
g2<-ggpie(ref,"Freq",lab.pos ="in",label="Freq",
          fill="Var1",color="white",palette=color)+
  theme(legend.position = "none")

#==================质谱识别========================#
ms<-MAP[,c("Epitopes","PlotType")]
ms<-ms[!duplicated(ms),]
ms<-data.frame(table(ms$PlotType))
index<-match(rev(Region),ms$Var1,nomatch = 0)
ms<-ms[index,]
index<-match(ms$Var1,ref$Var1)
data<-data.frame(Region=ms$Var1,ms=ms$Freq,ref=ref$Freq[index])

ms_ref<-data.frame()
for (i in 1:nrow(data)) {
  tt<-matrix(0,ncol = 2,nrow = 2)
  tt[,1]<-as.numeric(data[i,2:3])
  tt[,2]<-c(sum(data$ms[-i]),sum(data$ref[-i]))
  fish<-fisher.test(tt)
  OR<-round(fish$estimate,2)
  low<-round(fish$conf.int[1],2)
  high<-round(fish$conf.int[2],2)
  p<-format(fish$p.value,digits = 3)
  ms_ref<-rbind(ms_ref,c(as.character(data[i,1]),OR,low,high,p))
}
colnames(ms_ref)<-c("ORF Types","OR","low","high","P values")
ms_ref$CI<-paste0(ms_ref$OR,"(",ms_ref$low,"~",ms_ref$high,")")
ms_ref<-rbind(c("ORF Types",NA,NA,NA,"P values","OR(95% CI)"),ms_ref)

setwd("~/Desktop/22白血病/Runtime/2特征比较/9.ORF")
pdf("2C.pdf",width = 6,height = 5,onefile = F)
forestplot(as.matrix(ms_ref[,c(1,6,5)]),  #显示的文本
           mean=as.numeric(ms_ref$OR),#图形HR/OR部分
           lower=as.numeric(ms_ref$low),#95%CI下限
           upper=as.numeric(ms_ref$high),#95%CI上限
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           xticks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize =0.11, ##误差条中的圆心点大小
           col=fpColors(line = "#79976B", #误差条的线的颜色
                        box="#FFB511"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 2,   # 误差条的线的宽度
           ci.vertices.height = 0.07, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), xlab = gpar(cex = 0.7), cex = 0.7), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="")
dev.off()

#==================样本MAP频率========================#
MAP_num<-data.frame(table(MAP$Samples,MAP$PlotType))
MAP_num$Var2<-factor(MAP_num$Var2,levels = Region)
MAP_num$Var1<-factor(MAP_num$Var1,levels = paste0("AML",1:19))

g3<-ggplot(MAP_num,aes(x = Var1,y = Freq,fill = Var2))+
  geom_bar(stat="identity",position = "fill")+
  xlab("") +
  ylab("Fraction of ncMAPs across Different Source ORF Types")+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 12,color="black",angle = 90),
        axis.text.y = element_text(size = 12,color="black"),
        axis.title = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = color)

setwd("~/Desktop/22白血病/Runtime/2特征比较/9.ORF")
pdf("2A.pdf",width = 6,height = 5)
g1
dev.off()

pdf("2B.pdf",width = 6,height = 5)
g2
dev.off()

pdf("S2A.pdf",width = 10,height = 5)
g3
dev.off()

setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",check.names = F,fill = T)
aa<-MAP[,c("Epitopes","Type")]
aa<-aa[!duplicated(aa),]
aa<-data.frame(table(aa$Type))
aa$MHC<-"MHC"
aa$Var1<-factor(aa$Var1,levels = c("Non-Canonical","Canonical"))

setwd("~/Desktop/22白血病/Runtime/2特征比较/9.ORF")
pdf("bar.pdf",width = 4,height = 0.5)
ggplot(aa,aes(x = Freq,y = MHC,fill = Var1))+
  geom_bar(stat="identity",position = "fill")+
  xlab("") +
  ylab("")+
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"))
dev.off()
