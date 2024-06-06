library(ggplot2)
library(ggpubr)
library(motifStack)
library(HDMD)
library(ecodist)
library(dbscan)
library(fpc)
library(umap)
library(reshape2)
library(factoextra)
library(ggseqlogo)
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif")
data<-read.table("motif_umap.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

g1<-ggplot(data = data,aes(x = UMAP_1, y = UMAP_2,color = HLA,shape = HLA,size=Freq)) + 
  geom_point()+
  scale_shape_manual(values = c(1,1,1))+
  scale_color_manual(values = c("HLAA"="#F49B0C","HLAB"="#A3C78E","HLAC"="#F6D780"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12,color="black"),
        axis.text = element_text(size =12,color="black"),
        axis.title = element_text(size = 14),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.7),
        legend.position = "top",
        legend.title = element_blank())

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
pdf("Fig4B.pdf",height = 5,width = 6)
g1
dev.off()


setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

#--------------xLxxxxxxV----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
LV<-data[which(data$newcluster=="Cluster3"),]
LV<-LV[1:18,]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",LV$Allele,"nmds.txt")
out<-paste0(LV$Allele)
g<-list()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==LV$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==LV$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(LV$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==LV$Cluster[i]]))
  
  
}

p1<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],ncol = 6)
p2<-ggarrange(g[[7]],g[[8]],g[[9]],g[[10]],g[[11]],g[[12]],ncol = 6)
p3<-ggarrange(g[[13]],g[[14]],g[[15]],g[[16]],g[[17]],g[[18]],ncol = 6)
p<-ggarrange(p1,p2,p3,ncol = 1)
pdf("xLxxxxxxV.pdf",width = 24,height = 12)
p
dev.off()

#--------------xLxxxxxxL----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
LL<-data[which(data$newcluster=="Cluster2"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",LL$Allele,"nmds.txt")
out<-paste0(LL$Allele)
g<-list()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==LL$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==LL$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(LL$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==LL$Cluster[i]]))
  
  
}

empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
    ,plot.margin=unit(c(0.1, 0.1, 0, 0), "inches")
  )
p1<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],ncol = 6)
p2<-ggarrange(g[[7]],g[[8]],g[[9]],g[[10]],g[[11]],g[[12]],ncol = 6)
p3<-ggarrange(g[[13]],g[[14]],g[[15]],g[[16]],g[[17]],g[[18]],ncol = 6)
p4<-ggarrange(g[[19]],g[[20]],g[[21]],g[[22]],g[[23]],empty,ncol = 6)
p<-ggarrange(p1,p2,p3,p4,ncol = 1)
pdf("xLxxxxxxL.pdf",width = 24,height = 12)
p
dev.off()

#--------------xExxxxxxY----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
EY<-data[which(data$motif=="xExxxxxxY"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",EY$Allele,"nmds.txt")
out<-paste0(EY$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==EY$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==EY$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(EY$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==EY$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}
p<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol = 4,nrow = 1)

pdf("xExxxxxxY.pdf",width = 16,height = 4)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("EY.pdf",width = 6,height = 4)
fig1
dev.off()



#--------------xxDxxxxxY----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
DY<-data[which(data$motif=="xxDxxxxxY"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",DY$Allele,"nmds.txt")
out<-paste0(DY$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==DY$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==DY$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(DY$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==DY$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
    ,plot.margin=unit(c(0.1, 0.1, 0, 0), "inches")
  )
p1<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol = 4,nrow = 1)
p2<-ggarrange(g[[5]],g[[6]],g[[7]],empty,ncol = 4,nrow = 1)
p<-ggarrange(p1,p2,ncol = 1,nrow = 2)

pdf("xxDxxxxxY.pdf",width = 16,height = 8)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("DY.pdf",width = 6,height = 4)
fig1
dev.off()

#--------------KxxxxxxxF----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
KF<-data[which(data$motif=="KxxxxxxxF"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",KF$Allele,"nmds.txt")
out<-paste0(KF$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==KF$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==KF$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(KF$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==KF$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}

p<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol = 4,nrow = 1)

pdf("KxxxxxxxF.pdf",width = 16,height = 4)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("KF.pdf",width = 6,height = 4)
fig1
dev.off()
#--------------xxxxRxxxL----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
RL<-data[which(data$motif=="xxxxRxxxL"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",RL$Allele,"nmds.txt")
out<-paste0(RL$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==RL$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==RL$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(RL$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==RL$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}

p<-ggarrange(g[[1]],g[[2]],g[[3]],ncol = 3,nrow = 1)

pdf("xxxxRxxxL.pdf",width = 12,height = 4)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("RL.pdf",width = 6,height = 4)
fig1
dev.off()
#--------------xxxExxxxL----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
EL<-data[which(data$motif=="xxxExxxxL"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",EL$Allele,"nmds.txt")
out<-paste0(EL$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==EL$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==EL$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(EL$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==EL$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}

p<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol = 4,nrow = 1)
pdf("xxxExxxxL.pdf",width = 12,height = 4)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("EL.pdf",width = 6,height = 4)
fig1
dev.off()

#--------------xxxExxxxV----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
EV<-data[which(data$motif=="xxxExxxxV"),]
Fn<-paste0("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/1.submotif/",EV$Allele,"nmds.txt")
out<-paste0(EV$Allele)
g<-list()
pep<-c()
for (i in 1:length(Fn)) {
  nmds<-read.table(Fn[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(nmds$pep,MAP$Epitopes,nomatch = 0)
  nmds$Type<-MAP$Type[index]
  nmds$color<-ifelse(nmds$Cluster==EV$Cluster[i] & nmds$Type=="Canonical","#4DBBD5FF",
                     ifelse(nmds$Cluster==EV$Cluster[i] & nmds$Type=="Non-Canonical","#E64B35FF","#CCCCCC"))
  
  g[[i]]<-ggplot(data = nmds,aes(x = nmds1, y = nmds2,color=color)) + 
    geom_point()+
    scale_color_manual(values = c("#4DBBD5FF" = "#4DBBD5FF","#E64B35FF" = "#E64B35FF","#CCCCCC" = "#CCCCCC"))+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text= element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.title = element_blank())+
    annotate("text", x=min(nmds$nmds1), y=max(nmds$nmds2), 
             label=substr(EV$Allele[i],4,8))+
    annotate("text", x=max(nmds$nmds1), y=min(nmds$nmds2), 
             label=length(nmds$pep[nmds$Cluster==EV$Cluster[i]]))
  pep<-c(pep,nmds$pep[nmds$color!="#CCCCCC"])
}

p<-ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol = 4,nrow = 1)
pdf("xxxExxxxV.pdf",width = 12,height = 4)
p
dev.off()

fig1<-ggseqlogo(pep, seq_type="aa")+
  theme(axis.line = element_line(colour = "black",linewidth = 0.5),
        axis.ticks = element_line(colour = "black",linewidth = 0.5),
        legend.position = "top")
pdf("EV.pdf",width = 6,height = 4)
fig1
dev.off()

#----------------Fig4C----
num<-data.frame(table(data$motif))
num<-num[order(num$Freq,decreasing = T),]
data$motif<-factor(data$motif,levels = as.character(num$Var1))
g4<-ggplot(data = data,aes(x = motif,fill = HLA)) + 
  geom_bar()+
  scale_fill_manual(values = c("HLAA"="#F49B0C","HLAB"="#A3C78E","HLAC"="#F6D780"))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title = element_text(size = 14),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.7),
        legend.position = "top",
        legend.title = element_blank())
pdf("Fig4C.pdf",height = 5,width = 10)
g4
dev.off()

#--------------circ----
library(dendextend)
library(dendsort)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(motifStack)

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
MAP<-MAP[which(nchar(MAP$Epitopes)==9),]
Allele<-unique(MAP$Alleles)
Allele<-Allele[order(Allele)]

Res<-data.frame()
for (i in 1:length(Allele)) {
  c.AA<-unique(MAP$Epitopes[which(MAP$Type=="Canonical" & MAP$Alleles==Allele[i])])
  c.molecularEntropy<-MolecularEntropy(c.AA,type='AA')
  c.PPM<-c.molecularEntropy$freq
  nc.AA<-unique(MAP$Epitopes[which(MAP$Type!="Canonical" & MAP$Alleles==Allele[i])])
  nc.molecularEntropy<-MolecularEntropy(nc.AA,type='AA')
  nc.PPM<-nc.molecularEntropy$freq
  DE<-nc.PPM/c.PPM
  Res<-rbind(Res,t(DE))
  print(i)
}

Res<-as.matrix(Res)
Res<-log2(Res)
Res[is.nan(Res)]<-0
Res<-ifelse(Res=="Inf",max(Res[!is.infinite(Res)]),Res)
Res<-ifelse(Res=="-Inf",min(Res[!is.infinite(Res)]),Res)
rownames(Res)<-1:nrow(Res)

Allele<-gsub("HLA-","",Allele)
split<-rep(Allele,each = 9)
split<-factor(split, levels = Allele)
mycol<-colorRamp2(c(-4, 0,4),c("blue","white","red"))

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
pdf("FigS3D.pdf",height = 10,width = 10)
circos.heatmap(Res,col=mycol,split=split,track.height = 0.6,
               bg.border = T,rownames.side = "outside",cluster = FALSE)
lg=Legend(title="log2(Ratio)",col_fun=mycol,direction = c("vertical"))
grid.draw(lg)
#添加列名：
circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  if(CELL_META$sector.numeric.index==1){   #if(CELL_META$sector.numeric.index == 3) { # the last sector
    cn=colnames(Res)
    n=length(cn)
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.5,"mm"),#x坐标
                1:n,#调整y坐标
                cn,cex=0.6,adj=c(0,0.5),facing="inside")}
},bg.border=NA)
circos.clear()
dev.off()

rm(list = ls())
#--------------------motif
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
MAP<-MAP[which(nchar(MAP$Epitopes)==9),]
Allele<-unique(MAP$Alleles)
Allele<-Allele[order(Allele)]

fig<-list()
for (i in 1:length(Allele)) {
  HLA<-Allele[i]
  pep1<-unique(MAP$Epitopes[MAP$Alleles==HLA & MAP$Type=="Canonical"])
  pep2<-unique(MAP$Epitopes[MAP$Alleles==HLA & MAP$Type=="Non-Canonical"])
  fig1<-ggseqlogo(pep1, seq_type="aa")+
    theme(axis.line = element_line(colour = "black",linewidth = 0.5),
          axis.ticks = element_line(colour = "black",linewidth = 0.5),
          legend.position = "none")
  fig2<-ggseqlogo(pep2, seq_type="aa")+
    theme(axis.line = element_line(colour = "black",linewidth = 0.5),
          axis.ticks = element_line(colour = "black",linewidth = 0.5),
          legend.position = "none")
  
  fig[[i]]<-ggarrange(fig1,fig2,ncol = 2,labels = paste0(gsub("HLA-","",HLA),
                                                         c(": cMAPs",": ncMAPs")))
  print(i)
}
p1<-ggarrange(fig[[1]],fig[[2]],fig[[3]],fig[[4]],fig[[5]],ncol = 5,nrow = 1)
p2<-ggarrange(fig[[6]],fig[[7]],fig[[8]],fig[[9]],fig[[10]],ncol = 5,nrow = 1)
p3<-ggarrange(fig[[11]],fig[[12]],fig[[13]],fig[[14]],fig[[15]],ncol = 5,nrow = 1)
p4<-ggarrange(fig[[16]],fig[[17]],fig[[18]],fig[[19]],fig[[20]],ncol = 5,nrow = 1)
p5<-ggarrange(fig[[21]],fig[[22]],fig[[23]],fig[[24]],fig[[25]],ncol = 5,nrow = 1)
p6<-ggarrange(fig[[26]],fig[[27]],fig[[28]],fig[[29]],fig[[30]],ncol = 5,nrow = 1)
p7<-ggarrange(fig[[31]],fig[[32]],fig[[33]],fig[[34]],fig[[35]],ncol = 5,nrow = 1)
p8<-ggarrange(fig[[36]],fig[[37]],fig[[38]],fig[[39]],fig[[40]],ncol = 5,nrow = 1)
p9<-ggarrange(fig[[41]],fig[[42]],fig[[43]],fig[[44]],fig[[45]],ncol = 5,nrow = 1)
pp<-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol = 1,nrow = 9)


setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/4ncMAP呈递特征/2.motif/fig")
pdf("FigS3.pdf",height = 36,width = 60)
pp
dev.off()

tiff("FigS3.tif",height = 36,width = 60,res = 300,units = "in")
pp
dev.off()


