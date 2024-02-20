library(ggplot2)
library(ggpubr)
#==========Figure 1C==========#
setwd("~/Desktop/22白血病/Runtime/0rawdata/PXD018542/")
PEP<-read.table("peptides.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
PEP<-PEP[,c("Sequence","PEP")]

setwd("/Users/caiyangyang/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

index<-match(data$Epitopes,PEP$Sequence,nomatch = 0)
data$PEP<-PEP$PEP[index]
data<-data[,c("Epitopes","PEP","Type")]
data<-data[!duplicated(data),]

g1<-ggviolin(data, x="Type", y="PEP",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
              fill = "Type",short.panel.labs = F,outlier.shape = NA,
             add = "boxplot", add.params = list(fill="white"),
              ylab = "Posterior Error Probability(PEP)")+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title.y = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  stat_compare_means(aes(group = Type),method="wilcox.test",
                     label = "p.signif",tip.length=0,label.x = 1.5)

setwd("~/Desktop/22白血病/Runtime/2特征比较/")
pdf("1C.pdf",width = 6,height = 5)
g1
dev.off()
#==============================#


#==========Figure 1D==========#
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data<-data[,c("Epitopes","Type")]
data<-data[!duplicated(data),]

hydrophobicity<-c(1.8,2.5,-3.5,-3.5,2.8,-0.4,
                  -3.2,4.5,-3.9,3.8,1.9,-3.5,
                  -1.6,-3.5,-4.5,-0.8,-0.7,4.2,
                  -0.9,-1.3)
names(hydrophobicity)<-c("A","C","D","E","F","G",
                         "H","I","K","L","M","N",
                         "P","Q","R","S","T","V",
                         "W","Y")
data$hydrophobicity_score<-unlist(lapply(data$Epitopes,function(x){
  pep<-unlist(strsplit(x,""))
  index<-match(pep,names(hydrophobicity),nomatch = 0)
  h<-hydrophobicity[index]
  n<-length(pep)
  m<-4:(n-1)
  return(sum(h[m]))}))


g2<-ggviolin(data, x="Type", y="hydrophobicity_score",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
             fill = "Type",short.panel.labs = F,add = "boxplot", add.params = list(fill="white"),
             outlier.shape = NA,ylab = "T-cell Contact Hydrophobicity")+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title.y = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  stat_compare_means(aes(group = Type),method="wilcox.test",
                     label = "p.signif",tip.length=0,label.x = 1.5)

setwd("~/Desktop/22白血病/Runtime/2特征比较/")
pdf("1D.pdf",width = 6,height = 5)
g1
dev.off()
#==============================#

#==========Figure 1E==========#
setwd("/Users/caiyangyang/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data1<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

setwd("~/Desktop/22白血病/Runtime/0rawdata/DeepLC")
data<-read.csv("Non-Canonical+Canonical(bind+nonbind)_pre.csv")
data<-data[data$seq %in% data1$Epitopes,]

index<-match(data$seq,data1$Epitopes,nomatch = 0)
data$Type<-data1$Type[index]
colnames(data)[4:5]<-c("Observed Retention Time (min)","DeepLC Retention Time (min)")

g1<-ggscatterhist(data, x = "Observed Retention Time (min)", y = "DeepLC Retention Time (min)",
                  color = "Type", size = 1,
                  palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
                  add = "reg.line", alpha = 0.6,group = "Type",
                  margin.params = list(fill ="Type"),
                  add.params = list(color = "Type"), print = F,
                  cor.coef = TRUE, # 添加相关系数
                  cor.coeff.args = list(method = "pearson",label.x=3,label.y=115))

setwd("~/Desktop/22白血病/Runtime/2特征比较/")
pdf("1E.pdf",width = 6,height = 5,onefile = F)
g1
dev.off()
#==============================#


#==========Figure 1F==========#
setwd("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/0-fasta")
canonical<-read.table("uniprot.fa",header = F,sep = "\t",check.names = F,stringsAsFactors = F,fill = T,quote = "")
n<-nrow(canonical)
len<-nchar(canonical[,1])[seq(2,n,2)]
canonical<-data.frame(ORF = canonical[seq(1,n,2),1],len)
canonical$ORF<-gsub(">","",canonical$ORF)

setwd("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/1-SeqName/")
non_canonical<-read.table("5nuORFRename.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,fill = T,quote = "")
non_canonical<-data.frame(ORF=non_canonical$ORF_ID_Rename,len=nchar(non_canonical$AA_seq))

ORF<-rbind(canonical,non_canonical)

setwd("~/Desktop/22白血病/Runtime/1预处理/")
data<-read.table("All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

setwd("/Users/caiyangyang/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

index<-match(MAP$Epitopes,data$Epitopes,nomatch = 0)
MAP$ORF<-data$`Leading razor protein`[index]

index<-match(MAP$ORF,ORF$ORF,nomatch = 0)
MAP$source_protein_len<-ORF$len[index]
source<-MAP[,c("Epitopes","Type","source_protein_len")]
source<-source[!duplicated(source),]

g4<-gghistogram(source,x="source_protein_len",y="..density..",fill = "Type",rug = T,
                color = "Type",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
                ylab = "Density",xlab = "Source ORF Length(aa)")+
  scale_x_continuous(limits = c(0,1000),breaks = seq(0,1000,100))

setwd("~/Desktop/22白血病/Runtime/2特征比较")
pdf("1F.pdf",width = 6,height = 5)
g4
dev.off()
#==============================#


#==========Figure 1G==========#
setwd("/Users/caiyangyang/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")

data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data<-data[,c("Samples","Epitopes","Type")]
data<-data[!duplicated(data),]

s_e_f<-data.frame(table(data[,c("Epitopes","Type")]))
s_e_f<-s_e_f[which(s_e_f$Freq!=0),]
ff<-data.frame(table(s_e_f[,2:3]))
ff<-ff[which(ff$Freq.1!=0),]
sum_c<-sum(ff$Freq.1[which(ff$Type=="Canonical")])
sum_n<-sum(ff$Freq.1[which(ff$Type!="Canonical")])

ff$Freq.1<-ifelse(ff$Type=="Canonical",ff$Freq.1/sum_c*100,ff$Freq.1/sum_n*100)

g1<-ggplot(ff, aes(x = Freq, y = Freq.1,color= Type,shape = Type,group = Type))+
  geom_point(size=2)+
  geom_line(size=1)+
  theme(panel.background = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  labs(x = "Occure in Samples",y="% of MAPs")+
  scale_color_manual(values=c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"))

setwd("~/Desktop/22白血病/Runtime/2特征比较/")
pdf("1G.pdf",width = 6,height = 5)
g1
dev.off()
#==============================#
library(ggplot2)
library(ggpubr)
library(VennDiagram)
#==========Figure S1A==========#
setwd("~/Desktop/22白血病/Runtime/1预处理/")
data<-read.table("All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
sample<-unique(data$Samples)

AML_peptide<-data[,c("Samples","Epitopes","Type")]
AML_peptide<-AML_peptide[!duplicated(AML_peptide),]
AML_peptide_num<-data.frame(table(AML_peptide[,c("Samples","Type")]))

NetMHCpan<-data[which(data$BindLevel=="Binding"),c("Samples","Epitopes","Type")]
NetMHCpan<-NetMHCpan[!duplicated(NetMHCpan),]

MHCflurry<-data[which(data$MHCflurryBindLevel=="Binding"),c("Samples","Epitopes","Type")]
MHCflurry<-MHCflurry[!duplicated(MHCflurry),]

MixMHCpred<-data[which(data$MixMHCpredBindLevel=="Binding"),c("Samples","Epitopes","Type")]
MixMHCpred<-MixMHCpred[!duplicated(MixMHCpred),]

MAP_freq<-data.frame()
for (i in 1:19) {
  A<-NetMHCpan[which(NetMHCpan$Samples==sample[i]),]
  B<-MHCflurry[which(MHCflurry$Samples==sample[i]),]
  C<-MixMHCpred[which(MixMHCpred$Samples==sample[i]),]
  Epitopes<-rbind(A,B,C)
  Epitopes_Freq<-data.frame(table(Epitopes))
  Epitopes_Freq<-Epitopes_Freq[which(Epitopes_Freq$Freq!=0),]
  MAP_freq<-rbind(MAP_freq,Epitopes_Freq)
}

colnames(MAP_freq)[4]<-"Predict"
AML_support_num<-data.frame(table(MAP_freq[,c("Samples","Type","Predict")]))
AML_support_num$Predict<-ifelse(AML_support_num$Predict==1,"1 Predictor",
                                paste0(AML_support_num$Predict," Predictors"))

AML_MAP_num<-data.frame(table(MAP_freq[,c("Samples","Type")]))
aa<-paste0(AML_peptide_num$Samples,"_",AML_peptide_num$Type)
bb<-paste0(AML_MAP_num$Samples,"_",AML_MAP_num$Type)
index<-match(aa,bb,nomatch = 0)
AML_MAP_num<-AML_MAP_num[index,]

Non_MAP_num<-AML_peptide_num$Freq-AML_MAP_num$Freq
Non_MAP_num<-cbind(AML_MAP_num[,1:2],Predict = "Non-Binding",Freq = Non_MAP_num)

ALL<-rbind(AML_support_num,Non_MAP_num)
ALL$fill<-paste0(ALL$Type,"_",ALL$Predict)
ALL$fill<-factor(ALL$fill,levels = c("Canonical_Non-Binding","Canonical_1 Predictor",
                                     "Canonical_2 Predictors","Canonical_3 Predictors",
                                     "Non-Canonical_Non-Binding","Non-Canonical_1 Predictor",
                                     "Non-Canonical_2 Predictors","Non-Canonical_3 Predictors"))
ALL$Freq<-ifelse(ALL$Type=="Canonical",-ALL$Freq,ALL$Freq)


Canonical_num<-AML_peptide_num$Freq[which(AML_peptide_num$Type=="Canonical")]
NonCanonical_num<-AML_peptide_num$Freq[which(AML_peptide_num$Type!="Canonical")]
Peptide_num<-Canonical_num+NonCanonical_num
sample<-sample[order(Peptide_num,decreasing = T)]
ALL$Samples<-factor(ALL$Samples,levels = sample)

g1<-ggplot(ALL,aes(x = Samples,y = Freq,fill = fill))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = c("Canonical_1 Predictor" = "#C3E8F1",
                               "Canonical_2 Predictors" = "#88D1E3",
                               "Canonical_3 Predictors" = "#4DBBD5FF",
                               "Canonical_Non-Binding"= "#CCCCCC",
                               "Non-Canonical_1 Predictor" = "#F6C3BB",
                               "Non-Canonical_2 Predictors" = "#EE8778",
                               "Non-Canonical_3 Predictors" = "#E64B35FF",
                               "Non-Canonical_Non-Binding"= "#999999"))+
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        axis.ticks.length=unit(0.1,'cm'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 14,colour = "black"),
        legend.text = element_text(size = 12,colour = "black"),
        legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.background=element_rect(fill="white",colour="black",size=0.5))+
  ylab("Number of Peptides Identified by Immunopeptidome")+
  coord_flip()+
  facet_grid(~Type,scales = "free")

ALL$Samples<-factor(ALL$Samples,levels = paste0("AML",19:1))
g2<-ggplot(ALL,aes(x = Samples,y = Freq,fill = fill))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = c("Canonical_1 Predictor" = "#C3E8F1",
                               "Canonical_2 Predictors" = "#88D1E3",
                               "Canonical_3 Predictors" = "#4DBBD5FF",
                               "Canonical_Non-Binding"= "#CCCCCC",
                               "Non-Canonical_1 Predictor" = "#F6C3BB",
                               "Non-Canonical_2 Predictors" = "#EE8778",
                               "Non-Canonical_3 Predictors" = "#E64B35FF",
                               "Non-Canonical_Non-Binding"= "#999999"))+
  theme(axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        axis.ticks.length=unit(0.1,'cm'),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 14,colour = "black"),
        legend.text = element_text(size = 12,colour = "black"),
        legend.position = c(0.8,0.8),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.background=element_rect(fill="white",colour="black",size=0.5))+
  ylab("Number of Peptides Identified by Immunopeptidome")+
  coord_flip()+
  facet_grid(~Type,scales = "free")



setwd("~/Desktop/22白血病/Runtime/2特征比较")
write.table(MAP_freq,"MAP_support.txt",quote = F,sep = "\t",row.names = F)

pdf("S1A.pdf",width = 12,height = 5)
g1
dev.off()

pdf("比例.pdf",width = 12,height = 5,onefile = F)
g2
dev.off()


#==============================#
venn<-venn.diagram(x=list("NetMHCpan" = unique(NetMHCpan$Epitopes[which(NetMHCpan$Type=="Canonical")]),
                          "MHCflurry" = unique(MHCflurry$Epitopes[which(MHCflurry$Type=="Canonical")]),
                          "MixMHCpred" = unique(MixMHCpred$Epitopes[which(MixMHCpred$Type=="Canonical")])),
                   filename = NULL,scaled = F,
                   fill = c( "#EF7E0A","#2C7BB9","#3CAA39"))
pdf("Canonical_Venn.pdf")
grid.draw(venn)
dev.off()

venn<-venn.diagram(x=list("NetMHCpan" = unique(NetMHCpan$Epitopes[which(NetMHCpan$Type!="Canonical")]),
                          "MHCflurry" = unique(MHCflurry$Epitopes[which(MHCflurry$Type!="Canonical")]),
                          "MixMHCpred" = unique(MixMHCpred$Epitopes[which(MixMHCpred$Type!="Canonical")])),
                   filename = NULL,scaled = F,
                   fill = c( "#EF7E0A","#2C7BB9","#3CAA39"))
pdf("NonCanonical_Venn.pdf")
grid.draw(venn)
dev.off()

#==========Figure S1C==========#
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data<-data[,c("Epitopes","Type")]
data<-data[!duplicated(data),]

hydrophobicity_AA<-c("V","I","L","F","M","W","P","G","A")
data$Hydrophobicity_Fraction<-unlist(lapply(data$Epitopes,function(x){
  pep<-unlist(strsplit(x,""))
  index<-match(pep,hydrophobicity_AA,nomatch = 0)
  n<-length(which(pep %in% hydrophobicity_AA))
  m<-length(pep)
  return(n/m)}))

g2<-ggviolin(data, x="Type", y="Hydrophobicity_Fraction",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
             fill = "Type",short.panel.labs = F,add = "boxplot", add.params = list(fill="white"),
             outlier.shape = NA,ylab = "Hydrophobicity Fraction")+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title.y = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  stat_compare_means(aes(group = Type),method="wilcox.test",
                     label = "p.signif",tip.length=0,label.x = 1.5)

setwd("~/Desktop/22白血病/Runtime/2特征比较")
pdf("S1C.pdf",width = 6,height = 5)
g2
dev.off()

#==============================#

#==========Figure S1D==========#
setwd("/Users/caiyangyang/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data$Length<-nchar(data$Epitopes)

data<-data[,c("Epitopes","Length","Type")]
data<-data[!duplicated(data),]
data<-data.frame(table(data[,c("Length","Type")]))
sum_c<-sum(data$Freq[which(data$Type=="Canonical")])
sum_n<-sum(data$Freq[which(data$Type!="Canonical")])
data$percentage<-ifelse(data$Type=="Canonical",round(data$Freq/sum_c*100,1),round(data$Freq/sum_n*100,1))

g3<-ggplot(data, aes(x = Length, y = percentage,fill = Type))+
  geom_bar(stat="identity",position="dodge")+
  geom_text(aes(label=Freq),size=4,position = position_dodge(width = 0.9))+
  theme(panel.background = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  labs(x = "Length of MAPs",y="% of MAPs")+
  scale_fill_manual(values=c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"))

setwd("~/Desktop/22白血病/Runtime/2特征比较")
pdf("S1D.pdf",width = 6,height = 5)
g3
dev.off()
#==============================#

#==========Figure 2==========#
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
pdf("2B.pdf",width = 6,height = 5,onefile = F)
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

pdf("S2C.pdf",width = 6,height = 5)
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



