library(ggplot2)
library(ggpubr)
library(VennDiagram)
#==========1.样本分组条形图==========#
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



setwd("~/Desktop/22白血病/Runtime/2特征比较/1.样本分组条形图")
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

