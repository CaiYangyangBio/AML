library(ggplot2)
library(ggpubr)
library(VennDiagram)
#--------------------------------------1.样本分组条形图----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/1预处理/")
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




setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/1.样本分组条形图")
write.table(MAP_freq,"MAP_support.txt",quote = F,sep = "\t",row.names = F)

pdf("S1A.pdf",width = 12,height = 5)
g1
dev.off()

#------------------------Venn
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

#--------------------------------------2.至少被2个软件支持的肽----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/1预处理/")
data<-read.table("All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data<-data[,c("Samples","Epitopes","Alleles","Type","BindLevel","MHCflurryBindLevel","MixMHCpredBindLevel")]
data$Bindnum<-apply(data,1,function(x){length(which(x=="Binding"))})
data<-data[which(data$Bindnum>1),]

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图")
write.table(data,"MAP_Allele_support.txt",quote = F,sep = "\t",row.names = F)


#--------------------------------------3.PEP----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/0rawdata/PXD018542/")
PEP<-read.table("peptides.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
PEP<-PEP[,c("Sequence","PEP")]

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

index<-match(data$Epitopes,PEP$Sequence,nomatch = 0)
data$PEP<-PEP$PEP[index]
data<-data[,c("Epitopes","PEP","Type")]
data<-data[!duplicated(data),]

g1<-ggboxplot(data, x="Type", y="PEP",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
              fill = "Type",short.panel.labs = F,outlier.shape = NA,
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

pdf("1C.pdf",width = 6,height = 5)
g1
dev.off()

#--------------------------------------4.length----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data$Length<-nchar(data$Epitopes)

data<-data[,c("Epitopes","Length","Type")]
data<-data[!duplicated(data),]
data<-data.frame(table(data[,c("Length","Type")]))
sum_c<-sum(data$Freq[which(data$Type=="Canonical")])
sum_n<-sum(data$Freq[which(data$Type!="Canonical")])
data$percentage<-ifelse(data$Type=="Canonical",round(data$Freq/sum_c*100,1),round(data$Freq/sum_n*100,1))

g1<-ggplot(data, aes(x = Length, y = percentage,fill = Type))+
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

pdf("S1H.pdf",width = 6,height = 5)
g1
dev.off()

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/0rawdata/1-ncORF/0-fasta")
canonical<-read.table("uniprot.fa",header = F,sep = "\t",check.names = F,stringsAsFactors = F,fill = T,quote = "")
n<-nrow(canonical)
len<-nchar(canonical[,1])[seq(2,n,2)]
canonical<-data.frame(ORF = canonical[seq(1,n,2),1],len)
canonical$ORF<-gsub(">","",canonical$ORF)

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/0rawdata/1-ncORF/1-SeqName/")
non_canonical<-read.table("5nuORFRename.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,fill = T,quote = "")
non_canonical<-data.frame(ORF=non_canonical$ORF_ID_Rename,len=nchar(non_canonical$AA_seq))

ORF<-rbind(canonical,non_canonical)

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/1预处理/")
data<-read.table("All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

index<-match(MAP$Epitopes,data$Epitopes,nomatch = 0)
MAP$ORF<-data$`Leading razor protein`[index]

index<-match(MAP$ORF,ORF$ORF,nomatch = 0)
MAP$source_protein_len<-ORF$len[index]
source<-MAP[,c("Epitopes","Type","source_protein_len")]
source<-source[!duplicated(source),]

g2<-gghistogram(source,x="source_protein_len",y="..density..",fill = "Type",rug = T,
                color = "Type",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
                ylab = "Density",xlab = "Source ORF Length(aa)")+
  scale_x_continuous(limits = c(0,1000),breaks = seq(0,1000,100))

pdf("1F.pdf",width = 6,height = 5)
g2
dev.off()

#--------------------------------------5.hydrophobicity----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
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

g1<-ggboxplot(data, x="Type", y="Hydrophobicity_Fraction",palette = c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"),
              fill = "Type",short.panel.labs = F,
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

pdf("S1D.pdf",width = 6,height = 5)
g1
dev.off()

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
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

pdf("1D.pdf",width = 6,height = 5)
g2
dev.off()


#--------------------------------------5.RT----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
data1<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/0rawdata/DeepLC")
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

pdf("1E.pdf",width = 6,height = 5,onefile = F)
g1
dev.off()


#--------------------------------------6.Freq----
setwd("/Volumes/Seagate/work/2022-2-白血病/Runtime/2特征比较/2.等位分组条形图/")
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
  labs(x = "Number of Samples",y="% of MAPs")+
  scale_color_manual(values=c("Canonical"="#4DBBD5FF","Non-Canonical"="#E64B35FF"))


pdf("1G.pdf",width = 6,height = 5)
g1
dev.off()
