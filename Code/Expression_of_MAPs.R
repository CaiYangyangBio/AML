library(ggplot2)
library(ggpubr)
library(ggrepel)

setwd("~/Desktop/22白血病/Runtime/0rawdata/")
LAML<-read.table("LAML.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
setwd("~/Desktop/22白血病/Runtime/2特征比较/10表达/")
GRCh38<-read.table("gene_type_anno.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
coding_gene<-unique(GRCh38$gene_name[which(GRCh38$gene_type=="protein_coding")])
coding<-intersect(coding_gene,rownames(LAML))

noncoding_gene<-unique(GRCh38$gene_name[which(GRCh38$gene_type!="protein_coding")])
noncoding<-intersect(noncoding_gene,rownames(LAML))

#======#======#======背景的基因表达======#======#======#
#======#======经典的======#======#
setwd("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/1-SeqName")
cORF<-read.table("uniprot.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
cORF$Name<-gsub(">","",cORF$Name)
cORF$RawName<-unlist(lapply(strsplit(cORF$RawName,"GN="),function(x){x[2]}))
cORF<-cORF[!is.na(cORF$RawName),]
cORF$RawName<-unlist(lapply(strsplit(cORF$RawName," "),function(x){x[1]}))
cGene<-unique(cORF$RawName)
cGene<-data.frame(Gene = cGene,Type = "Protein-coding")

#======#======非经典的======#======#
setwd("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/1-SeqName/")
ncORF<-read.table("5nuORFRename.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
ncGene<-unique(ncORF$geneName)
ncGene<-data.frame(Gene = ncGene,Type = "Non-coding")

#======#======Gene======#======#
Gene<-rbind(cGene,ncGene)
Gene$Type<-ifelse(Gene$Gene%in%coding_gene,"Protein-coding",Gene$Type)
Gene<-Gene[!duplicated(Gene),]
int<-intersect(rownames(LAML),Gene$Gene)
index<-match(int,Gene$Gene,nomatch = 0)
Gene<-Gene[index,]
LAML<-LAML[int,]
Gene$Exp<-apply(LAML, 1, function(x){quantile(x,0.9)})
Gene<-Gene[order(Gene$Exp),]
Gene$GeneRank<-1:nrow(Gene)
Gene$Freq<-0

#======#======#======MAP的基因表达======#======#======#
setwd("~/Desktop/22白血病/Runtime/1预处理/")
Pep<-read.table("All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
MAP<-merge(MAP,Pep)

#======#======Canonical======#======#
Canonical<-MAP[which(MAP$Type=="Canonical"),]
index<-match(Canonical$`Leading razor protein`,cORF$Name,nomatch = 0)
Canonical<-Canonical[which(index!=0),]
Canonical$Symbol<-cORF$RawName[index]
Canonical<-Canonical[,c("Epitopes","Symbol")]
Canonical<-Canonical[!duplicated(Canonical),]
Canonical<-data.frame(table(Canonical$Symbol))

#======#======Non-Canonical======#======#
NonCanonical<-MAP[which(MAP$Type!="Canonical"),]
NonCanonical$Symbol<-unlist(lapply(strsplit(NonCanonical$`Leading razor protein`,"_"),function(x){x[1]}))
NonCanonical<-NonCanonical[,c("Epitopes","Symbol")]
NonCanonical<-NonCanonical[!duplicated(NonCanonical),]
NonCanonical<-data.frame(table(NonCanonical$Symbol))


#======#======Share======#======#
Share<-intersect(Canonical$Var1,NonCanonical$Var1)
index1<-match(Share,Canonical$Var1,nomatch = 0)
index2<-match(Share,NonCanonical$Var1,nomatch = 0)
Share<-data.frame(Gene = Share,
                  Freq = apply(cbind(Canonical$Freq[index1],NonCanonical$Freq[index2]),1,sum))

Canonical<-Canonical[!Canonical$Var1%in%Share$Gene,]
NonCanonical<-NonCanonical[!NonCanonical$Var1%in%Share$Gene,]

#======#======#======合并数据======#======#======#
index3<-match(Canonical$Var1,Gene$Gene,nomatch = 0)
Gene$Freq[index3]<-Canonical$Freq[which(index3!=0)]
Gene$Type[index3]<-"Protein-coding derived cMAP"

index4<-match(Share$Gene,Gene$Gene,nomatch = 0)
Gene$Freq[index4]<-Share$Freq[which(index4!=0)]
Gene$Type[index4]<-"Protein-coding derived MAP"

index5<-match(NonCanonical$Var1,Gene$Gene,nomatch = 0)
Gene$Freq[index5]<-NonCanonical$Freq[which(index5!=0)]
Gene$Type[index5]<-ifelse(Gene$Type[index5]=="Protein-coding","Protein-coding derived ncMAP","Non-coding derived ncMAP")

Gene$size<-ifelse(Gene$Freq>10,">10 peptides",
                  ifelse(Gene$Freq>5 & Gene$Freq<=10,"6-10 peptides",
                         ifelse(Gene$Freq>1 & Gene$Freq<=5,"2-5 peptides",
                                ifelse(Gene$Freq==1,"1 peptide","No peptide"))))
Gene$Type<-factor(Gene$Type,levels = c("Protein-coding","Non-coding",
                                       "Protein-coding derived cMAP",
                                       "Protein-coding derived MAP",
                                       "Protein-coding derived ncMAP",
                                       "Non-coding derived ncMAP"))

label<-Gene[which(Gene$size==">10 peptides"),]
g1<-ggscatter(Gene, x = "GeneRank", y = "Exp",color = "Type",size = "size",
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36",
                          "Protein-coding"="#ACADAD","Non-coding"="#D8D7D7"),
              xlab = "Gene Rank",ylab = "log2(TPM+1)",alpha = 1)+
  scale_y_continuous(limits = c(0,16))+
  theme(legend.position = "none")+
  scale_size_manual(values = c("No peptide"=0.5, "1 peptide"=1,"2-5 peptides"=1.5,
                               "6-10 peptides"=2,">10 peptides"=2.5))+
  geom_text_repel(data=label,aes(x=GeneRank,y=Exp,colour=Type,label=Gene),
                  size = 4,max.overlaps = 30)

g2<-ggdensity(Gene,x = "Exp",y = "..count..",fill = "Type",
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36",
                          "Protein-coding"="#ACADAD","Non-coding"="#D8D7D7"),
              xlab = "",ylab = "Gene frequency",alpha = 0.9)+
  scale_x_continuous(limits = c(0,16))+
  theme(legend.position = "none")+
  coord_flip()

MAPGene<-Gene[which(Gene$Type=="Protein-coding derived cMAP"|
                      Gene$Type=="Protein-coding derived ncMAP"|
                      Gene$Type=="Protein-coding derived MAP"|
                      Gene$Type=="Non-coding derived ncMAP"),]

MAPGene$Type<-factor(MAPGene$Type,levels = c("Protein-coding derived cMAP",
                                             "Protein-coding derived ncMAP",
                                             "Protein-coding derived MAP",
                                             "Non-coding derived ncMAP"))
# g3<-ggdensity(MAPGene,x = "Exp",y = "..count..",fill = "Type",
#               palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
#                           "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36"),
#               xlab = "",ylab = "Gene frequency",alpha = 0.9)+
#   theme(legend.position = "none")+
#   scale_x_continuous(limits = c(0,16))+
#   scale_y_continuous(limits = c(0,50))+
#   coord_flip()
Gene$Type<-factor(Gene$Type,levels = c("Protein-coding","Non-coding",
                                       "Protein-coding derived cMAP",
                                       "Protein-coding derived ncMAP",
                                       "Protein-coding derived MAP",
                                       "Non-coding derived ncMAP"))
g3<-ggdensity(Gene,x = "Exp",y = "..count..",fill = "Type",
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36",
                          "Protein-coding"="#ACADAD","Non-coding"="#D8D7D7"),
              xlab = "",ylab = "Gene frequency",alpha = 0.9)+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,16))+
  scale_y_continuous(limits = c(0,50))+
  coord_flip()
p1<-ggarrange(g1,g2,g3,ncol = 3,widths = c(1,1,1))

#======#======#======sharegene======#======#======#
#======#======Canonical======#======#
Canonical<-MAP[which(MAP$Type=="Canonical"),]
index<-match(Canonical$`Leading razor protein`,cORF$Name,nomatch = 0)
Canonical<-Canonical[which(index!=0),]
Canonical$Symbol<-cORF$RawName[index]
Canonical<-Canonical[,c("Epitopes","Symbol")]
Canonical<-Canonical[!duplicated(Canonical),]

#======#======Non-Canonical======#======#
NonCanonical<-MAP[which(MAP$Type!="Canonical"),]
NonCanonical$Symbol<-unlist(lapply(strsplit(NonCanonical$`Leading razor protein`,"_"),function(x){x[1]}))
NonCanonical<-NonCanonical[,c("Epitopes","Symbol")]
NonCanonical<-NonCanonical[!duplicated(NonCanonical),]

Share<-intersect(Canonical$Symbol,NonCanonical$Symbol)
data<-rbind(Canonical,NonCanonical)
data<-data[data$Symbol%in%Share,]
data<-merge(data,MAP)
data<-data[,c("Epitopes","Symbol","PlotType")]
data<-data[!duplicated(data),]
data1<-data.frame(table(data[,2:3]))
data1<-data1[which(data1$Freq!=0),]
data1$PlotType<-factor(data1$PlotType,
                       levels = rev(c("Canonical","lncRNA","3' Overlap dORF",
                                      "5' uORF","3' dORF","5' Overlap uORF",
                                      "Pseudogene","Out-of-Frame","Other")))
data<-data.frame(table(data$Symbol))
data<-data[order(data$Freq),]
data1$Symbol<-factor(data1$Symbol,levels = as.character(data$Var1))


color<-c("Canonical"="#4DBBD5FF","lncRNA"="#81BC7C","3' Overlap dORF"="#FFB9B4",
         "5' uORF"="#D8E289","3' dORF"="#FCAE3F","5' Overlap uORF"="#C3B5E6",
         "Pseudogene"="#FADD5B","Out-of-Frame"="#FFDDB8","Other"="#7466A3")

p2<-ggplot(data1, aes(x = Symbol, y = Freq,fill = PlotType))+
  geom_bar(stat="identity",position="fill")+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size =12,color="black",angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size =12,color="black"),
        axis.title = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "right",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  labs(x = "",y="Fraction of ORF")+
  scale_fill_manual(values=color)


#======#======#======表达比较======#======#======#
MAPGene$Type<-factor(MAPGene$Type,levels = c("Non-coding derived ncMAP",
                                             "Protein-coding derived ncMAP",
                                             "Protein-coding derived cMAP",
                                             "Protein-coding derived MAP"))
mycompare<-list(c("Non-coding derived ncMAP","Protein-coding derived ncMAP"),
                c("Non-coding derived ncMAP","Protein-coding derived cMAP"),
                c("Non-coding derived ncMAP","Protein-coding derived MAP"),
                c("Protein-coding derived ncMAP","Protein-coding derived cMAP"),
                c("Protein-coding derived ncMAP","Protein-coding derived MAP"),
                c("Protein-coding derived cMAP","Protein-coding derived MAP"))

p3<-ggboxplot(MAPGene,x = "Type",y = "Exp",fill = "Type",xlab = "",ylab = "log2(TPM+1)",outlier.shape = NA,
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36"))+
  stat_compare_means(comparisons = mycompare,method="wilcox.test",label = "p.signif",tip.length=0,label.x = 1.5)

#======#======#======原蛋白比较======#======#======#
#======#======Canonical======#======#
Canonical<-MAP[which(MAP$Type=="Canonical"),]
index<-match(Canonical$`Leading razor protein`,cORF$Name,nomatch = 0)
Canonical<-Canonical[which(index!=0),]
Canonical$Symbol<-cORF$RawName[index]
Canonical<-Canonical[,c("Epitopes","Symbol","Leading razor protein","Type","PlotType")]
Canonical<-Canonical[!duplicated(Canonical),]

#======#======Non-Canonical======#======#
NonCanonical<-MAP[which(MAP$Type!="Canonical"),]
NonCanonical$Symbol<-unlist(lapply(strsplit(NonCanonical$`Leading razor protein`,"_"),function(x){x[1]}))
NonCanonical<-NonCanonical[,c("Epitopes","Symbol","Leading razor protein","Type","PlotType")]
NonCanonical<-NonCanonical[!duplicated(NonCanonical),]

ALL<-rbind(Canonical,NonCanonical)

Gene<-rbind(cGene,ncGene)
Gene$Type<-ifelse(Gene$Gene%in%coding_gene,"Protein-coding",Gene$Type)
Gene<-Gene[!duplicated(Gene),]
ALL$GeneType<-""
index<-match(ALL$Symbol,Gene$Gene,nomatch = 0)
ALL$GeneType[which(index!=0)]<-Gene$Type[index]
ALL$GeneType[which(index==0)]<-"Non-coding"

ALL$MAPType<-""
ALL$MAPType<-ifelse(ALL$GeneType=="Protein-coding" & ALL$Type=="Canonical","Protein-coding derived cMAP",
                    ifelse(ALL$GeneType=="Protein-coding" & ALL$Type!="Canonical","Protein-coding derived ncMAP","Non-coding derived ncMAP"))
ALL$MAPType[ALL$Symbol%in%Share]<-"Protein-coding derived MAP"

#======#======Non-Canonical======#======#
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
source<-merge(source,ALL)
source$MAPType<-factor(source$MAPType,levels = c("Non-coding derived ncMAP",
                                                 "Protein-coding derived ncMAP",
                                                 "Protein-coding derived cMAP",
                                                 "Protein-coding derived MAP"))
p4<-ggboxplot(source, x="MAPType", y="source_protein_len",fill = "MAPType",xlab = "",ylab = "Source ORFs Length",outlier.shape = NA,
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36"))+
  stat_compare_means(comparisons = mycompare,method="wilcox.test",label = "p.signif",tip.length=0,label.x = 1.5)+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title.y = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())+
  scale_y_log10()

source$len<-nchar(source$Epitopes)/source$source_protein_len
p5<-ggboxplot(source, x="MAPType", y="len",fill = "MAPType",xlab = "",
              ylab = expression(frac("MAPs Covered Length","Source ORFs Length")),outlier.shape = NA,
              palette = c("Protein-coding derived cMAP"="#4DBBD5FF","Non-coding derived ncMAP"="#E64B35FF",
                          "Protein-coding derived ncMAP"="#FFC78B","Protein-coding derived MAP"="#F78E36"))+
  stat_compare_means(comparisons = mycompare,method="wilcox.test",label = "p.signif",
                     tip.length=0,label.x = 1.5,label.y = c(0.4,0.5,0.6,0.7,0.8,0.9))+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size =12,color="black"),
        axis.title.y = element_text(size = 14,color="black"),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size =0.5),
        legend.position = "top",
        legend.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank())


setwd("~/Desktop/22白血病/Runtime/2特征比较/10表达/")
pdf("2D.pdf",width = 18,height = 6)
p1
dev.off()

pdf("S2D.pdf",width = 10,height = 5)
p2
dev.off()

pdf("2E.pdf",width = 6,height = 5)
p3
dev.off()

pdf("2F.pdf",width = 6,height = 5)
p4
dev.off()

pdf("2G.pdf",width = 6,height = 5)
p5
dev.off()


