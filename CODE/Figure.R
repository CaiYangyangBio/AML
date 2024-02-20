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



