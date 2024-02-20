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

setwd("~/Desktop/22白血病/Runtime/2特征比较/5.PEP_Score")
pdf("1C.pdf",width = 6,height = 5)
g1
dev.off()

#==============================#

