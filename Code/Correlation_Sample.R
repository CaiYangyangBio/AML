library(ggplot2)
library(ggpubr)
library(factoextra)
library(reshape2)
library(readxl)
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


#==========1.样本MAP的辛普森系数==========#
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",check.names = F,fill = T)
Sample<-unique(MAP$Samples)

#=====All=====#
Simpson_A<-matrix(0,nrow = length(Sample),ncol = length(Sample))
colnames(Simpson_A)<-Sample
rownames(Simpson_A)<-Sample
for (i in 1:length(Sample)) {
  A_MAP<-unique(MAP$Epitopes[which(MAP$Samples==Sample[i])])
  for (j in 1:length(Sample)) {
    B_MAP<-unique(MAP$Epitopes[which(MAP$Samples==Sample[j])])
    n<-length(intersect(A_MAP,B_MAP))
    m<-min(length(A_MAP),length(B_MAP))
    Simpson_A[i,j]<-round(n/m,2)
  }
}

Simpson_All<-get_upper_tri(Simpson_A)
Simpson_All<- melt(as.matrix(Simpson_All), na.rm = T)
colnames(Simpson_All)<-c("Sample1","Sample2","Simpson_All")

#=====经典=====#
C<-MAP[which(MAP$Type=="Canonical"),]
Simpson_C<-matrix(0,nrow = length(Sample),ncol = length(Sample))
colnames(Simpson_C)<-Sample
rownames(Simpson_C)<-Sample
for (i in 1:length(Sample)) {
  A_MAP<-unique(C$Epitopes[which(C$Samples==Sample[i])])
  for (j in 1:length(Sample)) {
    B_MAP<-unique(C$Epitopes[which(C$Samples==Sample[j])])
    n<-length(intersect(A_MAP,B_MAP))
    m<-min(length(A_MAP),length(B_MAP))
    Simpson_C[i,j]<-round(n/m,2)
  }
}

# g1<-fviz_nbclust(Simpson_C, kmeans, method = "wss") + 
#   geom_vline(xintercept = 3, linetype = 2)
# 
# km_result <- kmeans(t(Simpson_C), centers = 3, nstart = 24)
# 
# g2<-fviz_cluster(km_result, data = t(Simpson_C),
#                  palette = c("#DC0000FF","#3C5488FF","#00A087FF"),
#                  ellipse.type = "euclid",main = "Canonical MAPs",
#                  star.plot = TRUE, 
#                  repel = TRUE,
#                  ggtheme = theme_bw())
# 
# p1<-ggarrange(g1,g2,ncol = 2)

Simpson_Canonical <- get_upper_tri(Simpson_C)
Simpson_Canonical <- melt(as.matrix(Simpson_Canonical), na.rm = T)
colnames(Simpson_Canonical)<-c("Sample1","Sample2","Simpson_C")

#=====非经典=====#
N<-MAP[which(MAP$Type!="Canonical"),]
Simpson_N<-matrix(0,nrow = length(Sample),ncol = length(Sample))
colnames(Simpson_N)<-Sample
rownames(Simpson_N)<-Sample
for (i in 1:length(Sample)) {
  A_MAP<-unique(N$Epitopes[which(N$Samples==Sample[i])])
  for (j in 1:length(Sample)) {
    B_MAP<-unique(N$Epitopes[which(N$Samples==Sample[j])])
    n<-length(intersect(A_MAP,B_MAP))
    m<-min(length(A_MAP),length(B_MAP))
    Simpson_N[i,j]<-round(n/m,2)
  }
}

# g3<-fviz_nbclust(Simpson_N, kmeans, method = "wss") + 
#   geom_vline(xintercept = 3, linetype = 2)
# 
# km_result <- kmeans(t(Simpson_N), centers = 3, nstart = 24)
# 
# g4<-fviz_cluster(km_result, data = t(Simpson_N),
#                  palette = c("#DC0000FF","#3C5488FF","#00A087FF"),
#                  ellipse.type = "euclid",
#                  star.plot = TRUE, main = "Non-Canonical MAPs",
#                  repel = TRUE,
#                  ggtheme = theme_bw())
# 
# p2<-ggarrange(g3,g4,ncol = 2)
# p<-ggarrange(p1,p2,ncol = 1)

Simpson_Non <- get_upper_tri(Simpson_N)
Simpson_Non <- melt(as.matrix(Simpson_Non), na.rm = T)
colnames(Simpson_Non)<-c("Sample1","Sample2","Simpson_N")


#==========2.不同样本共享的等位数==========#
setwd("~/Desktop/22白血病/Runtime/0rawdata/")
sample_info<-read_xlsx("PXD018542_AML.xlsx")
sample_info$allele<-gsub("\\*","",sample_info$allele)

Allele<-unlist(strsplit(sample_info$allele,","))
Allele<-paste0("HLA-",Allele)
Sample_Allele<-data.frame(Sample = rep(sample_info$sample,each = 6),Allele)
Sample_Allele<-Sample_Allele[!duplicated(Sample_Allele),]

Sample_Intersect<-matrix(0,nrow = length(Sample),ncol = length(Sample))
rownames(Sample_Intersect)<-Sample
colnames(Sample_Intersect)<-Sample

for (i in 1:length(Sample)) {
  A_Allele<-Sample_Allele$Allele[which(Sample_Allele$Sample==Sample[i])]
  for (j in 1:length(Sample)) {
    B_Allele<-Sample_Allele$Allele[which(Sample_Allele$Sample==Sample[j])]
    Sample_Intersect[i,j]<-length(intersect(A_Allele,B_Allele))
  }
}

Share_Allele <- get_upper_tri(Sample_Intersect)
Share_Allele <- melt(as.matrix(Share_Allele), na.rm = T)
colnames(Share_Allele)<-c("Sample1","Sample2","Share_Allele")


res<-merge(Simpson_All,Simpson_Canonical)
res<-merge(res,Simpson_Non)
res<-merge(res,Share_Allele)
index<-apply(res, 1, function(x){x[1]!=x[2]})
res<-res[index,]
res<-res[order(res$Share_Allele,res$Simpson_N,decreasing = T),]

Sample1<-as.numeric(gsub("AML","",res$Sample1))
Sample2<-as.numeric(gsub("AML","",res$Sample2))
index<-Sample1<Sample2
res[which(index==F),1:2]<-res[which(index==F),2:1]

setwd("~/Desktop/22白血病/Runtime/3样本等位相似性/1.样本相似性")
write.table(Simpson_A,"Simpson_All.txt",quote = F,sep = "\t")
write.table(Simpson_C,"Simpson_C.txt",quote = F,sep = "\t")
write.table(Simpson_N,"Simpson_N.txt",quote = F,sep = "\t")
write.table(Sample_Intersect,"Shared_Allele.txt",quote = F,sep = "\t")
write.table(res,"result.txt",quote = F,sep = "\t",row.names = F)

# pdf("聚类.pdf",width = 12,height = 10)
# p
# dev.off()








