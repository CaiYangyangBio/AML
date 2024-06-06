library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(ggpubr)
library(factoextra)

#==========1.样本信息加载==========#
setwd("~/Desktop/22白血病/Runtime/0rawdata/")
Sample_info<-read_xlsx("PXD018542_AML.xlsx")
Sample_info$allele<-gsub("\\*","",Sample_info$allele)

Sample<-unique(Sample_info$sample)

Allele<-unique(unlist(strsplit(Sample_info$allele,",")))
Allele<-paste0("HLA-",Allele)

#==========2.样本不同等位呈递MAPs的比例:Sample_Allele_MAP==========#
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

#=====MAP=====#
Sample_Allele_MAP<-matrix(0,nrow = length(Allele),ncol = length(Sample))
rownames(Sample_Allele_MAP)<-Allele
colnames(Sample_Allele_MAP)<-Sample
for (i in 1:length(Sample)) {
  AML_MAP<-data[which(data$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Alleles))
  Allele_MAP$Freq<-round(Allele_MAP$Freq/MAP_num*100)
  index<-match(Allele,Allele_MAP$Var1,nomatch = 0)
  Sample_Allele_MAP[which(index!=0),i]<-Allele_MAP$Freq[index]
}

#=====cMAP=====#
cMAP<-data[data$Type=="Canonical",]
Sample_Allele_cMAP<-matrix(0,nrow = length(Allele),ncol = length(Sample))
rownames(Sample_Allele_cMAP)<-Allele
colnames(Sample_Allele_cMAP)<-Sample

for (i in 1:length(Sample)) {
  AML_MAP<-cMAP[which(cMAP$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Alleles))
  Allele_MAP$Freq<-round(Allele_MAP$Freq/MAP_num*100)
  index<-match(Allele,Allele_MAP$Var1,nomatch = 0)
  Sample_Allele_cMAP[which(index!=0),i]<-Allele_MAP$Freq[index]
}

#=====ncMAP=====#
ncMAP<-data[data$Type!="Canonical",]
Sample_Allele_ncMAP<-matrix(0,nrow = length(Allele),ncol = length(Sample))
rownames(Sample_Allele_ncMAP)<-Allele
colnames(Sample_Allele_ncMAP)<-Sample

for (i in 1:length(Sample)) {
  AML_MAP<-ncMAP[which(ncMAP$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Alleles))
  Allele_MAP$Freq<-round(Allele_MAP$Freq/MAP_num*100)
  index<-match(Allele,Allele_MAP$Var1,nomatch = 0)
  Sample_Allele_ncMAP[which(index!=0),i]<-Allele_MAP$Freq[index]
}

#====================#


#==========3.MAP的等位呈递频率:Allele_MAP_freq==========#

#=====MAP=====#
Allele_MAP_freq<-matrix(0,nrow = 6,ncol = length(Sample))
rownames(Allele_MAP_freq)<-c("1 Allele",paste0(2:6," Alleles"))
colnames(Allele_MAP_freq)<-Sample

for (i in 1:length(Sample)) {
  AML_MAP<-data[which(data$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Epitopes))
  MAP_Freq<-data.frame(table(Allele_MAP$Freq))
  MAP_Freq$Freq<-MAP_Freq$Freq/MAP_num*100
  index<-match(1:6,MAP_Freq$Var1,nomatch = 0)
  Allele_MAP_freq[which(index!=0),i]<-MAP_Freq$Freq[index]
}
Allele_freq_sum<-apply(Allele_MAP_freq, 1, sum)
Allele_MAP_freq<-Allele_MAP_freq[which(Allele_freq_sum!=0),]

#=====cMAP=====#
Allele_cMAP_freq<-matrix(0,nrow = 6,ncol = length(Sample))
rownames(Allele_cMAP_freq)<-c("1 Allele",paste0(2:6," Alleles"))
colnames(Allele_cMAP_freq)<-Sample

for (i in 1:length(Sample)) {
  AML_MAP<-cMAP[which(cMAP$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Epitopes))
  MAP_Freq<-data.frame(table(Allele_MAP$Freq))
  MAP_Freq$Freq<-MAP_Freq$Freq/MAP_num*100
  index<-match(1:6,MAP_Freq$Var1,nomatch = 0)
  Allele_cMAP_freq[which(index!=0),i]<-MAP_Freq$Freq[index]
}
Allele_freq_sum<-apply(Allele_cMAP_freq, 1, sum)
Allele_cMAP_freq<-Allele_cMAP_freq[which(Allele_freq_sum!=0),]

#=====ncMAP=====#
Allele_ncMAP_freq<-matrix(0,nrow = 6,ncol = length(Sample))
rownames(Allele_ncMAP_freq)<-c("1 Allele",paste0(2:6," Alleles"))
colnames(Allele_ncMAP_freq)<-Sample

for (i in 1:length(Sample)) {
  AML_MAP<-ncMAP[which(ncMAP$Samples == Sample[i]),]
  MAP_num<-length(unique(AML_MAP$Epitopes))
  Allele_MAP<-data.frame(table(AML_MAP$Epitopes))
  MAP_Freq<-data.frame(table(Allele_MAP$Freq))
  MAP_Freq$Freq<-MAP_Freq$Freq/MAP_num*100
  index<-match(1:6,MAP_Freq$Var1,nomatch = 0)
  Allele_ncMAP_freq[which(index!=0),i]<-MAP_Freq$Freq[index]
}
Allele_freq_sum<-apply(Allele_ncMAP_freq, 1, sum)
Allele_ncMAP_freq<-Allele_ncMAP_freq[which(Allele_freq_sum!=0),]
#====================#


#==========4.不同等位呈递的MAP在样本中的出现次数:Sample_MAP_freq==========#
#=====MAP=====#
Sample_MAP_freq<-matrix(0,nrow = length(Allele),ncol = 19)
rownames(Sample_MAP_freq)<-Allele
colnames(Sample_MAP_freq)<-c("1 Sample",paste0(2:19," Samples"))

for (i in 1:length(Allele)) {
  Allele_MAP<-data[which(data$Alleles == Allele[i]),]
  MAP_num<-length(unique(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Sample_Freq$Freq))
  Sample_Freq$Freq<-Sample_Freq$Freq/MAP_num*100
  index<-match(1:19,Sample_Freq$Var1,nomatch = 0)
  Sample_MAP_freq[i,which(index!=0)]<-Sample_Freq$Freq[index]
}
Sample_freq_sum<-apply(Sample_MAP_freq, 2, sum)
Sample_MAP_freq<-Sample_MAP_freq[,which(Sample_freq_sum!=0)]

#=====cMAP=====#
Sample_cMAP_freq<-matrix(0,nrow = length(Allele),ncol = 19)
rownames(Sample_cMAP_freq)<-Allele
colnames(Sample_cMAP_freq)<-c("1 Sample",paste0(2:19," Samples"))

for (i in 1:length(Allele)) {
  Allele_MAP<-cMAP[which(cMAP$Alleles == Allele[i]),]
  MAP_num<-length(unique(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Sample_Freq$Freq))
  Sample_Freq$Freq<-Sample_Freq$Freq/MAP_num*100
  index<-match(1:19,Sample_Freq$Var1,nomatch = 0)
  Sample_cMAP_freq[i,which(index!=0)]<-Sample_Freq$Freq[index]
}
Sample_freq_sum<-apply(Sample_cMAP_freq, 2, sum)
Sample_cMAP_freq<-Sample_cMAP_freq[,which(Sample_freq_sum!=0)]

#=====ncMAP=====#
Sample_ncMAP_freq<-matrix(0,nrow = length(Allele),ncol = 19)
rownames(Sample_ncMAP_freq)<-Allele
colnames(Sample_ncMAP_freq)<-c("1 Sample",paste0(2:19," Samples"))

for (i in 1:length(Allele)) {
  Allele_MAP<-ncMAP[which(ncMAP$Alleles == Allele[i]),]
  MAP_num<-length(unique(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Allele_MAP$Epitopes))
  Sample_Freq<-data.frame(table(Sample_Freq$Freq))
  Sample_Freq$Freq<-Sample_Freq$Freq/MAP_num*100
  index<-match(1:19,Sample_Freq$Var1,nomatch = 0)
  Sample_ncMAP_freq[i,which(index!=0)]<-Sample_Freq$Freq[index]
}
Sample_freq_sum<-apply(Sample_ncMAP_freq, 2, sum)
Sample_ncMAP_freq<-Sample_ncMAP_freq[,which(Sample_freq_sum!=0)]
#====================#


#==========5.不同等位呈递的MAP在数目：Allele_MAP_num==========#
#=====MAP=====#
Allele_MAP<-data[,c("Epitopes","Alleles")]
Allele_MAP<-Allele_MAP[!duplicated(Allele_MAP),]
Allele_MAP_num<-data.frame(table(Allele_MAP$Alleles))

#=====cMAP=====#
Allele_cMAP<-data[which(data$Type=="Canonical"),c("Epitopes","Alleles")]
Allele_cMAP<-Allele_cMAP[!duplicated(Allele_cMAP),]
Allele_cMAP_num<-data.frame(table(Allele_cMAP$Alleles))

#=====ncMAP=====#
Allele_ncMAP<-data[which(data$Type!="Canonical"),c("Epitopes","Alleles")]
Allele_ncMAP<-Allele_ncMAP[!duplicated(Allele_ncMAP),]
Allele_ncMAP_num<-data.frame(table(Allele_ncMAP$Alleles))

#====================#


setwd("~/Desktop/22白血病/Runtime/3样本等位相似性/2.样本MAP")
write.table(Sample_Allele_MAP,"Sample_Allele_MAP.txt",quote = F,sep = "\t")
write.table(Sample_Allele_cMAP,"Sample_Allele_cMAP.txt",quote = F,sep = "\t")
write.table(Sample_Allele_ncMAP,"Sample_Allele_ncMAP.txt",quote = F,sep = "\t")

write.table(Allele_MAP_freq,"Allele_MAP_freq.txt",quote = F,sep = "\t")
write.table(Allele_cMAP_freq,"Allele_cMAP_freq.txt",quote = F,sep = "\t")
write.table(Allele_ncMAP_freq,"Allele_ncMAP_freq.txt",quote = F,sep = "\t")

write.table(Sample_MAP_freq,"Sample_MAP_freq.txt",quote = F,sep = "\t")
write.table(Sample_cMAP_freq,"Sample_cMAP_freq.txt",quote = F,sep = "\t")
write.table(Sample_ncMAP_freq,"Sample_ncMAP_freq.txt",quote = F,sep = "\t")

write.table(Allele_MAP_num,"Allele_MAP_num.txt",quote = F,sep = "\t")
write.table(Allele_cMAP_num,"Allele_cMAP_num.txt",quote = F,sep = "\t")
write.table(Allele_ncMAP_num,"Allele_ncMAP_num.txt",quote = F,sep = "\t")






