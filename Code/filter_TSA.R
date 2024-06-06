library(survival)
library(survminer)
LAML<-read.table("~/Desktop/22白血病/Runtime/0rawdata/LAML.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
#=======百分之九十分位数大于1
per90<-apply(LAML, 1, function(x){quantile(x,0.9)})
index<-which(per90>1)
LAML<-LAML[index,]


setwd("~/Desktop/22白血病/Runtime/5识别新抗原/1.90percentile/1Raw")
#========差异表达的
DEG<-read.table("Different.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
FC<-DEG[,grepl("logFC",colnames(DEG))]
FC_res<-apply(FC,1,function(x){length(which((x>1)==T))})

P<-DEG[,grepl("adj.P.Val",colnames(DEG))]
P_res<-apply(P,1,function(x){length(which((x<0.05)==T))})

index<-which(FC_res==30 & P_res==30)
ae<-rownames(DEG[index,])


GTEx<-read.table("GTEx_ALL.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
#=======百分之九十分位数小于1
num<-apply(GTEx, 1, function(x){length(which(x<1))})
index<-which(num==30)
TSA<-intersect(rownames(LAML),rownames(GTEx[index,]))

#======aeTSA
aeTSA<-intersect(ae,TSA)

#======TSA
TSA<-setdiff(TSA,aeTSA)

#======CTA
num<-apply(GTEx, 1, function(x){length(which(x[c(1:26,28:30)]<1))})
index<-which(num==29)
CTA<-intersect(rownames(LAML),rownames(GTEx[index,]))
CTA<-setdiff(CTA,union(TSA,aeTSA))

#======MAP
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
MAP<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,quote = "",check.names = F,fill = T)

#======MAP_gene
immune<-read.table("~/Desktop/22白血病/Runtime/1预处理/All_info.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
MAP_Gene<-merge(MAP,immune)

#======Gene
setwd("~/Desktop/22白血病/Runtime/0rawdata/1-ncORF/1-SeqName")
Gene<-read.table("uniprot.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
Gene$Name<-gsub(">","",Gene$Name)
Gene$RawName<-unlist(lapply(strsplit(Gene$RawName,"GN="),function(x){x[2]}))
Gene<-Gene[!is.na(Gene$RawName),]
Gene$RawName<-unlist(lapply(strsplit(Gene$RawName," "),function(x){x[1]}))

#======Canonical
Canonical<-MAP_Gene[which(MAP_Gene$Type=="Canonical"),]
index<-match(Canonical$`Leading razor protein`,Gene$Name,nomatch = 0)
Canonical<-Canonical[which(index!=0),]
Canonical$Symbol<-Gene$RawName[index]

#======Non-Canonical
NonCanonical<-MAP_Gene[which(MAP_Gene$Type!="Canonical"),]
NonCanonical$Symbol<-unlist(lapply(strsplit(NonCanonical$`Leading razor protein`,"_"),function(x){x[1]}))

#======pep
data<-rbind(Canonical,NonCanonical)
aeTSA<-data[data$Symbol%in%aeTSA,]
aeTSA$Class<-"aeTSA"

TSA<-data[data$Symbol%in%TSA,]
TSA$Class<-"TSA"

CTA<-data[data$Symbol%in%CTA,]
CTA$Class<-"CTA"
Epitopes<-rbind(aeTSA,TSA,CTA)

EpiGene<-Epitopes[,c("Epitopes","Symbol","Class","Type")]
EpiGene<-EpiGene[!duplicated(EpiGene),]

setwd("~/Desktop/22白血病/Runtime/5识别新抗原/1.90percentile")
write.table(EpiGene,"Epitopes.txt",quote = F,sep = "\t",row.names = F)
write.table(Epitopes,"Epitopes_info.txt",quote = F,sep = "\t",row.names = F)


