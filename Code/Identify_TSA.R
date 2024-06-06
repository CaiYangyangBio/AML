library(limma)
LAML<-read.table("~/Desktop/22白血病/Runtime/0rawdata/LAML.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

setwd("~/Desktop/22数据库/数据/原始数据/结果文件/GTEx/")
GTEx<-list.files(pattern = "txt")
GTEx<-GTEx[-29]
GTEx_tissue<-gsub(".txt","",GTEx)
res<-list()
for (i in 1:length(GTEx)) {
  Tissue<-read.table(GTEx[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  index<-match(rownames(LAML),rownames(Tissue),nomatch = 0)
  LAML<-LAML[which(index!=0),]
  Tissue<-Tissue[index,]
  data<-cbind(LAML,Tissue)
  
  group_list<-c(rep("LAML",ncol(LAML)),rep("GTEx",ncol(Tissue)))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(data)
  ##这个矩阵声明，我们要把High组跟Low进行差异分析比较
  contrast.matrix<-makeContrasts(paste0(c("LAML","GTEx"),collapse = "-"),levels = design)
  fit <- lmFit(data,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  index<-match(rownames(data),rownames(tempOutput),nomatch = 0)
  tempOutput<-tempOutput[index,]
  colnames(tempOutput)<-paste0(GTEx_tissue[i],"_",colnames(tempOutput))
  res[[i]]<-tempOutput[,c(1,3,4)]
  print(i)
}

res<-do.call("cbind",res)
rownames(res)<-rownames(data)
setwd("~/Desktop/22白血病/Runtime/5识别新抗原/1.90percentile/1Raw")
write.table(res,"Different.txt",quote = F,sep = "\t")






