setwd("~/Desktop/22数据库/数据/原始数据/结果文件/GTEx/")
GTEx<-list.files(pattern = "txt")
GTEx<-GTEx[-29]
GTEx_tissue<-gsub(".txt","",GTEx)

GTEx_m<-list()
for (i in 1:length(GTEx)) {
  Tissue<-read.table(GTEx[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
  GTEx_m[[i]]<-apply(Tissue, 1, function(x){quantile(x,0.9)})
}
GTEx_res<-do.call("rbind",GTEx_m)
GTEx_res<-t(GTEx_res)
colnames(GTEx_res)<-GTEx_tissue

setwd("~/Desktop/22白血病/Runtime/5识别新抗原/1.90percentile")
write.table(GTEx_res,"GTEx_ALL.txt",quote = F,sep = "\t")
