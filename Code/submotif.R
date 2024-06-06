library(ggplot2)
library(motifStack)
library(HDMD)
library(ecodist)
library(dbscan)
library(fpc)
library(factoextra)
getDistances <- function(peps, pos_weights=NULL, len=9) {
  if (is.null(pos_weights)) { pos_weights = rep(1, len) }
  n <- length(peps)
  dists <- matrix(nrow=n, ncol=n)
  . <<- lapply(1:n, function(i) { dists[i,i:n] <<- unlist(lapply(i:n, function(j) {
    pepDistPMBEC(strsplit(peps[i], '')[[1]], strsplit(peps[j], '')[[1]], len, pos_weights)
  }))})
  dists[lower.tri(dists)] <- t(dists)[lower.tri(dists)]
  colnames(dists) <- rownames(dists) <- peps
  return(dists)
}

pepDistPMBEC <- function(pepA, pepB, len, pos_weights) {
  return ((sum(unlist(lapply(1:len, function(i) {
    distPMBEC[pepA[i], pepB[i]]*pos_weights[i]
  })))) / len)
}

getLogo <- function(peptides) {
  freqs <- MolecularEntropy(peptides, type='AA')
  motif <- pcm2pfm(freqs$counts)
  motif <- new('pfm', mat=motif, name='',
               color=colorset(alphabet='AA', colorScheme='chemistry'))
  return (motif)
}
#==========1.数据加载==========#
setwd("~/Desktop/22白血病/Runtime/2特征比较/2.等位分组条形图/")
data<-read.table("MAP_Allele_support.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)
data$Length<-nchar(data$Epitopes)
data<-data[data$Length==9,]
data<-data[,c("Alleles","Epitopes")]
data<-data[!duplicated(data),]

Alleles<-unique(data$Alleles)
Alleles<-Alleles[order(Alleles)]

setwd("~/Desktop/22白血病/Runtime/4ncMAP呈递特征/1.submotif")
distPMBEC <- read.table('distPMBEC.txt', header=TRUE, stringsAsFactors=FALSE)

dist_output<-paste0(gsub("-|:","",Alleles),"dists.txt")
nmds_output<-paste0(gsub("-|:","",Alleles),"nmds.txt")
#=================HLA-A=================#
#==========1.HLA-A01:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[1])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[1],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.009, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[1],sep = "\t",quote = F,row.names = F)
#====================#


#==========2.HLA-A02:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[2])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[2],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.005, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[2],sep = "\t",quote = F,row.names = F)
#====================#


#==========3.HLA-A03:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[3])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[3],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.007, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[3],sep = "\t",quote = F,row.names = F)
#====================#


#==========4.HLA-A11:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[4])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[4],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.012, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[4],sep = "\t",quote = F,row.names = F)
#====================#


#==========5.HLA-A24:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[5])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[5],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0048, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[5],sep = "\t",quote = F,row.names = F)
#====================#


#==========6.HLA-A26:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[6])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[6],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 5)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0058, MinPts = 5)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[6],sep = "\t",quote = F,row.names = F)
#====================#


#==========7.HLA-A29:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[7])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[7],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0062, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[7],sep = "\t",quote = F,row.names = F)
#====================#


#==========8.HLA-A30:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[8])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[8],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.01, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[8],sep = "\t",quote = F,row.names = F)
#====================#


#==========9.HLA-A30:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[9])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[9],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0066, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[9],sep = "\t",quote = F,row.names = F)
#====================#


#==========10.HLA-A34:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[10])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[10],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0084, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[10],sep = "\t",quote = F,row.names = F)
#====================#


#==========11.HLA-A68:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[11])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[11],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0144, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[11],sep = "\t",quote = F,row.names = F)
#====================#

#=================HLA-B=================#
#==========12.HLA-B07:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[12])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[12],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0075, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[12],sep = "\t",quote = F,row.names = F)
#====================#


#==========13.HLA-B08:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[13])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[13],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.005, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[13],sep = "\t",quote = F,row.names = F)
#====================#


#==========14.HLA-B13:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[14])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[14],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.021, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[14],sep = "\t",quote = F,row.names = F)
#====================#


#==========15.HLA-B13:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[15])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[15],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0068, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[15],sep = "\t",quote = F,row.names = F)
#====================#


#==========16.HLA-B14:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[16])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[16],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.006, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[16],sep = "\t",quote = F,row.names = F)
#====================#


#==========17.HLA-B15:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[17])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[17],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.003, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[17],sep = "\t",quote = F,row.names = F)
#====================#


#==========18.HLA-B18:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[18])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[18],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.01, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[18],sep = "\t",quote = F,row.names = F)
#====================#


#==========19.HLA-B27:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[19])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[19],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0102, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[19],sep = "\t",quote = F,row.names = F)
#====================#


#==========20.HLA-B27:05的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[20])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[20],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.006, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[20],sep = "\t",quote = F,row.names = F)
#====================#


#==========21.HLA-B38:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[21])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[21],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.005, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[21],sep = "\t",quote = F,row.names = F)
#====================#


#==========22.HLA-B39:05的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[22])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[22],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0063, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[22],sep = "\t",quote = F,row.names = F)
#====================#


#==========23.HLA-B40:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[23])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[23],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.004, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[23],sep = "\t",quote = F,row.names = F)
#====================#


#==========24.HLA-B44:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[24])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[24],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.01, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[24],sep = "\t",quote = F,row.names = F)
#====================#


#==========25.HLA-B44:03的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[25])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[25],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0049, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[25],sep = "\t",quote = F,row.names = F)
#====================#


#==========26.HLA-B44:05的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[26])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[26],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.003, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[26],sep = "\t",quote = F,row.names = F)
#====================#


#==========27.HLA-B51:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[27])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[27],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0071, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[27],sep = "\t",quote = F,row.names = F)
#====================#


#==========28.HLA-B53:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[28])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[28],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0093, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[28],sep = "\t",quote = F,row.names = F)
#====================#


#==========29.HLA-B55:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[29])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[29],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.012, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[29],sep = "\t",quote = F,row.names = F)
#====================#


#==========30.HLA-B56:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[30])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[30],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0126, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[30],sep = "\t",quote = F,row.names = F)
#====================#


#==========31.HLA-B57:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[31])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[31],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0053, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[31],sep = "\t",quote = F,row.names = F)
#====================#


#==========32.HLA-B57:03的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[32])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[32],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0062, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[32],sep = "\t",quote = F,row.names = F)
#====================#


#=================HLA-C=================#
#==========33.HLA-C01:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[33])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[33],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0076, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[33],sep = "\t",quote = F,row.names = F)
#====================#


#==========34.HLA-C02:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[34])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[34],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0046, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[34],sep = "\t",quote = F,row.names = F)
#====================#


#==========35.HLA-C03:03的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[35])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[35],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0051, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[35],sep = "\t",quote = F,row.names = F)
#====================#


#==========36.HLA-C03:04的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[36])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[36],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0051, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[36],sep = "\t",quote = F,row.names = F)
#====================#


#==========37.HLA-C04:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[37])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[37],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0079, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[37],sep = "\t",quote = F,row.names = F)
#====================#


#==========38.HLA-C05:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[38])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[38],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.007, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[38],sep = "\t",quote = F,row.names = F)
#====================#


#==========39.HLA-C06:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[39])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[39],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0044, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[39],sep = "\t",quote = F,row.names = F)
#====================#


#==========40.HLA-C07:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[40])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[40],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.004, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[40],sep = "\t",quote = F,row.names = F)
#====================#


#==========41.HLA-C07:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[41])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[41],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.003, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[41],sep = "\t",quote = F,row.names = F)
#====================#


#==========42.HLA-C08:02的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[42])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[42],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0056, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[42],sep = "\t",quote = F,row.names = F)
#====================#


#==========43.HLA-C12:03的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[43])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[43],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0071, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[43],sep = "\t",quote = F,row.names = F)
#====================#


#==========44.HLA-C16:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[44])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[44],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.009, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[44],sep = "\t",quote = F,row.names = F)
#====================#


#==========45.HLA-C18:01的motif==========#
pep<-data$Epitopes[which(data$Alleles==Alleles[45])]
molecularEntropy<-MolecularEntropy(pep,type='AA')
pos_weights<-1-molecularEntropy$H
dists<-getDistances(pep,pos_weights)
write.table(dists,dist_output[45],sep = "\t",quote = F)
#=====降维=====#
rndseed = 12311
set.seed(rndseed)
nmds.tmp <- nmds(as.dist(dists), mindim=2, maxdim=2, nits=2)
nmds.dat <- nmds.tmp$conf[[which.min(nmds.tmp$stress)]]
nmds.dat<-data.frame(pep,nmds.dat)
colnames(nmds.dat)[2:3]<-c("nmds1","nmds2")
#=====聚类=====#
#选择最优的Eps值
set.seed(123456)
kNNdistplot(nmds.dat[,2:3],k = 3)
set.seed(123456)
Cluster <- fpc::dbscan(nmds.dat[,2:3], eps = 0.0062, MinPts = 3)
nmds.dat$Cluster<-paste0("Cluster",Cluster$cluster)
write.table(nmds.dat,nmds_output[45],sep = "\t",quote = F,row.names = F)
#====================#







