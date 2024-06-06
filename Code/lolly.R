library(ggplot2)
library(dplyr)
data<-data.frame(Features=c("PEP","Andromeda Scores","MHCflurry Affinity HLA-A","MHCflurry Affinity HLA-B","MHCflurry Affinity HLA-C",
                            "MHCflurry %rank HLA-A","MHCflurry %rank HLA-B","MHCflurry %rank HLA-C","NetMHCpan Affinity HLA-A","NetMHCpan Affinity HLA-B","NetMHCpan Affinity HLA-C",
                            "NetMHCpan %rank HLA-A","NetMHCpan %rank HLA-B","NetMHCpan %rank HLA-C","MixMHCpred %rank HLA-A","MixMHCpred %rank HLA-B","MixMHCpred %rank HLA-C",
                            "DeepLC","Length","Length of source ORF","Hydrophobicity Fraction","T-cell contact","sample"),
                 Canonical=c(0.02,102.71,53.49,73.72,134.66,0.20,0.14,0.61,144.49,835.05,2873.76,0.53,0.49,1.18,0.08,0.09,0.59,0.97,9,657,0.44,-3.3,2),
                 Non_canonical=c(0.021,100.76,48.20,74.72,124.66,0.17,0.15,0.57,105.44,923.01,2563.81,0.41,0.50,1.12,0.07,0.10,0.59,0.96,9,222,0.44,-2.9,2))

data<-data[-c(18,23),]
data<-data[c(3,4,5,9,10,11,6,7,8,12,13,14,15,16,17,1,2,18,19,20,21),]
data$Difference_canonical<-(data$Canonical-data$Non_canonical)/data$Canonical

ggplot(data,aes(x=factor(Features,level=data$Features),y=Difference_canonical))+geom_segment(aes(x=Features,xend=Features,y=0,yend=Difference_canonical),size=1.5,color="#C9CACA",linetype="solid")+
  geom_point(size=5,color="#74C6BE",fill="#74C6BE",shape=21)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1.1,hjust = 1.2))+scale_y_continuous(limits = c(-0.2, 0.8),breaks = seq(0.2,0.8,0.2))