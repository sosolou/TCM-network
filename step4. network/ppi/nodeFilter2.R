rm(list = ls()) #清除缓存
setwd("C:/Users/Administrator/Desktop/SLE/procedure/step7. network/mild/ppi")    #设置工作目录

library(dplyr)
library(tidyverse)
score2= read.table("score2.csv",sep=",",header=T,check.names=F) %>%
  select(name,Betweenness,Closeness,Degree,Eigenvector,LAC,Network) %>%
  column_to_rownames("name")
nCol=ncol(score2)      #多少个过滤条件

nRow=nrow(score2) #网络有多少个节点
nBian=sum(score2$Degree)/2 #网络有多少条边

#对score2每个过滤条件循环，找出每个条件大于中位值的基因
mat=data.frame()      #新建数据框
for(i in colnames(score2)){
  print(paste0("filter1,",i,": ",median(score2[,i])))
  value=ifelse(score2[,i]>median(score2[,i]),1,0)
  mat=rbind(mat,value)
}
mat=t(mat)
colnames(mat)=colnames(score2)
row.names(mat)=row.names(score2)
#统计score2哪些基因所有条件都大于中位值,结果保存在score2
geneName=row.names(mat[rowSums(mat)==nCol,])
score2=score2[geneName,]
score2out=cbind(name=row.names(score2),score2)
write.table(score2out, file="score2.out.txt", sep='\t', quote=F, row.names=F)
write.table(geneName, file="score2.out.gene.txt", sep='\t', quote=F, row.names=F,col.names=F)


