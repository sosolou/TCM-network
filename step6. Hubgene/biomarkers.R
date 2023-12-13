rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/SLE procedures/step6. Hubgene")

library(tidyverse)
library(dplyr)
library(limma)
library(reshape2)
library(ggplot2)
#轻型核心药物靶基因和疾病相关基因的交集基因
data1 <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step4. network/ppi/diseaseRelated-DrugTargets-intersectgenes.txt",header = F, sep="\t")

#GEO差异基因
data2 <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/results/diff.txt",header = T, sep="\t") 

jj <- list()
jj[["DrugTargets & DiseaseRelatedGenes"]]=unique(data1$V1)
jj[["GEO DEGs"]]=unique(data2$id)

library(venn)
mycol=c("#5D90BA","#D8D155")
pdf(file="./results/Figure 1.venn.pdf",width=5,height=5)
venn(jj,col=mycol,zcolor=mycol,box=F)
dev.off()

intersect_Genes=data2[which(data2$id %in% intersect(jj[[1]],jj[[2]])),] %>% arrange(desc(logFC))
write.table(file="./results/intersect_Genes.txt", intersect_Genes, sep="\t", quote=F, col.names=T, row.names=F)

table(intersect_Genes$regulate)

GO <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/results/GO.txt",sep="\t",header=T)

intersect_GO <- GO[grepl(paste(intersect_Genes$id, collapse = "|"), GO$geneID), ]
write.table(file="./results/intersect_GO.csv", intersect_GO, sep=",", quote=F, col.names=T, row.names=F)
nrow(GO)
nrow(intersect_GO)


#交集基因热图
exp <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/results/diffGeneExp.txt",header = T, sep="\t") %>%
  filter(id %in% intersect_Genes$id) %>%
  column_to_rownames("id") %>%
  as.data.frame()

group <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/group.txt",header = T, sep="\t") %>%
  column_to_rownames("GSM")

library(pheatmap)
pdf(file="./results/Figure 2.heatmap.pdf", width=8, height=6)
pheatmap(exp, 
         annotation=group, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()

#细胞死亡基因
celldeath <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/ref/celldeath.csv",sep=",",header=T,check.names=F) 
celldeath <- bind_rows(data.frame(ferr=celldeath[which(celldeath$ferrgene %in% intersect_Genes$id),]$ferrgene),
                          data.frame(pyro=celldeath[which(celldeath$pyrogene %in% intersect_Genes$id),]$pyrogene),
                          data.frame(anoikis=celldeath[which(celldeath$anoikisgene %in% intersect_Genes$id),]$anoikisgene),
                          data.frame(copper=celldeath[which(celldeath$coppergene %in% intersect_Genes$id),]$coppergene) )

celldeath <- as.data.frame(lapply(celldeath,function(x)c(x[!is.na(x)],x[is.na(x)])))
celldeath <- celldeath[rowSums(is.na(celldeath))!=ncol(celldeath),]
celldeath[is.na(celldeath)]<-""

write.table(file="./results/intersect_Genes.celldeath.txt", celldeath, sep="\t", quote=F, col.names=T, row.names=F)

#与免疫浸润细胞相关性
group_t <- group %>%
  filter(group=="Disease")
data <- exp[which(rownames(exp) %in% intersect_Genes$id),
               which(colnames(exp) %in% rownames(group_t))]  
data=t(data)

#读取免疫细胞结果文件，并对数据进行整理
immune.diff =read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/results/immuneDiff.xls",sep="\t",header=T,  check.names=F)
immune=read.table("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO/results/CIBERSORT-Results.txt" , header=T, sep="\t", check.names=F, row.names=1)
immune=immune[,which(colnames(immune) %in% immune.diff$Cell)]
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]

#相关性分析
outTab=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
  }
}

#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="./results/Figure 3.immun-cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    #去掉背景
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x轴字体
        axis.text.y = element_text(size = 8, face = "bold")) +       #y轴字体
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","correlation")) +   #设置图例
  scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()

#PPI交集基因
ppi <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step4. network/ppi/score2.out.gene.txt",sep="\t",header=F,  check.names=F)
hub <- intersect(intersect_Genes$id,ppi$V1)
hub
biomarker <- intersect_Genes[which(intersect_Genes$id %in% hub),]
write.table(file="./results/biomarkers.txt", biomarker, sep="\t", quote=F, col.names=T, row.names=F)







