rm(list = ls()) 


setwd("C:/Users/Administrator/Desktop/SLE procedures/step4. network")

library(tidyverse)
library(dplyr)

data1 <-  data.table::fread("C:/Users/Administrator/Desktop/SLE procedures/step2. TCMSP&Pubchem&SwissADME/CoreDrug.TargetPrediction.csv",data.table=F)
  
core <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step1. core drug/results/core drug.txt",header = T, sep="\t")
data1 <- data1[which(data1$drug_English %in% core$drug),] 
table(data1$drug_English)

allTargets.symbol <- data1 %>% 
  select(Drug=drug_English, MolId=MOL_ID,MolName=Molecule,Symbol=Common_name)
write.table(file="./network/allTargets.symbol.txt", allTargets.symbol, sep="\t", quote=F, col.names=T, row.names=F)

#疾病相关基因
data2 <-  read.table("C:/Users/Administrator/Desktop/SLE procedures/step3. disease-related genes/disease-related genes.txt",header = T, sep="\t") %>% 
  select(Total) %>% distinct()
write.table(file="./network/Disease.txt", data2, sep="\t", quote=F, col.names=F, row.names=F)

jj <- list()
jj[["Drug targets"]]=unique(data1$Common_name)
jj[["Disease-related genes"]]=unique(data2$Total)


library(venn)
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#91612D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="diseaseRelated-DrugTargets.venn.pdf",width=5,height=5)
venn(jj,col=mycol[1:length(jj)],zcolor=mycol[1:length(jj)],box=F)
dev.off()

#保存疾病和药物的交集基因
intersectGenes=intersect(jj[["Drug targets"]],jj[["Disease-related genes"]])
write.table(file="./ppi/diseaseRelated-DrugTargets-intersectgenes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)

