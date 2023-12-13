rm(list = ls()) #清除缓存
#修改1：设定工作目录
setwd("C:/Users/Administrator/Desktop/SLE procedures/step7.vina")

library(tidyverse)
library(dplyr)
hub <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step6. Hubgene/results/biomarkers.txt",header = T, sep="\t")


drug <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step1. core drug/results/core drug.txt",header = T,sep="\t")


SwissPrediction <-  data.table::fread("C:/Users/Administrator/Desktop/SLE procedures/step2. TCMSP&Pubchem&SwissADME/CoreDrug.TargetPrediction.csv",data.table=F) %>%
  mutate(id=paste0(paste0(MOL_ID,Molecule),Common_name))

TCMSP <- data.table::fread("C:/Users/Administrator/Desktop/SLE procedures/step2. TCMSP&Pubchem&SwissADME/TCMSP.csv",data.table=F,header = T,drop=1) %>%
  rename(drug=drug_English,drug_py=drug,MolId=`MOL ID`,MolName=MOlecule) %>%
  mutate(id=paste0(paste0(drug,MolId),MolName)) %>%
  select(id,MOlecule_url,InChIKey)

#network
vina <- data.table::fread("C:/Users/Administrator/Desktop/SLE procedures/step4. network/network/net.network.txt",data.table=F)  %>%
  filter(Node2 %in% hub$id) %>%
  arrange(Node2,Node1) %>%
  mutate(id=paste0(paste0(Node1,MolName),Node2)) %>%
  left_join(SwissPrediction,by="id") %>%
  distinct(Node1,Node2,MolName, .keep_all = TRUE) %>%
  select(Drug=drug_English,drug_py=drug,MolId=Node1,MolName,Target=Node2,TargetFullName=Target,Uniprot_ID) %>%
  mutate(id=paste0(paste0(drug,MolId),MolName)) %>%
  left_join(TCMSP,by="id") %>%
  select(Drug,drug_py,MolId,MolName,MOlecule_url,InChIKey,Target,TargetFullName,Uniprot_ID)  %>%
  filter(Drug %in% drug$drug) %>%
  distinct() %>%
  arrange(Target,Drug,MolId,MolName)

table(vina$drug_py) %>% as.data.frame()
table(vina$drug) %>% as.data.frame()
table(vina$Target) %>% as.data.frame()


write.table(file="vina.txt", vina, sep="\t", quote=F, col.names=T, row.names=F)

