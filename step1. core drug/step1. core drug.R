rm(list = ls())
setwd("C:/Users/Administrator/Desktop/SLE procedures/step1. core drug")

library(tidyverse)
library(dplyr)
library(arules) 
library(arulesViz) 
library(RColorBrewer)
library(shinythemes)
library(showtext)
library(cowplot)
library(ape)
library(factoextra)
showtext_auto(enable=T)

pre <- read.table("drug.csv",sep=",",header=T,check.names=F, row.names = 1)
N <- nrow(pre)

transdata0 <- read.transactions("drug.csv", format =c("basket"), 
                                sep = ",", cols = 1, header = TRUE, rm.duplicates = TRUE)

herbfreq <- as.data.frame(sort(itemFrequency(transdata0, "absolute"), decreasing = TRUE)) 
colnames(herbfreq)<-"freq"

herbfreq$'percent(%)' <- round(herbfreq$freq/N*100,digits = 2) 
herbfreq <- herbfreq %>% rownames_to_column("drug")

if(N > 50) {
  top_n <- round(N * 0.2, -1)
} else {
  top_n <- round(N * 0.2, 0)
}

top <- nrow(herbfreq[which(herbfreq$freq>=top_n),])
hmod <- herbfreq[which(herbfreq$freq>=top_n),]$drug

pdf("./results/Figure 1. Distribution chart of high-frequency TCM.pdf", height = 8, width = 14)
itemFrequencyPlot(transdata0, topN = top,
                        type ="absolute",
                        horiz = F, col = "#B0C4DE",
                        xlab="Herb",ylim=c(0,round(herbfreq[1,2],digits = -1)),
                        cex.names=0.8, cex.axis=0.9,
                        ylab="Frequency"
                        
)
dev.off()


strong_rules <-  apriori(transdata0, parameter = list( support = 0.1, confidence = 0.8,
                                                       maxlen = 10, minlen = 2, target = "rules"))
strong_rules <- sort(strong_rules, by = "support", decreasing = T ) 
strong_rules <-  subset(strong_rules, lift >= 3)
rules_s <- as(strong_rules,'data.frame') %>% as_tibble() %>%
  separate(rules, c("lhs","rhs"), sep = "=>") 

write.table(rules_s,file="./results/strong rules.txt",sep="\t",quote=F,col.names=T,row.names = T)

smod=unique(
  trimws(
    c(
      unique(unlist(strsplit(gsub("\\{|\\}","",rules_s$lhs),","))),
      unique(unlist(strsplit(gsub("\\{|\\}","",rules_s$rhs),",")))
    )))


pdf("./results/Figure 2. strong rules .pdf", height = 8, width = 14) 
par(cex = 1.5)
plot(strong_rules,method = "graph",measure = "confidence",
     control = list(shading="lift"),colors=c("#76EEC6","#FAEBD7"))
dev.off()


drug <- read.table("drug.csv",sep=",",header=T,check.names=F, row.names = 1)  %>% t()  %>% as_tibble() %>%
  sapply(table) %>% unlist() 

b <- names(drug) %>% as_tibble() %>% 
  separate(value, c("prescription", "drug"),sep="\\.") %>%
  filter(drug != " " ) %>%
  mutate(count=1) %>%
  spread(drug,count) 

b[is.na(b)]<-0
b<-data.frame(b,stringsAsFactors = F)  %>% column_to_rownames("prescription")

colsum <- colSums(b)

c<- rbind(b,colsum) 
cc <- c[,which((c[nrow(c),] > ceiling(N*0.05)))] 
database <- cc[-nrow(cc),] %>% t() 
rownames(database)=str_replace_all(rownames(database),"\\W+"," ")



nb1 <- fviz_nbclust(database, hcut, method = "wss")+
  labs(subtitle = "Elbow method",title = "")

nb2 <- fviz_nbclust(database, hcut, method = "silhouette")+
  labs(subtitle = "Silhouette method",title = "")


nb3 <- fviz_nbclust(database, hcut, method = "gap_stat")+
  labs(subtitle = "Gap statistic method",title = "")

pdf("./results/please choose a best K-value for cluster analysis.pdf", height = 4, width = 12)
plot_grid(nb1,nb2,nb3,
          nrow=1)
dev.off()


#请根据输出图形输入最佳K值，一般在2~10之间
best_K <- 3

nb1 <- fviz_nbclust(database, hcut, method = "wss")+
  geom_vline(xintercept = best_K, linetype = 2)+
  labs(subtitle = "Elbow method",title = "")

nb2 <- fviz_nbclust(database, hcut, method = "silhouette")+
  geom_vline(xintercept = best_K, linetype = 2)+
  labs(subtitle = "Silhouette method",title = "")


nb3 <- fviz_nbclust(database, hcut, method = "gap_stat")+
  geom_vline(xintercept = best_K, linetype = 2)+
  labs(subtitle = "Gap statistic method",title = "")

pdf("./results/Figure 3. best K-value for cluster analysis.pdf", height = 4, width = 12)
plot_grid(nb1,nb2,nb3,
          nrow=1)
dev.off()

hc <- hclust(dist(database,method = "binary"),method = "complete")
kmod <- as.data.frame(cutree(hc, k=best_K)) 
colnames(kmod)="cluster"
kmod <- kmod %>% rownames_to_column("drug") %>% arrange(cluster) 



mypal = c("#556270","#4ECDC4","#FF6B6B","lightgoldenrod","lightgreen","rosybrown","skyblue","violet","peachpuff","yellowgreen","paleturquoise")
clus3 = cutree(hc,best_K) 
pdf("./results/Figure 4. Cluster Dendrogram.pdf", height = 10, width = 10) 
plot(as.phylo(hc),cex=0.8,type="cladogram",tip.color = mypal[clus3],  col = "red")
dev.off()

jmod=kmod[which(kmod$cluster==1),]$drug

cmod <- intersect(jmod,intersect(smod,hmod)) %>% as.data.frame() 
colnames(cmod) <- "drug"
cmod1<-merge(cmod,herbfreq,by = 'drug',all.x=TRUE,all.y=FALSE) %>% arrange(desc(freq))
cmod1
write.table(cmod1,file="./results/core drug.txt",sep="\t",quote=F,col.names=T,row.names = T)


library(ggVennDiagram)
hx <- list(
  A = hmod, 
  B = smod,
  C = jmod
)

pdf("./results/Figure 5.core drug.pdf", height = 8, width = 8) 
ggVennDiagram(hx, 
                    label_alpha = 0,  
                    label = "count", 
                    edge_size = 0.1, edge_lty = "solid",
                    category.names = c(paste0("Inclusion criteria 1:","\n","High-frequency herbs"),paste0("Inclusion criteria 2:","\n","In accord with strong association rules"),paste0("Inclusion criteria 3:","\n","In maximum cluster clustering")),
              set_size = 4)+
  scale_fill_distiller(palette = "Blues", direction = -3) +
  scale_color_brewer(palette = "Paired")+
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  guides(fill=FALSE) 
dev.off()




