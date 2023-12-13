rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/SLE procedures/step5. GEO")
#定义显著性
logFoldChange=0.58
adjustP=0.05

library(Biobase)
library(limma) 
library(stringr)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(tibble)
library(reshape)
library(vioplot) 

group <- read.table("group.txt",sep="\t",header=T,check.names=F) %>% arrange(group)
table(group$group)

exp <- read.table("GEOexp.txt",sep="\t",header=T,check.names=F,row.names = 1) 
order_vector <- as.vector(group$GSM)
order_index <- match(order_vector, colnames(exp))
exp <- exp[, order_index]

rt<-as.matrix(exp)
class <- as.character(group$group)
design <- model.matrix(~0+factor(class)) 
colnames(design) <- c("Control","Disease")
table(group$group)

col <- c(rep("skyblue",table(group$group)[1]),rep("pink",table(group$group)[2]))
pdf("./results/Figure1. Homogeneity of variance in sample data.pdf", height = 6, width = 6) 
par(mar=c(6,5,4,3) + 0.2)
boxplot(exp,outline=FALSE,notch=T,col=col,las=2,cex.lab=0.5)
dev.off()

fit <- lmFit(rt,design)
cont<-makeContrasts(Disease-Control,levels=design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
allDiff <- topTable(fit2,adjust='fdr',number=200000)

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSig <- diffSig[order(-diffSig$logFC),] %>%
  mutate(regulate=case_when(logFC>0~"up",logFC<0~"down"))
table(diffSig$regulate)
diffSigOut <- rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="./results/diff.txt", sep="\t", quote=F, col.names=F)

diffGeneExp<- rt[rownames(diffSig),]
diffGeneExpOut<- rbind(id=colnames(diffGeneExp),diffGeneExp)
write.table(diffGeneExpOut,file="./results/diffGeneExp.txt",sep="\t",quote=F,col.names=F)


Sig=ifelse((allDiff$adj.P.Val<adjustP) & (abs(allDiff$logFC)>logFoldChange), ifelse(allDiff$logFC>logFoldChange,"Up","Down"), "Stable")
hs=cbind(allDiff, Sig=Sig)
p <- ggplot(
  hs, aes(x = logFC, y = -log10(adj.P.Val), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-logFoldChange,logFoldChange),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(adjustP),lty=4,col="black",lwd=0.8) +
  labs(x="log2(Fold Change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf(file="./results/Figure2. Volcano plot of differentially expressed genes.pdf", width=6, height=6)
print(p)
dev.off()

rt=diffSig %>% rownames_to_column("id")

pvalueFilter=0.05       
qvalueFilter=0.05       
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        


kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="./results/GO.txt",sep="\t",quote=F,row.names = F)


ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

pdf("./results/Figure3. GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                    })
circos.clear()

middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)

main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)

logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


celldeath <- read.table("./ref/celldeath.csv",sep=",",header=T,check.names=F) 
celldeath <- bind_rows(data.frame(ferr=celldeath[which(celldeath$ferrgene %in% rownames(diffSig)),]$ferrgene),
                      data.frame(pyro=celldeath[which(celldeath$pyrogene %in% rownames(diffSig)),]$pyrogene),
                      data.frame(anoikis=celldeath[which(celldeath$anoikisgene %in% rownames(diffSig)),]$anoikisgene),
                      data.frame(copper=celldeath[which(celldeath$coppergene %in% rownames(diffSig)),]$coppergene) )

celldeath <- as.data.frame(lapply(celldeath,function(x)c(x[!is.na(x)],x[is.na(x)])))
celldeath <- celldeath[rowSums(is.na(celldeath))!=ncol(celldeath),]
celldeath[is.na(celldeath)]<-""

write.table(file="./results/celldeath.txt", celldeath, sep="\t", quote=F, col.names=T, row.names=F)



#########################mild SLE immunoinfiltration######################################
source("./ref/geoCRG.CIBERSORT.R")       

outTab=CIBERSORT("./ref/ref.txt", "./GEOexp.txt", perm=1000, QN=T)

conNum <- as.numeric(table(group$group)[1])
treatNum <- as.numeric(table(group$group)[2])

outTab1 <- outTab %>% as.data.frame() %>% rownames_to_column("GSM") %>%
  left_join(group,by="GSM") %>% column_to_rownames("GSM")
table(outTab1$group) #过滤前样本数

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
outTab2 <- outTab1[outTab1[,"P-value"]<0.05,]
table(outTab2$group) #过滤后样本数

outTab3=as.matrix(outTab2[,1:(ncol(outTab2)-4)])
outTab3.out=rbind(id=colnames(outTab3),outTab3)
write.table(outTab3.out, file="./results/CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


data <- t(outTab3)
rt <- outTab3


Type <- outTab2$group
data1=cbind(as.data.frame(t(data)), Type)
data1=melt(data1, id.vars=c("Type"))
colnames(data1)=c("Type", "Immune", "Expression")

#输出小提琴图
outTab2=data.frame()
pdf(file="./results/Figure4. mild.immune.vioplot.pdf", height=8, width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，正常组用蓝色表示，肿瘤组用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  conData=rt[1:conNum,i]
  treatData=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'forestgreen')
  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'orange')
  wilcoxTest=wilcox.test(conData,treatData)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab2=rbind(outTab2,cellPvalue)
  }
  mx=max(c(conData,treatData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topleft", 
       c("Control", "Disease"),
       lwd=3,bty="n",cex=1,
       col=c("forestgreen","orange"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

library(dplyr)

mean = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),mean)
sd = aggregate(x=data1$Expression, by=list(data1$Type,data1$Immune),sd)
immuneDiff = left_join(mean,sd,by=c("Group.1","Group.2"))

immuneDiff <- immuneDiff[which(immuneDiff$Group.2 %in% outTab2$Cell),] %>%
  group_by(Group.2) %>%
  dplyr::select(Cell=Group.2, group=Group.1,mean=x.x,sd=x.y) %>%
  left_join(outTab2,by="Cell")
write.table(immuneDiff,file="./results/immuneDiff.xls",sep="\t",row.names=F,quote=F)

pdf(file="./results/Figure5. barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"normal",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"patient",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

