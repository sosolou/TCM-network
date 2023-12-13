rm(list = ls()) #清除缓存

setwd("C:/Users/Administrator/Desktop/SLE/procedure/step11. DGIdb")
library(statnet)
library(circlize)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

gene <- read.table("C:/Users/Administrator/Desktop/SLE procedures/step6. Hubgene/results/biomarkers.txt", header=T, sep="\t", check.names=F)
interaction <- read.table("interactions.csv", header=T, sep=",", check.names=F) %>%
  filter(Name %in% gene$id) %>%
  arrange(Name,desc(`Interaction Score`)) %>%
  group_by(Name) %>%
  filter(row_number()<=4) %>%
  select(-`Interaction Type & Directionality`) %>%
  spread(key = "Name",
         value = `Interaction Score`)  %>% 
  arrange(Drug) %>%
  column_to_rownames("Drug") 

# 按照所有列排序
interaction<- interaction[do.call(order, lapply(interaction,function(x)-x)),] %>%
  as.matrix()



nm <- unique(unlist(dimnames(interaction)))
group <- c(rep("Drug",nrow(interaction)),rep("Gene",ncol(interaction)))
names(group) <- nm

grid.col = NULL
grid.col[rownames(interaction)] = brewer.pal(nrow(interaction),"Set3") 
grid.col[colnames(interaction)] = brewer.pal(nrow(interaction),"Set2") 


pdf(file="Figure 1.mild_gene_drug.circ.pdf", width=10, height=10)
chordDiagram(interaction, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col,group=group)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .2, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.8)
  circos.axis(h = "top", labels.cex = 0.6, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


