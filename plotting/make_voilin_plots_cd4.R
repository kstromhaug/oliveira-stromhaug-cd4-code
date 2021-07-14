library(Seurat)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(Matrix)


setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cd4_analysis/')

tils <- readRDS('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/Seurat_Objects/tils.CD4.harmonized.20200514.rds')


Idents(tils) <- tils$seurat_clusters

VlnPlot(tils, split.by='seurat_clusters', features=c('LEF1','CCR7','IL7R'), pt.size=0)

cluster_order = c('4','0','2','9','5','1','7','8','3','6')
vgenes <- c('LEF1','CCR7','IL7R','TCF7','GZMM','GNLY','KLRB1','HSPA1B','NR4A1','FOXP3','IL2RA','CTLA4','MKI67','TOX','PDCD1','HAVCR2','ENTPD1')
### make table for ggplot

setdiff(vgenes, rownames(tils))
genexp <- data.frame(t(tils@assays$RNA@data[vgenes, ]))

genexp <- merge(tils@meta.data[,c('seurat_clusters','patient')], genexp, by=0); dim(genexp)

genexp.m <- melt(genexp, id.vars=c('Row.names','seurat_clusters','patient'))
genexp.m <- genexp.m[genexp.m$seurat_clusters != 10, ]
genexp.m$cluster <- factor(genexp.m$seurat_clusters, levels = cluster_order)
gcolors <- c('0'= 'dodgerblue','1'='magenta','2'='green3','3'='darkgoldenrod1','4'='blue','5'='turquoise3','6'='red','7'='darkmagenta','8'='darkorange1', '9'='blueviolet')
## epanechnickov, biweight, cosine
## nrd, ucv, (bcv), JS-ste, SJ-dpi
ggplot(genexp.m) + geom_violin(aes(x=variable, y=value, fill=seurat_clusters), color='gray60', scale='width', kernel='gaussian', adjust=1) + 
  scale_y_continuous(name="Expression", limits=c(-2, 10), breaks=c(0,5,10)) +
  theme(axis.ticks = element_blank()) +
  facet_grid(rows='cluster') + 
  theme_classic() + scale_fill_manual(values=gcolors) +
  theme(axis.text.x = element_text(angle = 90)) + xlab('Gene') + NoLegend() + 
  theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_blank())

pdf("figures/tils_cd4_violins.pdf", width=5, height=5)
# insert ggplot code
dev.off()


#################
liu.sigs <- readRDS('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/additional_analysis/liu_et_al/liu_et_al_all_sigs.rds')

sigs_list = liu.sigs[c('Ttol_NR4A1_reg_up', 'Ttol_rel_up')]

tils <- AddModuleScore(tils, features=sigs_list, name='Liu')

liu.mod <- tils@meta.data[,paste0('Liu',1:2)]
colnames(liu.mod) <- names(sigs_list)
write.table(liu.mod, 'liu_modulescores_on_cd4.txt', sep='\t')

tils.meta <- tils@meta.data[,c("seurat_clusters", "patient")]
tils.meta <- tils.meta[sample(1:nrow(tils.meta)), ]
tils.meta <- tils.meta[order(tils.meta$seurat_clusters), ]
tils.meta <- tils.meta[tils.meta$seurat_clusters!=10, ]


### make violing plots with the liu signatures
liu.mod.u <- liu.mod[,c("Ttol_NR4A1_reg_up", "Ttol_rel_up")]
liu.mod.u <- merge(tils.meta, liu.mod.u, by=0); dim(liu.mod.u)
liu.mod.u <- melt(liu.mod.u, id.vars=c('Row.names','seurat_clusters','patient'))
liu.mod.u$cluster <- factor(liu.mod.u$seurat_clusters, levels = cluster_order)
liu.mod.u$variable <- factor(liu.mod.u$variable, levels = c('Ttol_rel_up','Ttol_NR4A1_reg_up'))
ggplot(liu.mod.u) + geom_violin(aes(x=cluster, y=value, fill=cluster), color='gray60', scale='width', kernel='gaussian', adjust=1) + # , scale='width'
  scale_y_continuous(name="Signature Score", limits=c(-0.25, 1.5), breaks=c(0,0.5,1,1.5)) +
  facet_grid(rows='variable') + 
  theme_classic() + scale_fill_manual(values=gcolors) + NoLegend() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab(label='') +
  theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_blank())

pdf("figures/liu_signatures_violin_cd4.pdf", width=5, height=2.5)
# insert ggplot code
dev.off()  
  
  
  
  
