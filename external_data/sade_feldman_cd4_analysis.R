library(tidyverse)
library(Seurat)

setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/Sade-Feldman/')
source('~/Documents/source_code/seurat_functions.R')

sfmeta <- read.delim('sf_clusters.txt', sep='\t', stringsAsFactors = F); dim(sfmeta)
rownames(sfmeta) <- sfmeta$Cell.Name


sf = readRDS('s_f_seurat_unprocessed.rds')
sf <- AddMito(sf)
sf <- AddMetaData(sf, sfmeta)

Idents(sf) <- sf$patient
VlnPlot(sf, features=c('nCount_RNA'), pt.size = 0)
VlnPlot(sf, features=c('nFeature_RNA'), pt.size = 0)
VlnPlot(sf, features=c('percent_mito'),pt.size = 0)

sf <- DoSeurat(sf)
DimPlot(sf, group.by='Cluster.number')
DimPlot(sf, group.by = 'seurat_clusters')

FeaturePlot(sf, c('CD8A','CD3E','CD3D','CD4'))
FeaturePlot(sf, 'CD8A', max.cutoff = 0.25)

removecd8 <- colnames(sf@assays$RNA@scale.data[, sf@assays$RNA@scale.data['CD8A', ] > 0.25])

cd4 <- subset(sf, idents=)


sf2 <- read.delim('GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt', sep='\t')
sf2 <- sf2[2:nrow(sf2), ]
colnames(sf2) <- colnames(sf)
expcd8 <- t(sf2['CD8A',])
expcd8 <- as.numeric(expcd8)
names(expcd8) <- colnames(sf2)
remcd8 <- expcd8[expcd8 == 0]; length(remcd8)

notcd8 <- subset(sf, cells=names(remcd8))

FeaturePlot(notcd8, features=c('CD8A','CD4','CD3E','CD3D'))
FeaturePlot(sf, features=c('CD8A','CD4','CD3E','CD3D'))

umap = data.frame(notcd8@reductions$umap@cell.embeddings)
umap = umap[umap$UMAP_1 > -7, ]; dim(umap)

notcd8 <- subset(notcd8, cells=rownames(umap))
FeaturePlot(notcd8, features=c('CD8A','CD4','CD3E','CD3D'))


notcd8 <- DoSeurat(notcd8)
DimPlot(notcd8)
FeaturePlot(notcd8, features=c('CD8A','CD4','CD3E','CD3D'))

umap = data.frame(notcd8@reductions$umap@cell.embeddings)
umap = umap[umap$UMAP_1 > -5, ]; dim(umap)

cd4 <- subset(notcd8, cells=rownames(umap))

cd4 <- DoSeurat(cd4)
DimPlot(cd4)
FeaturePlot(cd4, features=c('CD8A','CD4','CD3E','CD3D'))

counts <- cd4@assays$RNA@counts
counts <- counts[, counts['CD8A', ] > 0]

cd4 <- subset(cd4, cells = setdiff(colnames(cd4), colnames(counts)))

cd4 <- DoSeurat(cd4)
cd4 <- AddMetaData(cd4, sfmeta)
DimPlot(cd4, label=T)
DimPlot(cd4, group.by='Cluster.number')
FeaturePlot(cd4, features=c('CD8A','CD4','CD3E','CD3D'))
FeaturePlot(cd4, features=l_immune)
View()

FeaturePlot(cd4, features=c('CD8A','CD8B','CD4','CD3D','CD3E','CD69','MS4A1','HLA-DRB1','ITGAM','GNLY','CD19'))

FeaturePlot(cd4, features=c('ITGAM','CD3D'))

smoothScatter()


########################################################################################################################
########################################################################################################################


sf = readRDS('s_f_seurat_unprocessed.rds')

sfmeta <- read.delim('sf_clusters.txt', sep='\t', stringsAsFactors = F); dim(sfmeta)
rownames(sfmeta) <- sfmeta$Cell.Name
sfmeta <- sfmeta[colnames(sf), ]
sfmeta$Cluster.number <- paste0('G', sfmeta$Cluster.number)
keep <- sfmeta[sfmeta$Cluster.number %in% c('G5','G6','G7','G8','G9','G10','G11'),]; dim(keep)

sf <- AddMito(sf)
sf <- AddMetaData(sf, sfmeta)
table(sf$Cluster.number)

sf <- subset(sf, cells=rownames(keep))
table(sf$Cluster.number)
sf <- DoSeurat(sf)

DimPlot(sf, label=T)
DimPlot(sf, group.by='Cluster.number', label=T)

par(mfrow=c(1,2))
# plot.new()
smoothScatter(sf@assays$RNA@data["CD19", ], sf@assays$RNA@data["CD3D", ],nrpoints=Inf, xlab = "CD19", ylab = "CD3D", xlim = c(-2, 4), ylim = c(0, 5))
smoothScatter(sf@assays$RNA@data["CD19", ], sf@assays$RNA@data["CD3E", ],nrpoints=Inf, xlab = "CD19", ylab = "CD3E", xlim = c(-2, 4), ylim = c(0, 5))

par(mfrow=c(1,2))
smoothScatter(sf@assays$RNA@data["CD8A", ], sf@assays$RNA@data["CD3E", ],nrpoints=Inf, xlab = "CD8A", ylab = "CD3E", xlim = c(0, 4), ylim = c(0, 5))
smoothScatter(sf@assays$RNA@data["CD8B", ], sf@assays$RNA@data["CD3E", ],nrpoints=Inf, xlab = "CD8B", ylab = "CD3E", xlim = c(0, 4), ylim = c(0, 5))

par(mfrow=c(1,2))
smoothScatter(sf@assays$RNA@data["CD3E", ], sf@assays$RNA@data["CD4", ],nrpoints=Inf, xlab = "CD3E", ylab = "CD4", xlim = c(0, 4), ylim = c(0, 5))
smoothScatter(sf@assays$RNA@data["CD8A", ], sf@assays$RNA@data["CD4", ],nrpoints=Inf, xlab = "CD8A", ylab = "CD4", xlim = c(0, 4), ylim = c(0, 5))


FeaturePlot(sf, features=c('CD8A','CD8B','CD4','CD3D','CD3E','CD69','MS4A1','HLA-DRB1','ITGAM','GNLY','CD19'))

scale <- sf@assays$RNA@scale.data
counts <- sf@assays$RNA@counts
data <- sf@assays$RNA@data

cd8.remove <- data[, data['CD8A',] > 0.5 | data['CD8B',] > 0.5]; dim(cd8.remove)
cd3.remove <- data[, data['CD3E',] + data['CD3D',] == 0]; dim(cd3.remove)
cd19.remove <- data[, data['CD19',] > 0 | data['ITGAM',] > 0]; dim(cd19.remove)
remove <- unique(c(colnames(cd8.remove), colnames(cd3.remove), colnames(cd19.remove))); length(remove)
keep <- setdiff(colnames(sf), remove); length(keep)

sf2 <- subset(sf, cells = setdiff(colnames(sf), remove)); dim(sf2)

sf2 <- DoSeurat(sf2, res=0.45)

DimPlot(sf2, group.by='seurat_clusters', label=T)
DimPlot(sf2, group.by='Cluster.number')
FeaturePlot(sf2, features=c('CD4','CD3D','CD3E','CD69','MS4A1','HLA-DRB1','ITGAM','GNLY','CD19'))
FeaturePlot(sf2, features=c('CD4','CD3D','CD3E'))

sf2 <- subset(sf2, idents=6, invert=T)
sf2 <- DoSeurat(sf2, res=0.45)

DimPlot(sf2, group.by='seurat_clusters', label=T)
DimPlot(sf2, group.by='Cluster.number')
# FeaturePlot(sf2, features=c('CD4','CD3D','CD3E','CD69','MS4A1','HLA-DRB1','ITGAM','GNLY','CD19'))
FeaturePlot(sf2, features=c('CD4','CD3D','CD3E','GNLY'))


cd4.0 <- sf2@assays$RNA@data[,sf2@assays$RNA@data['CD4',] == 0]; dim(cd4.0)
cd4.l <- sf2@assays$RNA@data[,sf2@assays$RNA@data['CD4',] > 0]; dim(cd4.l)



sf2 <- FindClusters(sf2, resolution = 0.5)
sf2$clusters_0.5 <- Idents(sf2)
p5 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.55)
sf2$clusters_0.55 <- Idents(sf2)
p55 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.6)
sf2$clusters_0.6 <- Idents(sf2)
p6 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.7)
sf2$clusters_0.7 <- Idents(sf2)
p7 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.8)
sf2$clusters_0.8 <- Idents(sf2)
p8 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.85)
sf2$clusters_0.85 <- Idents(sf2)
p85 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.9)
sf2$clusters_0.9 <- Idents(sf2)
p9 <- DimPlot(sf2, label=T)
sf2 <- FindClusters(sf2, resolution = 0.95)
sf2$clusters_0.95 <- Idents(sf2)
p95 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.0)
sf2$clusters_1.0 <- Idents(sf2)
p10 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.05)
sf2$clusters_1.05 <- Idents(sf2)
p105 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.1)
sf2$clusters_1.1 <- Idents(sf2)
p11 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.2)
sf2$clusters_1.2 <- Idents(sf2)
p12 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.3)
sf2$clusters_1.3 <- Idents(sf2)
p13 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.4)
sf2$clusters_1.4 <- Idents(sf2)
p14 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.5)
sf2$clusters_1.5 <- Idents(sf2)
p15 <- DimPlot(sf2, label=T)

sf2 <- FindClusters(sf2, resolution = 1.6)
sf2$clusters_1.6 <- Idents(sf2)
p16 <- DimPlot(sf2, label=T)
# 
figure <- ggarrange(p3, p35,p4,p45,p5,p55,p6,p7,p8,
                    labels = c('0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.7','0.8'),
                    ncol = 3, nrow = 3)
figure

figure <- ggarrange(p7,p8,p85,p9,p95,p105,
                    labels = c('0.7','0.8','0.85','0.9','0.95','1.05'),
                    ncol = 2, nrow = 3)
figure

ggarrange(p13,p14,p15,p16,f1,f2,
                    labels = c('1.3','1.4','1.5','1.6','IL2RA','FOXP3'),
                    ncol = 3, nrow = 2)

FeaturePlot(sf2, features=c('IL2RA','FOXP3'))

f1 <-FeaturePlot(sf2, 'IL2RA')
f2 <- FeaturePlot(sf2, 'FOXP3')

DimPlot(sf2, group.by='patient')

ggarrange(p85,p9,p95,p105,p11,p12,p13,f1, f2,
          labels=c('0.85','0.9','0.95','1.05','1.1','1.2','1.3', 'IL2RA','FOXP3'),
          ncol=3, nrow=3)
saveRDS(sf2, 'cd4/sf_cd4_clustered_resolutions.rds')

meta <- sf2@meta.data

for.g <- meta[,c('patient', 'clusters_1.0', "Cluster.number")]
for.g$cell.barcode <- rownames(for.g)

write.table(for.g, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/Sade-Feldman/cd4/cd4_cells_clusters.txt', sep='\t', row.names = F)


sf <- readRDS('sf_cd4_clustered_resolutions.rds')
sf_markers <- FindAllMarkers(sf, logfc.threshold=0.1, min.pct=0.05)
sf_markers_reg <- FindAllMarkers(sf)

write.table(sf_markers, 'sf_cd4_markers_inclusive.txt', sep='\t', row.names=F)
write.table(sf_markers_reg, 'sf_cd4_markers_normal.txt', sep='\t', row.names=F)



tils.labels <- read.delim('sf_cd4_with_tilscd4_singler_labels.txt', sep='\t')

sf <- AddMetaData(sf, metadata = tils.labels)
DimPlot(sf, group.by='oliveira.labels')

tils.table <- sf@meta.data[,c('seurat_clusters','oliveira.labels')]
tils.table$seurat_clusters <- paste0('sf_', tils.table$seurat_clusters)
tils.table$oliveira.labels <- paste0('oliveira_', tils.table$oliveira.labels)
tils.table <- table(tils.table)
tils.table <- data.frame(tils.table)
tils.table <- dcast(tils.table, seurat_clusters ~ oliveira.labels)
rownames(tils.table) <- tils.table$seurat_clusters
tils.table$seurat_clusters <- NULL

for (i in 1:nrow(tils.table)) {
  tils.table[i,] <- tils.table[i,] / sum(tils.table[i,])
}

for (i in 1:ncol(tils.table)) {
  tils.table[,i] <- tils.table[,i] / sum(tils.table[,i])
}

tils.table$sf_clusters <- rownames(tils.table)
tils.table <- tils.table[,c('sf_clusters',paste0('oliveira_',0:10))]

write.table(tils.table, 'sf_oliveira_cd4_singler_labels_20210125.txt', sep='\t', row.names=F)



tils.with.sf.labels <- read.delim('tils_with_sf_labels_20210125.txt', sep='\t')
tils <- readRDS('../../tils_combined/cd4_analysis/tils_cd4_labeltransfer_20200902.rds')

tils <- AddMetaData(tils, metadata=tils.with.sf.labels)
tils.table <- tils@meta.data[,c('oliveira.labels','seurat_clusters')]
tils.table$seurat_clusters <- paste0('sf_', tils.table$seurat_clusters)
tils.table$oliveira.labels <- paste0('oliveira_', tils.table$oliveira.labels)
tils.table <- table(tils.table)
tils.table <- data.frame(tils.table)
tils.table <- dcast(tils.table, oliveira.labels ~ seurat_clusters)
rownames(tils.table) <- tils.table$oliveira.labels
tils.table$seurat_clusters <- NULL

sf <- readRDS('sf_cd4_clustered_resolutions.rds')
tcrs.r <- read.xlsx('cd4/sf_cd4_clonotype_families_revised.xlsx')
names = read.delim('cd4_cell_with_terra_barcode.txt', sep='\t', stringsAsFactors = F)

cols.to.remove <- c('HLA.A...allele.1','HLA.A...allele.2','HLA.B...allele.1','HLA.B...allele.2','HLA.C...allele.1','HLA.C...allele.2',
                    'TRAV','TRAJ','TRBV','TRBD','TRBJ','alpha.aa','beta.aa','alpha.n','beta.n',
                    'clone.freq','clonotype','clone','found',
                    'CDR3...alpha.or.gamma.chain','CDR3...beta.or.delta.chain',
                    'CDR3..AA....alpha.or.gamma.chain','CDR3..AA....beta.or.delta.chain',
                    'alpha.gamma.V','alpha.gamma.J','beta.delta.V','beta.delta.J','beta.delta.D',
                    'TCRalpha.chain.ID','TCRbeta.chain.ID',
                    'clonotype.family.size','Clonotype.family.number')
sf@meta.data[,colnames(sf@meta.data) %in% cols.to.remove] <- NULL

head(names)
tcrs.add <- merge(tcrs, names[,c('cell.barcode','barcode')]); dim(tcrs); dim(tcrs.add)
length(intersect(tcrs.add$cell.barcode, colnames(sf)))
rownames(tcrs.add) <- tcrs.add$cell.barcode

duped <- tcrs.add[duplicated(tcrs.add$cell.barcode) | duplicated(tcrs.add$cell.barcode, fromLast=T),]
View(duped[,c('barcode','alpha_1','alpha_2','beta_1','beta_2', "p_clonotype_size")])
indices_to_remove <- c(26,83,126,129,650,681,812,848)

notcr <- sf@meta.data[! rownames(sf@meta.data) %in% tcrs.add$cell.barcode, ]; dim(notcr)
notcr$timepoint <- notcr$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.
notcr$characteristics..patinet.ID..Pre.baseline..Post..on.treatment. <- NULL
notcr <- notcr[,c('Row.names','seurat_clusters','patient','timepoint')]
write.table(notcr, 'cells_with_no_tcr.txt', sep='\t', row.names=F)

