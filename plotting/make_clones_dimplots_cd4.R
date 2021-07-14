library(Seurat)
library(tidyverse)
library(openxlsx)

# t.met <- tils@meta.data[,c("seurat_clusters", "final.clonotype.family")]
# t.met$cell.barcode <- rownames(t.met)
# 
# met <- read.xlsx('../TIL-TCR antigen specificity.xlsx')
# dim(met)
# length(intersect(met$TCR.clonotype.family, t.met$final.clonotype.family))
# 
# met <- merge(t.met, met, by.x='final.clonotype.family',by.y='TCR.clonotype.family'); dim(t.met); dim(met)
# 
# met$final.clonotype.family <- NULL
# met$patient <- NULL; met$seurat_clusters <- NULL
# rownames(met) <- met$cell.barcode; met$cell.barcode <- NULL
# write.table(met, 'tspec_clones_metadata_for_plotting.txt', sep='\t')

setwd('~/Dropbox (Partners HealthCare)/Melanoma_CD4/')

tils <- readRDS('Data_Objects/tils.CD4.harmonized.20200514.rds')
tcrs <- read.xlsx('CD4 TCRs.xlsx')
## MAKE THE TABLE TO USE FOR PLOTTING
met <- tils@reductions$umap@cell.embeddings
meta <- tils@meta.data[,c('seurat_clusters','patient','final.clonotype.family')]
meta <- merge(meta, met, by=0); dim(meta)

length(unique(tcrs$TCR.ID))
length(intersect(unique(tcrs$TCR.ID), tils$final.clonotype.family))

tab <- merge(tcrs, meta, by.x='TCR.ID', by.y='final.clonotype.family', all=T); dim(tils); dim(tab)
# tab$Color <- ifelse(is.na(tab$Color), 'lightgrey', as.character(tab$Color))

#### PREP THE DATASET FOR PLOTTING
tab$neocol <- ifelse(tab$Category=='NeoAg' & !is.na(tab$Category), tab$Color, 'lightgrey')
tab$neocat <- ifelse(tab$Category=='NeoAg' & !is.na(tab$Category), tab$Antigen, 'other')
tab$neo.order <- ifelse(tab$neocat=='other', 1, sample(2:3000))

tab$viralcol <- ifelse(tab$Category=='Viral' & !is.na(tab$Category), tab$Color, 'lightgrey')
tab$viralcat <- ifelse(tab$Category=='Viral' & !is.na(tab$Category), tab$Antigen, 'other')
tab$viral.order <- ifelse(tab$viralcat=='other', 1, sample(2:3000))

tab$maacol <- ifelse(tab$Category=='MAA' & !is.na(tab$Category), tab$Color, 'lightgrey')
tab$maacat <- ifelse(tab$Category=='MAA' & !is.na(tab$Category), tab$Antigen, 'other')
tab$maa.order <- ifelse(tab$maacat=='other', 1, sample(2:3000))

tab$Tumorviralcol <- ifelse(tab$Category=='Tumorviral' & !is.na(tab$Category), tab$Color, 'lightgrey')
tab$Tumorviralcat <- ifelse(tab$Category=='Tumorviral' & !is.na(tab$Category), tab$Antigen, 'other')
tab$Tumorviral.order <- ifelse(tab$Tumorviralcat=='other', 1, sample(2:3000))


############ ACTUALLY DO THE PLOTTING
## first have to order the cells so the NA's are first, and the rest are random
## had to manually go in and assign the Colors. if they are in the wrong order for you, just re-order the Colors in the order they appear on the plot

## You can adjust the point size with the 'size' parameter.
# tab[,c("Antigen", "Color")] %>% unique()

### IF YOU WANT TO DOUBLE CHECK THAT THE COLORS ARE CORRECT, THEN COMMENT OUT THE 'NoLegend' 
###   BY CHANGEING '+NoLegend()' to '#+NoLegend()'

tab <- tab[order(tab$neo.order), ]
ggplot(tab) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=neocat), size=1.5) + theme_classic() + 
  scale_colour_manual(values=c('dodgerblue','cyan','seagreen1','steelblue4','lightseagreen',
                               'green','blue','grey86','purple2','darkblue','slateblue1',
                               'lightskyblue','darkcyan','forestgreen')) #+NoLegend()

tab <- tab[order(tab$viral.order), ]
ggplot(tab) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=viralcat), size=1.5) + theme_classic() + 
  scale_colour_manual(values=c('grey90','black','grey86'))#+NoLegend()


tab <- tab[order(tab$maa.order), ]
ggplot(tab) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=maacat), size=1.5) + theme_classic() + 
  scale_colour_manual(values=c('red','grey86','darkorange'))#+NoLegend()


tab <- tab[order(tab$Tumorviral.order), ]
ggplot(tab) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Tumorviralcat), size=1.5) + theme_classic() + 
  scale_colour_manual(values=c('grey90','black', 'grey86'))#+NoLegend()


############# BELOW IS FOR PLOTTING THE BACKGROUND AND DOTS SEPARATELY

### background plots
tab <- tab[order(tab$neo.order), ]
ggplot(tab[tab$neocol=='lightgrey',]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=neocat), size=0.2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('grey86'))+NoLegend()


tab <- tab[order(tab$viral.order), ]
ggplot(tab[tab$viralcol=='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=viralcat), size=0.2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('grey86'))+NoLegend()


tab <- tab[order(tab$maa.order), ]
ggplot(tab[tab$maacol=='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=maacat), size=0.2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('grey86'))+NoLegend()


tab <- tab[order(tab$Tumorviral.order), ]
ggplot(tab[tab$Tumorviralcol=='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Tumorviralcat), size=0.2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('grey86'))+NoLegend()


### overlay plots
tabplot <- tab[!is.na(tab$Antigen), ]; dim(tabplot)

tabplot <- tabplot[order(tabplot$neo.order), ]
ggplot(tabplot[tabplot$neocol != 'lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=neocat), size=2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('dodgerblue','cyan','seagreen1','steelblue4','lightseagreen',
                               'green','blue','purple2','darkblue','slateblue1',
                               'lightskyblue','darkcyan','forestgreen')) + NoLegend()

tabplot <- tabplot[order(tabplot$viral.order), ]
ggplot(tabplot[tabplot$viralcol!='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=viralcat), size=2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5)  +
  scale_colour_manual(values=c('grey60','black'))+NoLegend()


tabplot <- tabplot[order(tabplot$maa.order), ]
ggplot(tabplot[tabplot$maacol!='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=maacat), size=2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('red','darkorange'))+NoLegend()


tabplot <- tabplot[order(tabplot$Tumorviral.order), ]
ggplot(tabplot[tabplot$Tumorviralcol!='lightgrey', ]) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=Tumorviralcat), size=2) + theme_classic() + xlim(-7,8.5) + ylim(-6, 6.5) +
  scale_colour_manual(values=c('grey60','black'))+NoLegend()






