library(Seurat)
library(tidyverse)
library(openxlsx)

percent <- function(row) {
  rsum <- sum(row)
  perc <- row / rsum * 100
  return(perc)
}

DoCorr <- function(data, method) {
  rho.table <- data.table::data.table()
  p.table <- data.table::data.table()
  stat.table <- data.table::data.table()
  for (col1 in colnames(data)) {
    rowp <- c()
    rowrho <- c()
    rowstat <- c()
    show(col1)
    for (col2 in colnames(data)) {
      corr <- cor.test(data[,col1], data[,col2], method = method)
      rowp <- c(rowp, corr$p.value)
      rowrho <- c(rowrho, corr$estimate)
      rowstat <- c(rowstat, corr$statistic)
    }
    rho.table[,col1] <- rowrho
    p.table[,col1] <- rowp
    stat.table[,col1] <- rowstat
  }
  rho.table <- data.frame(rho.table); rownames(rho.table) <- colnames(rho.table)
  p.table <- data.frame(p.table); rownames(p.table) <- colnames(p.table)
  stat.table <- data.frame(stat.table); rownames(stat.table) <- colnames(stat.table)
  p.log <- -log10(p.table)
  p.rev <- 1-p.table
  
  return(list('rho.table'=rho.table, 'p.table'=p.table, 'stat.table'=stat.table))
}

MakeExcel <- function(corrtable, ptable, ptrans, filename) {
  wb <- createWorkbook()
  name <- "Correlation"
  addWorksheet(wb, name)
  writeData(wb, sheet = name, corrtable)
  
  name <- "Pvalues"
  addWorksheet(wb, name)
  writeData(wb, sheet = name, ptable)
  
  name <- "P_transformed"
  addWorksheet(wb, name)
  writeData(wb, sheet = name, ptrans)
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}


setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/')
tils <- readRDS('Seurat_Objects/tils.CD4.harmonized.20200514.rds')

clons <- read.xlsx('../CD4TIL clones.xlsx')
clones <- clons %>% separate(Clonotype.ID, '-->', into=c('clone', 'matched.clone'), remove=F)
clones <- clones %>% separate(clone, '-', into=c('patient', 'clone.1', 'clone.2'), remove=F)
clones$clone.2[is.na(clones$clone.2)] <- ""
clones$clonotype.family.name <- ifelse(clones$clone.2=="", clones$clone.1, paste(clones$clone.1, clones$clone.2, sep='-'))

clones.4 <- subset(clones, Sum >= 4); nrow(clones); nrow(clones.4)

require(ggpubr)
require(corrplot)


clorder <- 1:10
cn <- paste0('cluster_', clorder)

############ ON EVERYTHING ##############
clones.4 <- subset(clones, Sum >= 4); nrow(clones); nrow(clones.4)
c4.forcor <- clones.4[,as.character(c(0:10, 'clone'))]
rownames(c4.forcor) <- c4.forcor$clone; c4.forcor$clone <- NULL
colnames(c4.forcor) <- paste0('cluster_', colnames(c4.forcor))
c4.forcor[is.na(c4.forcor[])] <- 0

c4.docorr <- DoCorr(c4.forcor, 'spearman')
ptable <- (1 - as.matrix(c4.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/spearman_raw_n4_tables.xlsx')

c4.p <- data.frame(t(apply(c4.forcor, 1, percent)))
c4.p.docorr <- DoCorr(c4.p, 'spearman')
ptable <- (1 - as.matrix(c4.p.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.p.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.p.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/spearman_percent_n4_tables.xlsx')

############### ON PATIENT 15 ################
p15.clones.4 <- clones.4[clones.4$patient=='p15',]; dim(clones.4); dim(p15.clones.4)
p15.c4.forcor <- p15.clones.4[,as.character(c(0:10, 'clone'))]
rownames(p15.c4.forcor) <- p15.c4.forcor$clone; p15.c4.forcor$clone <- NULL
colnames(p15.c4.forcor) <- paste0('cluster_', colnames(p15.c4.forcor))
p15.c4.forcor[is.na(p15.c4.forcor[])] <- 0

c4.docorr <- DoCorr(p15.c4.forcor, 'spearman')
ptable <- (1 - as.matrix(c4.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p15_spearman_raw_n4_tables.xlsx')

c4.p <- data.frame(t(apply(p15.c4.forcor, 1, percent)))
c4.p.docorr <- DoCorr(c4.p, 'spearman')
ptable <- (1 - as.matrix(c4.p.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.p.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.p.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p15_spearman_percent_n4_tables.xlsx')

############### ON PATIENT 2 ################
p2.clones.4 <- clones.4[clones.4$patient=='p2',]; dim(clones.4); dim(p2.clones.4)
p2.c4.forcor <- p2.clones.4[,as.character(c(0:10, 'clone'))]
rownames(p2.c4.forcor) <- p2.c4.forcor$clone; p2.c4.forcor$clone <- NULL
colnames(p2.c4.forcor) <- paste0('cluster_', colnames(p2.c4.forcor))
p2.c4.forcor[is.na(p2.c4.forcor[])] <- 0

c4.docorr <- DoCorr(p2.c4.forcor, 'spearman')
ptable <- (1 - as.matrix(c4.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p2_spearman_raw_n4_tables.xlsx')

c4.p <- data.frame(t(apply(p2.c4.forcor, 1, percent)))
c4.p.docorr <- DoCorr(c4.p, 'spearman')
ptable <- (1 - as.matrix(c4.p.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.p.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.p.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p2_spearman_percent_n4_tables.xlsx')

############### ON PATIENT 11-pre ################
p11.clones.4 <- clones.4[clones.4$patient=='p11pre',]; dim(clones.4); dim(p11.clones.4)
p11.c4.forcor <- p11.clones.4[,as.character(c(0:10, 'clone'))]
rownames(p11.c4.forcor) <- p11.c4.forcor$clone; p11.c4.forcor$clone <- NULL
colnames(p11.c4.forcor) <- paste0('cluster_', colnames(p11.c4.forcor))
p11.c4.forcor[is.na(p11.c4.forcor[])] <- 0

c4.docorr <- DoCorr(p11.c4.forcor, 'spearman')
ptable <- (1 - as.matrix(c4.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p11_spearman_raw_n4_tables.xlsx')

c4.p <- data.frame(t(apply(p11.c4.forcor, 1, percent)))
c4.p.docorr <- DoCorr(c4.p, 'spearman')
ptable <- (1 - as.matrix(c4.p.docorr[['p.table']])^(1/3))[clorder, clorder]
corrtable <- as.matrix(c4.p.docorr[['rho.table']])[clorder, clorder]
corrplot(corrtable, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper')
corrplot(ptable, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black', type='upper', col='black', cl.pos='n')
MakeExcel(corrtable, as.matrix(c4.p.docorr[['p.table']])[clorder, clorder], ptable, 'cd4_analysis/p11_spearman_percent_n4_tables.xlsx')

