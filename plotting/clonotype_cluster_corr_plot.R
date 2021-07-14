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


setwd('~/Dropbox (Partners HealthCare)/Melanoma_CD4/')
tils <- readRDS('Data_Objects/tils.CD4.harmonized.20200514.rds')

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





### disregard below, I think

################################################################################################################################
################################################################################################################################
################################################################################################################################

cluster.table <- clones[,c('clone', as.character(0:10))]
## find duplicates in the clone names
duprows <- cluster.table$clone[duplicated(cluster.table$clone)]
duped <- cluster.table[cluster.table$clone %in% duprows,]
cluster.table <- cluster.table[setdiff(rownames(cluster.table), c('4413','3984')),]

rownames(cluster.table) <- cluster.table$clone
cluster.table[is.na(cluster.table[])] <- 0
cluster.table$clone <- NULL
colnames(cluster.table) <- paste0('cluster_', colnames(cluster.table))

ct.sub <- cluster.table[1:10,]

cluster.table.p <- data.frame(t(apply(cluster.table, 1, percent)))


require(ggpubr)

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

tables <- DoCorr(cluster.table.p, 'kendall')
pearson.t <- DoCorr(cluster.table.p, 'pearson')
spearman.t <- DoCorr(cluster.table.p, 'spearman')

pearson.og <- DoCorr(cluster.table, 'pearson')
pearson.og2 <- cor(cluster.table, method='pearson')

ps <- log10(0.99999999999999-tables[['p.table']])
ps = scale(tables[['p.table']])

percent <- function(row) {
  rsum <- sum(row)
  perc <- row / rsum * 100
  return(perc)
}

spear <- cor(cluster.table.p, method='spearman')
corrplot(spear, method="circle")
corrplot(as.matrix(spearman.t$rho.table), method='circle')

################################################################################################################################
############################################ GRAPHING PEARSON CORRELATION ######################################################
################################################################################################################################

pearson.og <- DoCorr(cluster.table, 'pearson')
pearson.og2 <- cor(cluster.table, method='pearson')
sqorder <- paste0('cluster_', c(12,2,1,7,10,3,6,9,4,5,8,11,0))

require(corrplot)
corrplot(pearson.og2, is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 0.6, tl.col = 'black')
corrplot(pearson.og2[sqorder, sqorder], is.corr = FALSE, method =c('circle'), addgrid.col = 'grey', tl.cex = 1, tl.col = 'black')

################################################################################################################################
################################################################################################################################
################################################################################################################################

plot_size <- as.matrix(read.table('plot_size.txt', sep = '\t', header = TRUE))

plot_significance <- as.matrix(read.table('plot_significance.txt', sep = '\t', header = TRUE))

fdr_template <- as.matrix(read.table('fdr_template.txt', sep = '\t', header = TRUE))

fdr_significance <- as.matrix(read.table('fdr_significance.txt', sep = '\t', header = TRUE))

## Define colors
col <- colorRampPalette(c(rep('dodgerblue', 5), 'white', rep('firebrick', 5)))

require(corrplot)
##Create significance plot
corrplot(as.matrix(1-tables[['p.table']]), method =c('circle'), col=('black'), tl.cex = 0.6, tl.col = ('black'), cl.pos='n')
corrplot(as.matrix(1-pearson.t[['p.table']]), method =c('circle'), col=('black'), tl.cex = 0.6, tl.col = ('black'), cl.pos='n')
corrplot(as.matrix(1-spearman.t[['p.table']]), method =c('circle'), col=('black'), tl.cex = 0.6, tl.col = ('black'), cl.pos='n')
## Create size measure plot
corrplot(as.matrix(tables[['rho.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')
corrplot(as.matrix(tables[['stat.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')
corrplot(as.matrix(pearson.t[['rho.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 0.6, tl.col = 'black', cl.pos='n')
corrplot(as.matrix(spearman.t[['rho.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')
##Create FDR plot	
corrplot(fdr_template, method =c('circle'), col=('black'), tl.cex = 0.6, tl.col = ('black'), p.mat = fdr_significance, insig = 'blank', sig.level = 0.009, cl.pos='n')

## Create legends	
corrplot(legend_significance, method =c('circle'), col=('grey'), tl.cex = 0.6, tl.col = ('black'), cl.pos='n')
corrplot(plot_size, is.corr = FALSE, method =c('color'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='r', cl.lim = c(-3,3), cl.ratio=0.4, cl.length=7)



corrplot(as.matrix(tables[['rho.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')
corrplot(as.matrix(tables[['rho.table']]), method =c('circle'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')
corrplot(as.matrix(tables[['rho.table']]), method =c('square'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')

tbl <- tables[['rho.table']]
tbl[tbl[]==1e0] <- 0.5
tbl$cluster_0[tbl$cluster_0==1.0] <- 0.5
corrplot(as.matrix(tbl), method =c('square'), addgrid.col = 'grey', col = col(200), tl.cex = 0.6, tl.col = 'black', cl.pos='n')






for.g <- tables[['rho.table']]


write.table(for.g, 'Tables/kendall_corr_values.txt', sep='\t')

f.p <- read.delim('Tables/kendall_corr_values.txt', sep='\t')
##12-2-1-7-10 (little space/white line) -3-6-9 (little space/white line) 5-4-8-11-0
sqorder <- paste0('cluster_', c(12,2,1,7,10,3,6,9,4,5,8,11,0))
f.p <- f.p[sqorder,sqorder]

f.p[f.p[]==1.0] <- 0.25
f.p[] <- f.p[] * 4

pal <- colorRampPalette(c(rep('dodgerblue', 5), 'white', rep('firebrick', 5)))

corrplot(as.matrix(f.p), method =c('circle'), addgrid.col = 'grey', tl.cex = 0.8, tl.col = 'black',
         title='Cluster Correlations')#, p.mat=as.matrix(tables$p.table))

corrplot(as.matrix(f.p), method="circle", addgrid.col='grey', tl.cex=0.8, tl.col='black')

pt <- pearson.t[['rho.table']]
corrplot(as.matrix(pearson.t[['rho.table']]), is.corr = FALSE, method =c('color'), addgrid.col = 'grey', tl.cex = 0.6, tl.col = 'black', cl.pos='n')


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

head(cluster.table.p)
count.table <- cluster.table.p

res.table <- data.table::data.table()
sd.table <- data.table::data.table()
for (col in colnames(count.table)) {
  colvals <- c()
  sdvals <- c()
  for (col2 in colnames(count.table)) {
    numnonzero <- nrow(count.table[count.table[,col] > 0, ])
    numbothnonzero <- nrow(count.table[count.table[,col] > 0 & count.table[,col2] > 0, ])
    # cat(col, col2, 'counts:', numnonzero, numbothnonzero, '\n\n')
    appeartogether <- numbothnonzero / numnonzero
    background <- RandomizeTable(count.table[,col], count.table[,col2])
    numsd <- (appeartogether  - mean(background)) / sd(background)
    colvals <- c(colvals, appeartogether)
    sdvals <- c(sdvals, numsd)
  }
  res.table[,col] <- colvals
  sd.table[,col] <- sdvals
}
res.table <- data.frame(res.table)
rownames(res.table) <- colnames(res.table)


sqorder <- paste0('cluster_', c(12,2,1,7,10,3,6,9,4,5,8,11,0))
res.table <- res.table[sqorder, sqorder]
corrplot(as.matrix(res.table), method =c('color'), addgrid.col = 'grey', tl.cex = 0.6, tl.col = 'black')

require(pheatmap)
sd.table <- data.frame(sd.table)
rownames(sd.table) <- colnames(sd.table)
sd.table[is.na(sd.table[])] <- 0
sd.table[sd.table[] == 'Inf'] <- 100 
sd.table[sd.table[] > -min(sd.table)] <- -min(sd.table)
ppal <- colorRampPalette(colors=c( 'dodgerblue', 'white', rep('firebrick3', 6)))
ppal2 <- colorRampPalette(colors = c('dodgerblue', 'white', 'firebrick3'))
pheatmap(sd.table[sqorder,sqorder], cluster_rows=F, cluster_cols=F, color=ppal2(75))

mod.table <- res.table
mod.table <- mod.table - 0.5
corrplot(as.matrix(mod.table), method =c('color'), addgrid.col = 'grey', tl.cex = 0.6, tl.col = 'black')

res.scaled <- scale(res.table)



RandomizeTable <- function(c1, c2) {
  backgroundscores <- c()
  for (i in 1:10) {
    r1 <- sample(c1)
    r2 <- sample(c2)
    df <- data.frame('c1'=r1, 'c2'=r2)
    numnonzero <- nrow(df[df$c1 > 0, ])
    numbothnonzero <- nrow(df[df$c1 > 0 & df$c2 > 0, ])
    appeartogether <- numbothnonzero / numnonzero
    backgroundscores <- c(backgroundscores, appeartogether)
    c1 <- r1
    c2 <- r2
  }
  return(backgroundscores)
}

numnonzero <- nrow(count.table[count.table$cluster_5 > 0, ])
numbothnonzero <- nrow(count.table[count.table$cluster_5 > 0 & count.table$cluster_4 > 0, ])
# cat(col, col2, 'counts:', numnonzero, numbothnonzero, '\n\n')
appeartogether <- numbothnonzero / numnonzero; appeartogether
back <- RandomizeTable(count.table$cluster_5, count.table$cluster_4); back
sd(back)

numsd <- (appeartogether  - mean(back)) / sd(back)

