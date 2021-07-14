library(SingleR)
library(Seurat)
library(scater)

trimmed_mean <- function(row) {
	row = row[order(row)]
	ind = length(row) * 0.1
	mn = mean(row[(ind+1):(length(row)-ind)])
	return(mn)
}


prep_data_to_score <- function(obj, cell_categories, meta_othercol, genes_to_use) {
	objtest = as.matrix(obj@assays$RNA@counts[rownames(obj@assays$RNA@counts) %in% genes_to_use,])
	dim(objtest)
	tcolD = obj@meta.data[,c(cell_categories, meta_othercol)]
	objtest <- SummarizedExperiment(list(counts=objtest), colData=tcolD)
	objtest <- scater::logNormCounts(objtest)

	return(objtest)
}


prep_data_for_training <- function(counts, obj, cell_categories, meta_othercol) {
	Idents(obj)<-obj@meta.data[,cell_categories]

	tcell.types = unique(obj@meta.data[,cell_categories])

	sigs = data.frame(matrix(nrow=nrow(counts), ncol=1))
	rownames(sigs) <- rownames(counts)
	colnames(sigs) <- 'gene'
	sigs$gene <- rownames(sigs)

	show('making SummarizedExperiment object')
	tcolD = obj@meta.data[,c(cell_categories, meta_othercol)]
	sumobj <- SummarizedExperiment(list(counts=counts), colData=tcolD)

	sumobj <- scater::logNormCounts(sumobj)
	logcounts <- sumobj@assays@data$logcounts
	show('normalized data, making lists')
	## interate through the cell types
	for (cell in tcell.types) {
		print(cell)

		expsub <- subset(obj, idents=cell)
		logcountssub <- logcounts[,colnames(expsub)]

		show('applying trimmed means')
		trms <- data.frame('trms'=apply(logcountssub, 1, trimmed_mean))
		colnames(trms)<-c(cell)

		sigs <- cbind(sigs, trms)
	}
	return(sigs)
}


genesintersect = intersect(rownames(test), rownames(train)); length(genesintersect)

testexp = test@assays$RNA@counts
testexp = testexp[genesintersect,]
drop = which(apply(testexp, 1, max)<1)
testexp = testexp[-drop,]; dim(testexp) ## 15225 genes left now

trainexp = train@assays$RNA@counts
trainexp = trainexp[genesintersect,]
drop = which(apply(trainexp, 1, max)<1)
trainexp = trainexp[-drop,]; dim(trainexp) ## 13608 genes left now

finalgenes = intersect(rownames(testexp), rownames(trainexp)); length(finalgenes)
## 13365 genes total to use now

trainexp = trainexp[finalgenes, ]

### train on train
trainsigs <- prep_data_for_training(trainexp, train, train_variable, train_meta_1)
# write.table(trainsigs, 'trainsigs_cd4_training_202008902.txt', sep='\t')
trainsigs$gene <- NULL
tcell.types <- colnames(trainsigs)
trained.object = trainSingleR(ref=as.matrix(trainsigs), labels=tcell.types)

### train labels to our object
test.object <- prep_data_to_score(test, test_meta_1, test_meta_2, finalgenes)

test.scored = classifySingleR(test=test.object, trained=trained.object, fine.tune=TRUE)
colnames(test.scored) <- paste0(assigned.label.name, colnames(test.scored))

test.assigned.labels <- test.scored[,paste0(assigned.label.name, c('labels', 'pruned.labels'))]
superf









