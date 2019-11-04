# Packages
require(mlr)
require(tweeDEseq)

source('~/RNAseqML/R/functions/functions.r')
load('~/RNAseqML/data/filter_genes.RData')

print('Load and match datasets ...')
data = create.data(path.brca = '~/RNAseqML/data/brca.rds', path.luad = '~/RNAseqML/data/luad.rds')
target = as.character(data$target)
rnames = colnames(data)

i = intersect(colnames(data), genes)
xdata = data[,i]

# Normalize by TMM approach
data_tmm = normalizeCounts(xdata, method = 'TMM')

xdata = as.data.frame(cbind(data_tmm, target))
xdata$target = as.character(xdata$target)
names(xdata) = make.names(names(xdata))
xdata[] <- lapply(xdata, type.convert, as.is = TRUE)

print('Running Feature Selection. Method: Kruskal Test ...')
# Feature Selection (Filter univariate)
task = makeClassifTask(data = xdata, target = 'target')

# nfeat = 20
nfeat = c(10, 20, 40, 80)
fs.type = 'kruskal.test'
tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))

for (i in 1:length(nfeat)) {
  tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
}

tdata = list()
for (i in 1:length(tasks)) {
  tdata[[i]] = tasks[[i]]$env$data
  names(tdata)[[i]] = paste(fs.type, nfeat[i], sep = '_')
}

# Machine Learning (Random Forest and Glmnet)
source('~/RNAseqML/R/functions/machineLearning.r')
execute.ml(list.data = tdata, path = '~/git/RNAseqML/results/', filename = 'ML-Benchmark-result.rds', win = F) # si se corre en windows hay que poner win = T para la paralelización


bmr = readRDS('~/RNAseqML/results/ML-Benchmark-result.rds')

# Plot the results
png(filename = '~/RNAseqML/plots/BMRsummary.png')
plotBMRSummary(bmr)
dev.off()

png(filename = '~/RNAseqML/plots/BMRboxplot.png')
plotBMRBoxplots(bmr, style = 'violin')
dev.off()

png(filename = '~/RNAseqML/plots/varImp-glmnet.png')
varImp.glmnet(bmr, n.model = 2)
dev.off()

png(filename = '~/RNAseqML/plots/varImp-rf.png')
varImp.rf(bmr, n.model = 2)
dev.off()

# Analisis de los resultados
# bmr = readRDS('~/tmp/XoveTIC/XoveTIC_ml_bmr.rds')
# bmr
# boxplot = plotBMRBoxplots(bmr)
# boxplot + theme(axis.text.x =  element_blank(),
#                 axis.ticks.x = element_blank())
#
# plotBMRSummary(bmr)
#
# source('~/r-package-tcga-methodology/5.Anal_Results/results_analysis.r')
# vi = varImp.glmnet(bmr, 1)
#
# d = data.frame(Importance = vi)
# genes = rownames(d)
# ggplot2::ggplot(d, aes(x = genes, y = Importance)) + geom_bar(stat = 'identity', fill = 'steelblue') + theme(axis.title.x = element_blank())
