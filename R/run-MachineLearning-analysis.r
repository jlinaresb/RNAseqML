# Packages
require(mlr)
require(tweeDEseq)

source('~/git/RNAseqML/R/functions/functions.r')
load('~/git/RNAseqML/data/filter_genes.RData')

print('Load and match datasets ...')
data = create.data(path.brca = '~/git/RNAseqML/data/brca.rds', path.luad = '~/git/RNAseqML/data/luad.rds')
target = as.character(data$target)
rnames = colnames(data)

i = intersect(colnames(data), genes)
xdata = data[,i]

# Normalize by TMM approach
data_tmm = normalizeCounts(xdata, method = 'TMM')

xdata = as.data.frame(cbind(data_tmm, target))
xdata$target = as.character(xdata$target)
names(xdata) = make.names(names(xdata))

print('Running Feature Selection. Method: Kruskal Test ...')
# Feature Selection (Filter univariate)
task = makeClassifTask(data = xdata, target = 'target')

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
source('~/git/RNAseqML/R/functions/machineLearning.r')
execute.ml(list.data = tdata, path = '~/git/RNAseqML/results/', filename = 'ML-Benchmark-result.rds')

# Analisis de los resultados
# bmr = readRDS('~/tmp/XoveTIC/XoveTIC_ml_bmr.rds')
# bmr
# boxplot = plotBMRBoxplots(bmr)
# boxplot + theme(axis.text.x =  element_blank(),
#                 axis.ticks.x = element_blank())
#
# plotBMRSummary(bmr)
#
# source('~/git/r-package-tcga-methodology/5.Anal_Results/results_analysis.r')
# vi = varImp.glmnet(bmr, 1)
#
# d = data.frame(Importance = vi)
# genes = rownames(d)
# ggplot2::ggplot(d, aes(x = genes, y = Importance)) + geom_bar(stat = 'identity', fill = 'steelblue') + theme(axis.title.x = element_blank())
