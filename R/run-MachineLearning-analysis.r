# Packages
require(mlr)
require(tweeDEseq)

source('functions/functions.r')
load('../data/filter_genes.RData')

print('Load and match datasets ...')
data = create.data(path.brca = '../data/brca.rds', path.luad = '../data/luad.rds')
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
source('functions/machineLearning.r')
bmr=execute.ml(list.data = tdata, win = F) # si se corre en windows hay que poner win = T para la paralelización

# Plot the results
png(filename = '../plots/BMRsummary.png')
plotBMRSummary(bmr)
dev.off()

png(filename = '../plots/BMRboxplot.png')
plotBMRBoxplots(bmr, style = 'violin')
dev.off()

png(filename = '../plots/varImp-glmnet.png')
varImp.glmnet(bmr, n.model = 2)
dev.off()

png(filename = '../plots/varImp-rf.png')
varImp.rf(bmr, n.model = 2)
dev.off()
