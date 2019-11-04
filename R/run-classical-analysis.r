# Require packages
require(tweeDEseq)
require(edgeR)
require(vegan)
require(ggplot2)
library(ggbiplot)

source('~/git/RNAseqML/R/functions/functions.r')
load('~/git/RNAseqML/data/filter_genes.RData')

data = create.data(path.brca = '~/git/RNAseqML/data/brca.rds', path.luad = '~/git/RNAseqML/data/luad.rds')
target = data$target
rnames = colnames(data)

i = intersect(colnames(data), genes)
xdata = data[,i]

# Normalize by TMM approach
data_tmm = normalizeCounts(xdata, method = 'TMM')

xdata = as.data.frame(cbind(data_tmm, target), row.names = rnames)
xdata$target = as.character(xdata$target)

xdata.t <- as.data.frame(t(xdata), stringsAsFactors = FALSE)
xdata.t[] <- lapply(xdata.t, type.convert, as.is = TRUE)
xdata.t = as.matrix(xdata.t)

# Differential expression analysis with edgeR
d<- DGEList(counts= xdata.t, group=target)
d<- estimateCommonDisp(d)
d<- estimateTagwiseDisp(d)
res.edgeR.common<- exactTest(d, pair=c("brca", "luad"), dispersion="common")
res.edgeR.tagwise<- exactTest(d, pair=c("brca", "luad"), dispersion="tagwise")

par(mfrow = c(2,1))
hist(res.edgeR.common$table$PValue, main = 'Histogram P-values CommonDisp')
abline(v=0.05,col="red")
hist(res.edgeR.tagwise$table$PValue, main = 'Histogram P-values TagwiseDisp')
abline(v=0.05,col="red")

plotBCV(d)
plotMD(res.edgeR.tagwise, main = 'Mean-Difference Plot of Expression Data')

x = topTags(res.edgeR.tagwise, n = 20)
top20_edger = rownames(x$table)

save(top20_edger, file = '~/git/RNAseqML/results/top20_DE.RData')

