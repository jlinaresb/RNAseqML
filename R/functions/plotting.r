# Carga de los datos
require(tweeDEseq)
require(edgeR)
require(vegan)
require(ggplot2)
require(ggbiplot)

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


# Visualization
############################

# Perfil genético medio de la población por cada gen
brca.mean = apply(data_tmm[which(target == 'brca'),], 2, mean)
luad.mean = apply(data_tmm[which(target == 'luad'),], 2, mean)

png(filename = '~/RNAseqML/plots/expression_profile.png')
par(mfrow = c(2,1))
barplot(brca.mean, main = 'BRCA population', axisnames = F)
barplot(luad.mean, main = 'LUAD population', axisnames = F)
dev.off()


#MAplot
png(filename = '~/RNAseqML/plots/MA-plot.png')
maPlot(data_tmm[52,], data_tmm[2,],
       pch = 19, cex = .5, ylim = c(-8,8),
       allCol = "darkgrey", lowess = T,
       xlab = "Average expression",
       ylab = "Mean fold change")
grid(col="black")
title("MA-plot RNAseq v2 TMM")
dev.off()


#MDS plot
mds = metaMDS(comm = data_tmm, distance = 'bray')
mds_xy = data.frame(mds$points)
png(filename = '~/RNAseqML/plots/metaMDSplot.png')
ggplot(mds_xy, aes(MDS1, MDS2, color = target)) +geom_point() + theme_bw()
dev.off()


# PCA plot
pca = prcomp(data_tmm, center = T, scale. = T)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
              groups = target, ellipse = TRUE,
              circle = TRUE, var.axes = F, alpha = 0.5)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
png(filename = '~/RNAseqML/plots/PCAplot.png')
g
dev.off()