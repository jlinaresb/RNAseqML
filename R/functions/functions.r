# Download RNASEQ V2 data from TCGA firehose
get.rnaseq.data = function(disease){

  require(XML)
  require(gdata)
  require(RCurl)

  dataset = disease
  url = 'http://gdac.broadinstitute.org/runs/stddata__latest/'

  ldoc = htmlTreeParse(getURL(url), useInternalNodes = T)
  datasets = xpathSApply(ldoc, "//a[contains(@href, 'Standardized+Data+Run+Release+Notes')]", xmlValue)
  llinks = unlist(XML::xpathApply(ldoc, "//a[@href]", XML::xmlGetAttr,"href"))
  dlinks = llinks[grepl(paste("/data/", dataset, "/", sep = ""), llinks)]
  ddoc = XML::htmlTreeParse(dlinks, useInternalNodes = T)

  keyword = paste("", "Level_3__RSEM_genes_normalized__data.Level_3", sep = "") # LA PALABRA NORMALIZED MUY IMPORTANTE! HABRÍA QUE MIRAR SI NOS INTERESA BAJARNOS LOS NO NORMALIZADOS
  keyword = paste("//a[contains(@href, '", keyword, "')]", sep = "")
  plinks = XML::xpathSApply(ddoc, keyword, XML::xmlGetAttr, 'href')
  plinks = plinks[grepl(paste("*", dataset, ".Merge_rnaseqv2__illuminahiseq*._rnaseqv2__.*.tar[.]gz$", sep = ""), plinks)]

  timestamp = unlist(strsplit(dlinks, "/"))
  timestamp = timestamp[length(timestamp)]

  download_link = paste(dlinks, trim(plinks[1]), sep = "/")

  utils::download.file(url = download_link, destfile = paste(dataset, "-RNAseq2GeneNorm.tar.gz", sep=""), method = "auto", quiet = T, mode = 'w')

  system(paste('tar xzvf ', dataset, "-RNAseq2GeneNorm.tar.gz", sep = ""))
  grepSearch = paste( dataset, ".rnaseqv2__", sep = "")
  fileList = dir(path = paste(path, 'gdac.broadinstitute.org_', dataset, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/', sep = ''), pattern = grepSearch)

  fname = paste(dataset, "__", timestamp, "-RNAseq2GeneNorm.txt", sep = "")
  setwd(paste(path, 'gdac.broadinstitute.org_', dataset, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/', sep = ''))
  file.rename(from = fileList, to = fname)
  delFodler = paste('~/projects/r-package-tcga-methodology/01_datos', "/", strsplit(fileList, "/")[[1]][1], sep = "")
  #unlink(delFodler, recursive = T)

  tmpCols = utils::read.delim(fname, nrows = 1, colClasses = "character")
  tmpdat = utils::read.delim(fname, skip = 1, sep ="\t")

  colOrder = 1:ncol(tmpCols)
  colOrder = colOrder[tmpCols[1,] == "normalized_count"]

  gene.id = as.vector(as.character(tmpdat$gene_id))

  gnames = sapply(gene.id, function(s) unlist(strsplit(s, "\\|"))[1])
  badg = which(gnames == "?" | duplicated(gnames))

  gdat = tmpdat[-badg, colOrder]
  colnames(gdat) = colnames(tmpCols)[colOrder]
  colnames(gdat) = gsub("\\.", "-", colnames(gdat))
  rownames(gdat) = gnames[-badg]

  gdat = as.data.frame(t(gdat))

  setwd('~/')

  return(gdat)
}



# Sample Selection
sample.selection = function(data, samples = '01'){

  t = sapply(rownames(data), function(s) unlist(strsplit(s,"-"))[4])
  sample.select = data[grep(samples, t),]

  return(sample.select)
}





# Create data
create.data = function(path.brca, path.luad){

  require(dplyr)
  
  print('Cargando datasets ...')
  brca = readRDS(path.brca)
  luad = readRDS(path.luad)

  luad = sample.selection(luad, '01')
  brca = sample.selection(brca, '01')

  luad = cbind(luad, target = rep('luad', nrow(luad)))
  brca = cbind(brca, target = rep('brca', nrow(brca)))
  brca = sample_n(brca, nrow(luad))

  print('Uniendo datasets ...')
  data = rbind(luad, brca)
  data = sample_n(data, nrow(data))

  return(data)

}



# Univariate Filtering with mlr package
fs.abs = function(data, fs.type, nfeat){

  stopifnot('target' %in% names(data))

  require(mlr)
  task = makeClassifTask(data = data, target = 'target')

  tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))

  for (i in 1:length(nfeat)) {
    tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
  }

  t = list()
  for (i in 1:length(tasks)) {
    t[[i]] = tasks[[i]]$env$data
    names(t)[[i]] = paste(fs.type, nfeat[i], sep = '_')
  }

  return(t)
}



# Variable Importance
# El nmodel se refiere al número del modelo. 
# Hay que primero saber lo que se tiene en el objeto bmr para hacer ejecutar la función
varImp.glmnet = function(bmr, n.model){  
  
  models = getBMRModels(bmr)
  glmnet = models[[n.model]]$classif.glmnet.tuned
  
  sum.betas = rep(0)
  for (model in 1:length(glmnet)) {
    
    learner.models = getLearnerModel(glmnet[[model]])
    sum.betas = sum.betas + as.vector((learner.models$learner.model$beta))
  }
  
  names(sum.betas) = glmnet[[n.model]]$features
  barplot(abs(sum.betas), main = 'Glmnet Variable Importance', xaxt = "n", ylim = c(0,50), col = topo.colors(length(sum.betas)))
  
  return(sum.betas)
}


varImp.rf = function(bmr, n.model){
  
  models = getBMRModels(bmr)
  rf = models[[n.model]]$classif.randomForest.tuned
  
  sum.imp = rep(0)
  for (i in 1:length(rf)) {
    
    iters = getFeatureImportance(rf[[i]])
    res = as.matrix(iters$res)
    sum.imp =+ sum.imp + res
  }
  print(length(sum.imp))
  barplot(sum.imp, main = 'Random Forest Variable Importance', xaxt = "n", col = topo.colors(length(sum.imp)))
  return(sum.imp)
  
}

