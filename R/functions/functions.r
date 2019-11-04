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
  barplot(abs(sum.betas), main = 'Glmnet Variable Importance', xaxt = "n", ylim = c(0,10), col = topo.colors(length(sum.betas)))
  
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

