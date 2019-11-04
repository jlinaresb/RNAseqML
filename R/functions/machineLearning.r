# Machine Learning with mlr package

execute.ml = function(list.datasets, path = '', filename = '', win = F){

  require(mlr)
  print('Making task')
  l = list()
  for (i in 1:length(list.datasets)) {
    l[[i]] = makeClassifTask(id = paste('nfeat_', ncol(list.datasets[[i]])-1), data = list.datasets[[i]], target = 'target')
  }
  print('Removing Constant Features')
  f = list()
  for (i in 1:length(l)) {
    f[[i]] = removeConstantFeatures(l[[i]])
  }
  print('Normalizing Features')
  n = list()
  for (i in 1:length(f)) {
    n[[i]] = normalizeFeatures(f[[i]])
  }


  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()  #Hacer que se pueda modificar la búsqueda!!!!
  inner<-makeResampleDesc("Holdout")

  # GLMNET
  psglmnet = makeParamSet(
    makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
    makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1))
  )
  l<-makeLearner("classif.glmnet", predict.type = "prob")
  lrn_glmnet<-makeTuneWrapper(l, inner, psglmnet, measures = auc, ctrl, show.info=T)

  # Random Forest
  psrf<-makeParamSet(
    makeDiscreteParam("mtry", values = sqrt(ncol(n[[i]]$env$data))),
    makeDiscreteParam("ntree", values= 1000L),
    makeDiscreteParam("nodesize", values= c(1:3))
  )
  l<-makeLearner("classif.randomForest", predict.type = "prob")
  lrn_rf<-makeTuneWrapper(l,  resampling = inner, par.set = psrf, measures = auc, control=ctrl,  show.info = T)


  learners = list(lrn_glmnet, lrn_rf)

  # Outer Cross-Validation
  # outer = rep(list(makeResampleDesc('RepCV' , reps = 5, folds = 10 , stratify = T)), length(list.datasets))
  outer = rep(list(makeResampleDesc('CV' , iters = 3,  stratify = T)), length(list.datasets)) 

  print('Training the model')


  if (win == TRUE){
    require(parallelMap)
    parallelStartSocket(2, level = 'mlr.tuneParams')
  } else{
    library(parallelMap)
    parallelStartMulticore(2L , level = 'mlr.tuneParams')
  }


  # Benchmarking
  bmr = benchmark(learners, n , outer , measures =  list(acc , auc, mmce) , show.info = T , models = T)
  
  saveRDS(bmr, file = paste0(path, filename, sep = ''))

  parallelStop()

}
