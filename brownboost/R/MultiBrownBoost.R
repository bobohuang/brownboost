source("brownboost/R/BinaryBrownBoost.R")
library(RWeka)

#multiboost!

binaryDataFrame <- function (data, class) {
  data[which(data$Class != class),]$Class <- -1
  data[which(data$Class == class),]$Class <- 1
  return(data)
}



bbBuildMultiEnsemble <- function (trainingData, c) {
  allClasses <- unique(trainingData$Class)      # classes in data 
  multiEnsemble <- list()                       # list of ensembles

  for (class in allClasses) {
    cat("Working on class: ", class, "\n")
    className <- as.character(class)
    d <- binaryDataFrame(trainingData, class)
    multiEnsemble[[className]] <- bbBuildEnsemble(d, c)
  }

  return (multiEnsemble)
}


bbRunMultiEnsemble <- function (testData, multiEnsemble) {
  allClasses <- unique(testData$Class)      # classes in data
  resultFrame <- data.frame()
  for (class in allClasses) {
    className <- as.character(class)
    d <- binaryDataFrame(testData, class)
    ensemble <- multiEnsemble[[as.character(className)]]
    results  <- bbRunEnsemble(ensemble, d)
    cat("\n\nFor class: ", class, "\n")
    cat("Number of Stumps used: ", length(ensemble[[1]]), "\n")
    cat("Accuracy: ", bbAccuracy(results, d$Class), "\n")
    print("Confusion Matrix - rows are true classes, cols are predictions")
    print(bbConfusionMatrix(results, d$Class))
  }
  return(resultFrame)
}


runMultiBrownBoost <- function (trainingData, testData, c) {  
  multiensemble <- bbBuildMultiEnsemble(trainingData, c)
  bbRunEnsemble(testData, multiensemble)
}
