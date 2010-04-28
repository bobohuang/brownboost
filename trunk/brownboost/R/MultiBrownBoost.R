source("brownboost/R/BinaryBrownBoost.R")


#multi boost!

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
  allClasses <- unique(trainingData$Class)      # classes in data
  resultFrame <- data.frame()
  for (class in allClasses) {
    className <- as.character(class)
    d <- binaryDataFrame(trainingData, class)
    ensemble <- multiEnsemble[[as.character(className)]]
    resultFrame <- cbind(resultFrame, className=runBinaryEnsemble(ensemble, d))
  }
  return(resultFrame)
}


