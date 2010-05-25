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
  multiEnsemble <- list()                       # list of ensembles, named by class
  for (class in allClasses) {
    cat("Working on class: ", class, "\n")
    className <- as.character(class)
    d <- binaryDataFrame(trainingData, class)
    multiEnsemble[[className]] <- bbBuildEnsemble(d, c)
  }
  return (multiEnsemble)
}


bbMultiResult <- function (example, multiEnsemble) {
  results <- c()
  for (classifier in multiEnsemble) {
    results <- c(results, bbRunEnsemble(classifier, example))
  }
  print("EXAMPLE")
  print(example$Class)
  print(results)
  print(names(multiEnsemble)[which(max(results) == results)])
  bestResult = names(multiEnsemble)[which(max(results) == results)]
                                        # index of the maximum
  return(bestResult)
}


bbMultiLeaveOneOut <- function (data, c) {
  resultVector <- c()                              # prediction for each example .. the one left out
  for (i in 1:nrow(data)) {
    example <- data[i,]
    newdata <- data[-i, ]
    multiEnsemble <- bbBuildMultiEnsemble(newdata, c)
    bestResult <- bbMultiResult(example, multiEnsemble)  # index to ensemble
    resultVector <- c(resultVector, bestResult)
    cat(example$Class, "\t", names(multiEnsemble[[bestResult]]), "\n")
  }
  return(resultVector)
}

bbMultiTestSet <- function (trainset, testset, c) {
  results <- c()
  multiEnsemble <- bbBuildMultiEnsemble(trainset, c)
  print(names(multiEnsemble))
  for (i in 1:nrow(testset)) {
    results <- c(results, bbMultiResult(testset[i,], multiEnsemble))
  }
  return(results)
}



