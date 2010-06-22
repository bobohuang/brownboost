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
    #cat("Working on class: ", class, "\n")
    className <- as.character(class)
    d <- binaryDataFrame(trainingData, class)
    multiEnsemble[[className]] <- bbBuildEnsemble(d, c)
  }
  return (multiEnsemble)
}


bbMultiResult <- function (exampleSet, multiEnsemble) {
  results <- data.frame(rows=1:length(exampleSet[,1]))
  cs <- length(multiEnsemble)
  for (c in 1:cs) {
    results <- cbind(results, class=(bbRunEnsemble(multiEnsemble[[c]], exampleSet)))
  }
  return(results[,-1])
}


bbMultiLeaveOneOut <- function (data, c) {
  resultVector <- c()                              # prediction for each example .. the one left out
  for (i in 1:nrow(data)) {
    example <- data[i,]
    newdata <- data[-i, ]
    multiEnsemble <- bbBuildMultiEnsemble(newdata, c)
    results <- bbMultiResult(example, multiEnsemble)
    bestResults <- apply(X=results, FUN=max, MARGIN=1)
    classIndex <- sapply(1:length(bestResults),
                         function(i, x, y) which(x[i] == y[i,]),
                         x = bestResults, y = results)
    classNames <- names(multiEnsemble)[classIndex]
    cat(example$Class, "\t", classNames, "\n")
    resultVector <- c(resultVector, classNames[1])
  }
  return(resultVector)
}

bbTreeCount <- function(multiEnsemble) {
  numTrees <- c()
  for (m in multiEnsemble) {
    numTrees <- c(numTrees, length(m[[1]]))
  }
  return(numTrees)
}

bbMultiRun <- function (trainset, testset, c) {
  multiEnsemble <- bbBuildMultiEnsemble(trainset, c)
  numTrees <- bbTreeCount(multiEnsemble)
  results <- bbMultiResult(testset, multiEnsemble)
  bestResults <- apply(X=results, FUN=max, MARGIN=1)
  classIndex <- sapply(1:length(bestResults),
                      function(i, x, y) which(x[i] == y[i,]),
                      x = bestResults, y = results)
  classNames <- names(multiEnsemble)[classIndex]
  cat("Trees used: ", numTrees, "\n")
  cat("True Class: ", testset$Class, "\n")
  cat("Pred Class: ", classNames, "\n")
  return(as.numeric(classNames))
}


