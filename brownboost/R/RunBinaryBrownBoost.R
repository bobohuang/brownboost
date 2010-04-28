source("brownboost/R/BinaryBrownBoost.R", local=T)

runBinaryBrownBoost <- function (trainingData, testData, c) {
  
  ensemble <- bbBuildEnsemble(trainingData, c)
  results  <- bbRunEnsemble(ensemble, testData)
  cat("Number of Stumps used: ", length(ensemble[[1]]), "\n")
  cat("Accuracy: ", bbAccuracy(results, testData$Class), "\n")
  print("Confusion Matrix - rows are true classes, cols are predictions")
  print(bbConfusionMatrix(results, testData$Class))

}
