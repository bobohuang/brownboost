source("brownboost/R/BinaryBrownBoost.R", local=T)

runBinaryBrownBoost <- function (trainingData, testData, c) {

  print("testdata")
  print(testData$Class)
  print("traindata")
  print(trainingData$Class)
  
  ensemble <- bbBuildEnsemble(trainingData, c)
  results  <- bbRunEnsemble(ensemble, testData)

  cat("Number of Stumps used: ", length(ensemble[[1]]), "\n")
  cat("Accuracy: ", bbAccuracy(results, testData$Class), "\n")
  print("Confusion Matrix - rows are true classes, cols are predictions")
  print(bbConfusionMatrix(results, testData$Class))

}


randomBinaryBrownBoost <- function (data, c) {
  v <- 1:length(data$Class)
  sampleSize <- 0.80 * length(v)
  trainDex <- sample(x=v, replace=F, size=sampleSize)
  testDex <- which(!v %in% trainDex)
  runBinaryBrownBoost(data[trainDex,], data[testDex,], c)
}
