source("brownboost/R/BinaryBrownBoost.R", local=T)

runBinaryBrownBoost <- function (trainingData, testData, c) {
  
  ensemble <- bbBuildEnsemble(trainingData, c)
  results  <- bbRunEnsemble(ensemble, testData)
  cat("Number of Stumps used: ", length(ensemble[[1]]), "\n")
  cat("Accuracy: ", bbAccuracy(results, testData$Class), "\n")
  print("Confusion Matrix - rows are true classes, cols are predictions")
  print(bbConfusionMatrix(results, testData$Class))

}


randomBinaryBrownBoost <- function (data, c, v=NULL) {
  # v is the sample to sample from
  if(is.null(v)) {
    print("creating sample vector")
    v <- 1:length(data[,1])
  }
  sampleSize <- 0.75 * length(v)
  trainDex <- sample(v, replace=F, size=sampleSize)
  print(trainDex)
  testDex <- which(!v %in% trainDex)
  print(testDex)
  runBinaryBrownBoost(data[trainDex,], data[testDex,], c)
}
