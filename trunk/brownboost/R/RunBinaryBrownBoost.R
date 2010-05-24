source("brownboost/R/BinaryBrownBoost.R", local=T)

runBinaryBrownBoost <- function (trainingData, testData, c) {

  #print("testdata")
  #print(testData$Class)
  #print("traindata")
  #print(trainingData$Class)
  
  ensemble <- bbBuildEnsemble(trainingData, c)
  results  <- bbRunEnsemble(ensemble, testData)
  # in this case we need the sign of the results
  results <- sign(results)
  #cat("Number of Stumps used: ", length(ensemble[[1]]), "\n")
  #cat("Accuracy: ", bbAccuracy(results, testData$Class), "\n")
  #print("Confusion Matrix - rows are true classes, cols are predictions")
  #print(bbConfusionMatrix(results, testData$Class))
  cat(bbAccuracy(results, testData$Class), "\t", length(ensemble[[1]]), "\t", c, "\n")
  #return( bbAccuracy(results, testData$Class))
}


randomBinaryBrownBoost <- function (data, c) {
  v <- 1:length(data$Class)
  sampleSize <- 0.80 * length(v)
  trainDex <- sample(x=v, replace=F, size=sampleSize)
  testDex <- which(!v %in% trainDex)
  runBinaryBrownBoost(data[trainDex,], data[testDex,], c)
}


varyCBinaryBrownBoost <- function (data) {
  v <- 1:length(data$Class)
  sampleSize <- 0.80 * length(v)
  varyc <- c()
  trainDex <- sample(x=v, replace=F, size=sampleSize)
  testDex <- which(!v %in% trainDex)

  for (i in seq(from=1, to=8.25, by=0.25)) {
    a <- runBinaryBrownBoost(data[trainDex,], data[testDex,], i)
    varyc <- c(varyc, a)
  }

  return(varyc)
}
