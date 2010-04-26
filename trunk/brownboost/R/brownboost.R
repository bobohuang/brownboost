library(RWeka)
source("~/School/aml/brownboost/brownboost/R/deSolver.R")

# This is going to return a collection of weights and classifiers.

# dataFrame contains the true class vector in a named column Class
#    which is in values of -1, 1 ... numeric!
# c is a positive valued parameter
# v > 0 is a small constant used to avoid degenerate cases

brownboost <- function (trainingData, c, v, testData) {

  rows <- length(trainingData[,1])
  cols <- length(names(trainingData))
  y    <- trainingData[,"Class"]
  #trainingData[,"Class"] <- as.factor(trainingData[,"Class"])

  resultsTable <- data.frame(row.names=seq(1, length(testData[,1])))
  weightsTable <- data.frame(row.names=seq(1, rows))
  
  #browns <- c()
  #alphas <- c()
  
  # The margin
  r <- rep(0, rows)

  # The step
  s <- c
  
  # The iteration
  i <- 1
  while (s > 0) {

    # each example gets a positive weight
    weights <- exp(-(r + s)^2 / c)
    sampProb <- weights / sum(weights)
    weightsTable <- cbind(weightsTable, weights)
    
    # sample from the data to train the classifier
    sampledData <- sample(seq(1, rows), replace=T, prob=sampProb, size=0.25*rows)
    classifier <- DecisionStump(Class ~ ., data=trainingData, subset=sampledData)
    #print(summary(classifier))

    # run the classifier on all the training data to get h(x)
    h <- predict(classifier, newdata=trainingData)
    #h <- (as.numeric(as.character(h))) # comes back as a factor

    # solve for alpha and t
    alphaAndTee <- solvede(r, s, h, y, c)
    alpha <- alphaAndTee[1]
    tee <- alphaAndTee[2]

    #update the margin
    r <- r + alpha * h * y

    #update the time remaining
    s <- s - tee
    print(s)

    resultsTable <- cbind(resultsTable, alpha*sign(predict(classifier, newdata=testData)))
    
    #store the classifier and alpha for the ensemble
    #alphas <- c(alphas, alpha)
    #browns <- c(browns, classifier)
  }

  rs <- seq(1, length(weightsTable[1,]))
  plot (x = rs, y = weightsTable[1,])

  
  #return(list(ensemble=browns, alphas=alphas))
  return(resultsTable)
}



brownAccuracy <- function (resultsTable, classes) {

  confMatrix <- matrix(data=c(0,0,0,0), nrow=2)
  l <- length(classes)
  
  resultsSign <- as.vector(sign(apply(FUN=sum, X=t1, MARGIN=c(1))))
  wrong <- which(resultsSign * classes == -1)


  for (i in seq(1,l)) {
    r <- convertToR(classes[i])
    c <- convertToR(resultsSign[i])
    confMatrix[r][c] <- confMatrix[r][c] + 1
  }

  print("Confusion Matrix ... rows = true")
  print(confMatrix)
  
  confMatrix <- matrix(data=c(0,0,0,0), nrow=2)
  for (i in wrong) {
    r <- convertToR(classes[i])
    c <- convertToR(resultsSign[i])
    confMatrix[r][c] <- confMatrix[r][c] + 1
  }
  print("Confusion Matrix of errors ... rows = true")
  print(confMatrix)

  cat("Accuracy = ", length(wrong) / length(classes), "\n")
  return(wrong)
}


convertToR <- function(x) {
  if (x == 1) {
    return(1)
  } else if (x == -1) {
    return(2)
  } else {
    return(0)
  }
}
