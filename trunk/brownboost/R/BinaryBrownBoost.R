library(RWeka)
source("brownboost/R/SolverWrapper.R", local=T)

# This is going to return a collection of weights and classifiers.

# dataFrame contains the true class vector in a named column Class
#    which is in values of -1, 1 ... numeric OK!
# c is a positive valued parameter
# v > 0 is a small constant used to avoid degenerate cases .. not used..

bbBuildEnsemble <- function (trainingData, c) {
  print("building binary ensemble")
  rows <- length(trainingData[,1])      # number of examples in the data
  y    <- trainingData$Class            # the true class
  browns <- list()                      # list of classifiers
  length(browns) <- 100                 # start with space for 100 elements...
  alphas <- list()                      # list of weights
  length(alphas) <- 100
  r <- rep(0, rows)                     # The margin
  s <- c                                # The step
  i <- 1                                # The iteration
  index <- 1
  listlen <- 100
  
  while (s > 0) {
                                        # each example gets a positive weight
    weights <- exp(-(r + s)^2 / c)
    sampProb <- weights / sum(weights)
                                        # sample from the data to train the classifier
    sampledData <- sample(seq(1, rows), replace=T, prob=sampProb, size=0.25*rows)
    classifier <- DecisionStump(Class ~ ., data=trainingData, subset=sampledData)
                                        # run the classifier on all the training data to get h(x)
    h <- predict(classifier, newdata=trainingData)
                                        # solve for alpha and t
    alphaAndTee <- solvede(r, s, h, y, c)
    alpha <- alphaAndTee[1]
    tee <- alphaAndTee[2]
                                        #update the margin
    r <- r + alpha * h * y
                                        #update the time remaining
    s <- s - tee
    print(s)
                                        #store the classifier and alpha for the ensemble
    alphas[[index]] <- alpha
    browns[[index]] <- classifier
    index <- index + 1
                                        # need more memory!
    if (index >= listlen) {
      listlen <- listlen + 100
      length(alphas) <- listlen
      length(browns) <- listlen
    }    
  }

  # remove any extra null elements in the list
  browns <- browns[!sapply(browns, is.null)]
  alphas <- alphas[!sapply(alphas, is.null)]
  return(list(ensemble=browns, alphas=alphas))
}



# This function takes the ensemble and test data,
# and resurns a vector of predictions for each
# example.  
bbRunEnsemble <- function (ensemble, data) {

  alphas <- ensemble[[2]]
  browns <- ensemble[[1]]
  results <- rep(0, length(data$Class))
                                        # list of predictions, as a vector, for each classifier
  xs <- lapply(browns, function(x) predict(x, newdata=data))

                                        # multiply the alpha to each classifier ..
  ys <- lapply(1:length(alphas), function(i, x, y) x[[i]] * y[[i]],
               x = xs, y = alphas)
                                          # sum across vectors to get the prediction for each example
  for (y in ys) {
    results <- results + y
  }

  return( sign(results) )
}

  

bbConfusionMatrix <- function (results, classes) {

  print(results)
  print(classes)
  
  confMatrix <- matrix(data=c(0,0,0,0), nrow=2)

  for (i in 1:length(classes)) {
    c <- convertToRIndex(results[i])
    r <- convertToRIndex(classes[i])
    confMatrix[r,c] <- confMatrix[r,c] + 1
  }

  return(confMatrix)
}


bbAccuracy <- function (results, classes) {
  total <- length(results)
  wrong <- length(which(results * classes == -1))
  return((total-wrong) / total)
}


convertToRIndex <- function(x) {
  if (x == 1) {
    return(1)
  } else if (x == -1) {
    return(2)
  } else {
    return(0)
  }
}

