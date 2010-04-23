library(RWeka)
source("~/School/aml/brownboost/brownboost/R/deSolver.R")

# This is going to return a collection of weights and classifiers.

# dataFrame contains the true vector in a named column TrueClass
#     which is the last column
# c is a positive valued parameter
# v >0 is a small constant used to avoid degenerate cases

brownboost <- function (trainingData, c, v) {

  rows <- length(trainingData[,1])
  cols <- length(names(trainingData))
  data <- trainingData[,1:cols-1]
  y    <- trainingData[,cols]

  browns <- c()
  alphas <- c()
  
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

    # sample from the data to train the classifier
    sampledData <- sample(seq(0, l), replace=F, prob=sampProb, size=0.25*l)
    classifier <- J48(formula=TrueClass ~ ., data=data[sampledData,])

    # run the classifier on all the training data to get h(x)
    h <- predict(classifier, newdata=trainingData)

    # solve for alpha and t
    alphaAndTee <- solvede(r, s, h, y, c)
    alpha <- aAndt[1]
    tee <- aAndt[2]

    #update the margin
    r <- r + alpha * h * y

    #update the time remaining
    s <- s - tee

    #store the classifier and alpha for the ensemble
    alphas <- c(alphas, alpha)
    browns <- c(browns, classifier)
  }

  return(list(ensemble=browns, alphas=alphas))
}
