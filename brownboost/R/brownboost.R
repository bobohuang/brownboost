library(RWeka)

# This is going to return a collection of weights and classifiers.

# dataFrame contains the true vector in a named column TrueClass
# c is a positive valued parameter
# v >0 is a small constant used to avoid degenerate cases

brownboost <- function (trainingData, c, v) {

  l <- length(trainingData[,1])

  # The margin
  r <- rep(0, l)
  s <- c
  i <- 1
  while (s > 0) {

    # each example gets a positive weight
    weights <- exp(-(r + s)^2 / c)
    sampProb <- weights / sum(weights)
    sampledData <- sample(seq(0, l), replace=F, prob=sampProb, size=0.25*l)
    classifier <- J48(formula=TrueClass ~ ., data=sampledData)
