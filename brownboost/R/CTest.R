dyn.load("brownboost/src/bisection.so")
source("brownboost/R/BisectionDESolver.R")
rm(.Random.seed)

desolverTest1 <- function(r, s, h, y, c) {
  try (
    result <- .C("solvede", as.double(r), as.double(s), as.double(h),
                 as.double(y), as.double(c), as.integer(length(r)), output=numeric(2))
       )
    return(result$output)
}

#void bigfun (double* a, double* b, #
#	     double* v, double* x,
#	     double* c,
#	     int* sign, int n)
bigfunTest <- function(a, b, v, z, c) {
  result <- .C("bigfun", as.double(a), as.double(b), as.double(-1.0),
               as.double(z), as.double(c), output=integer(2), as.integer(length(b)))
  return(result$output)
}


# This is just testing that the error functions
# between R and the C code agree.... and they do.

test1 <-  function() {

  print("starting test...")
  rm(.Random.seed)
  
  for ( i in 1:10) {
    
    c <- 3                                      # c can be defined from epsilon
    s <- 3                                        # step
    r <- rep(0.0, times=500)
    h <- sample(c(-1, 1), 500, replace=T)
    y <- h
    ydex <- sample(1:500, 50, replace=T)
    y[ydex] <- y[ydex] * -1
    a <- r + s;
    b <- h * y;
    v <- list(b, -1)
    z <- sample(seq(-1, 1, by=0.0001), replace=T, size=2)
    
    r_result <- f(a, b, v, z, c)
    c_result <- bigfunTest(a, b, v, z, c)
    print(r_result)
    print(sign(r_result))
    print(c_result)
  }
}


# This is the main test of the de solver.

test2 <- function() {
  c <- 4                                        # c can be defined from epsilon
  s <- c                                        # step
  r <- rep(0, times=500)                        # This is the margin r(x)
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  y <- h                                        # This is the true class
  sampsize <- sample(1:500, replace=F, size=1)
  i = sample(x=seq(1, 500), replace=T, size=sampsize)
  y[i] <- -1 * y[i]                             # Give 30% error to y(x)
  result <- desolverTest1(r, s, h, y, c)
  print(result)
}
