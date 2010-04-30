
# contains the function erf(a)
library(NORMT3)

# dot product of one vector and one list
dotl <- function (v, l) {
    return(v[1]*l[[1]] + v[2]*l[[2]])
}

# This function removes the imaginary part from
# the result of erf .... which is usually 0 in
# my exerience... Also the erf is defined differently
# in the papers ... it's 2/pi instead of 2/sqrt(pi)...

erfd <- function(a) {
  #cat("min a: ", min(a), "\n")
  #return((1/pi)*as.double(erf(a)))
  return(2*pnorm(a*sqrt(2)) - 1)
}

# when the solution is found, this result should be c(0,0)
# returns true if each function returns zero
# v = [b, -1] and z = [a, t]

boundryCondition <- function (a, b, v, z, c) {
  #print("BOUNDRY")
  #cat("a ", a[1], "\n")
  #cat("b ", b[1], "\n")
  #cat("v ", v[[1]][1], "\n")
  #cat("z ", z, "\n")
  #cat("c ", c, "\n")
  #cat("err1 ", (erfd((a + dotl(z, v))/sqrt(c)))[1], "\n")
  #cat("err2 ", (erfd(a/sqrt(c)))[1], "\n")
            
  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  f2 <- sum(erfd((a + dotl(z, v))/sqrt(c)) - erfd(a/sqrt(c)))

  #cat("f1 ", f1, "\n")
  #cat("f2 ", f2, "\n")
  #print("END BOUNDRY")
  return(data.frame(f1=f1,f2=f2))
}

# when calculating the jacobian matrix of the boundry condition
# functions, we need elements W, U, B, and V.
# when calculating the updates for alpha and tee we need E
# v = [b, -1] and z = [alpha, tee]
jacobianElements <- function(a, b, v, z, c) {
  d <- a + dotl(z,v)
  w <- exp(-d^2/c)
  W <- sum(w)
  U <- sum(w*d*b)
  B <- sum(w*b)
  V <- sum(w*d*b*b)
  E <- sum(erfd(d/sqrt(c)) - erfd(a/sqrt(c)))
  return(data.frame(W=W,U=U,B=B,V=V,E=E))
}


updateStep <- function(alpha, tee, c, W, U, B, V, E) {
  if (W != 0) {
    Bk <- B/W
    Vk <- V/W
    Uk <- U/W
    a2 <- alpha + (c*Bk + sqrt(pi*c)*Uk*(E/W)) / (2*(Vk - Uk*Bk))
    t2 <- tee + (c*Bk^2 + sqrt(pi*c)*Vk*(E/W)) / (2*(Vk - Uk*Bk))
  } else {
    a2 = alpha
    t2 = tee
  }
  return(data.frame(alpha=a2, tee=t2))
}

#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
solvede <- function(r, s, h, y, c) {
  alpha <- 0
  tee <- 0
  a <- r + s
  b <- h * y
  z <- c(alpha, tee)
  v <- list(b, -1)
  boundry <- boundryCondition(a, b, v, z, c)
  lastBoundry <- 0
  # run until the variables stop changing...
  while (abs(sum(lastBoundry) - sum(boundry)) > 1e-11) {
    j <- jacobianElements(a, b, v, z, c)
    update <- updateStep(alpha, tee, c, j$W, j$U, j$B, j$V, j$E)
    z_last <- z
    z <- c(update$alpha, update$tee)
    #v <- list(b, -1)
    lastBoundry <- boundry
    boundry <- boundryCondition(a, b, v, z, c)
    print("last boundry")
    print((lastBoundry))
    print("new boundry")
    print(boundry)
    # Check boundry condition ... time to quit?
  }
  #print("Done with solving")
  #print(boundry)
  #print(z)
  return(c(alpha, tee))
}


# OK so here,
#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
#  run by solvede(r, s, h, y, c) 
runSolverTest1 <- function() {
  epsilon <- 0.005                              # Estimated error of data set
  c <- 4                                        # c can be defined from epsilon
  s <- c                                        # step
  r <- rep(0, times=500)                        # This is the margin r(x)
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  y <- h                                        # This is the true class
  i = sample(x=seq(1, 500), replace=F, size=150)
  y[i] <- -1 * y[i]                             # Give 30% error to y(x)
  result <- solvede(r, s, h, y, c)
  print(result)
}


#  Test Case:  What if hypothesis and prediction match perfectly.
#  .........
#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
runSolverTest2 <- function() {
  alpha <- 0.0
  tee <- 0.0
  epsilon <- 0.1                               # Estimated error of data set
  c <- (1/as.double(erfc(1-epsilon)))^2         # c defined from epsilon
  c <- 32
  s <- c                               # s == c at t == 0
  a <-rep(0, times=500) + s                     # This is the margin r(x) + step
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  y <- h                                        # This is the true class
  i = sample(x=seq(1, 500), replace=F, size=40)
  y[i] <- -1 * y[i]                             # Give some error to y(x)
  b <- h * y                                    # agreement == 1, disagreement == -1
  result <- solvede(alpha, tee, a, b, c)
  print(result)
}
