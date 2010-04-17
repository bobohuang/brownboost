
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
  return((1/pi)*as.double(erf(a)))
}

# when the solution is found, this result should be c(0,0)
# returns true if each function returns zero
# v = [b, -1] and z = [a, t]

boundryCondition <- function (a, b, v, z, c) {
  #print(a)
  #print(b)
  #print(v)
  #print(z)
  #print(c)
  #print(erfd((a + dot(v, z))))
  #print(erfd(a/sqrt(c)))
            
  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  f2 <- sum(erfd((a + dotl(z, v))/sqrt(c)) - erfd(a/sqrt(c)))
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
    a2 = t2 = 0
  }
  return(data.frame(alpha=a2, tee=t2))
}


solvede <- function(alpha, tee, a, b, c) {
  z <- c(alpha, tee)
  v <- list(b, -1)
  boundry <- boundryCondition(a, b, v, z, c)
  # Check boundry condition .... quit?
  for (i in 1:20) {
    j <- jacobianElements(a, b, v, z, c)
    update <- updateStep(alpha, tee, c, j$W, j$U, j$B, j$V, j$E)
    alpha <- update$alpha
    tee <- update$tee
    z <- c(alpha, tee)
    v <- list(b, -1)
    boundry <- boundryCondition(a, b, v, z, c)
    # Check boundry condition ... time to quit?
  }
  print(boundry)
  return(c(alpha, tee))
}


# OK so here,
#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
runSolver <- function() {
  alpha <- 0.0
  tee <- 0.0
  epsilon <- 0.005                              # Estimated error of data set
  c <- (1/as.double(erfc(1-epsilon)))^2         # c defined from epsilon
  s <- c                                        # s == c at t == 0
  a <-rep(0, times=500) + s                     # This is the margin r(x) + step
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  y <- h                                        # This is the true class
  i = sample(x=seq(1, 500), replace=F, size=150)
  y[i] <- -1 * y[i]                             # Give 30% error to y(x)
  b <- h * y                                    # agreement == 1, disagreement == -1
  result <- solvede(alpha, tee, a, b, c)
  print(result)
}
