
# contains the function erf(a)
library(NORMT3)

# dot product of two vectors
dot <- function (v1, v2)	{
    return(sqrt(sum(v1*v2)))
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
            
  f1 <- sum(b * exp( -(1/c) * (a + dot(v, z))^2))
  f2 <- sum(erfd((a + dot(v, z))/sqrt(c)) - erfd(a/sqrt(c)))
  return(data.frame(f1=f1,f2=f2))
}

# when calculating the jacobian matrix of the boundry condition
# functions, we need elements W, U, B, and V.
# when calculating the updates for alpha and t we need E
# v = [b, -1] and z = [a, t]
jacobianElements <- function(a, b, v, z, c) {
  d <- a + dot(v,z)
  w <- exp(-d^2/c)
  W <- sum(w)
  U <- sum(w*d*b)
  B <- sum(w*b)
  V <- sum(w*d*b*b)
  E <- sum(erfd(d/sqrt(c)) - erfd(a/sqrt(c)))
  return(data.frame(W=W,U=U,B=B,V=V,E=E))
}


updateStep <- function(alpha, t, c, W, U, B, V, E) {
  if (W != 0) {
    Bk <- B/W
    Vk <- V/W
    Uk <- U/W
    a2 <- alpha + (c*Bk + sqrt(pi*c)*Uk*(E/W)) / (2*(Vk - Uk*Bk))
    t2 <- t + (c*Bk^2 + sqrt(pi*c)*Vk*(E/W)) / (2*(Vk - Uk*Bk))
  } else {
    a2 = t2 = 0
  }

  return(data.frame(alpha=a2, t=t2))
}


solvede <- function(alpha, t, a, b, c) {
  z <- c(alpha, t)
  v <- c(b, rep(-1, times=length(b)))
  for (i in 1:10) {
    print(i)
    j <- jacobianElements(a, b, v, z, c)
    print(j)
    update <- updateStep(alpha, t, c, j$W, j$U, j$B, j$V, j$E)
    print(update)
    alpha <- update$alpha
    t <- update$t
    z <- c(alpha, t)
    v <- c(b, rep(-1, times=length(b)))
    boundry <- boundryCondition(a, b, v, z, c)
    print(boundry)
  }
}


# OK so here,
# a_j == r_i(x_j)
# b_j == h_i(x_j)y_j + s
# .... at least I think the s is there...

runSolver <- function() {
  alpha <- 0.1
  t <- alpha*alpha/3
  epsilon <- 0.05
  c <- (1/as.double(erfc(1-epsilon)))^2
  a <- abs(rnorm(50, 10, 1))
  b <- abs(rnorm(50, 30, 1))
  solvede(alpha, t, a, b, c)
}
