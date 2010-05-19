
# dot product of one vector and one list
dotl <- function (v, l) {
  #print("dotl")
  #print(v[1]*l[[1]])
  #print( v[2]*l[[2]])
  #print(v[1]*l[[1]] + v[2]*l[[2]])
    return(v[1]*l[[1]] + v[2]*l[[2]])
}

erf <- function(a) {
  return(2*pnorm(a*sqrt(2)) - 1)
}

# when the solution is found, this result should be c(0,0)
# returns true if each function returns zero
# v = [b, -1] and z = [a, t]

f <- function (a, b, v, z, c) {
  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  f2 <- sum(erf((a + dotl(z, v))/sqrt(c)) - erf(a/sqrt(c)))
  return(c(f1,f2))
}



fprint <- function (a, b, v, z, c) {
  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  print("fprint f1")
  print("z")
  print(z)
  print("v")
  print(v)
  print(dotl(z,v))
  print(a + dotl(z,v))
  print( (a + dotl(z, v))^2)
  print( -(1/c) * (a + dotl(z, v))^2)
  print(exp( -(1/c) * (a + dotl(z, v))^2))
  print(b * exp( -(1/c) * (a + dotl(z, v))^2))
  print(sum(b * exp( -(1/c) * (a + dotl(z, v))^2)))
  print("fprint f1 done....")

#####################################
  
  f2 <- sum(erf((a + dotl(z, v))/sqrt(c)) - erf(a/sqrt(c)))
  print(" ")
  print("NOW f2")
  print(sqrt(c))
  print(a/sqrt(c))
  print(erf(a/sqrt(c)))
  print(dotl(z, v))
  print(a + dotl(z, v))
  print((a + dotl(z, v))/sqrt(c))
  print(erf((a + dotl(z, v))/sqrt(c)))
  print(erf((a + dotl(z, v))/sqrt(c)) - erf(a/sqrt(c)))
  print(sum(erf((a + dotl(z, v))/sqrt(c)) - erf(a/sqrt(c))))
  return(c(f1,f2))
}


sameSign <- function (someList) {

  i <- c()
  for (j in 1:length(someList)) {
    if (sign(someList[[j]][1]) == sign(someList[[j]][2])) {
      i <- c(i, j)
    }
  }
  return(i)
}


closestToZero <- function(errors, pointIndex) {

  minDistance <- 1000000
  minIndex <- 0
  for (i in pointIndex) {
    if (sqrt(errors[[i]][1]^2 + errors[[i]][2]^2) <= minDistance) {
      minIndex <- i
      minDistance <- sqrt(errors[[i]][1]^2 + errors[[i]][2]^2)
    }
  }
  return(minIndex)
}


decideOnNewPoint <- function (x1, x2, a, b, v, c) {
  # Assume that the errors of x1 and x2 are both either (+, +) or (-, -)

  mid <- (x1+x2) / 2
  c1 <- c(mid[1], x2[2])
  c2 <- c(x2[1], mid[2])
  c3 <- c(mid[1], x1[2])
  c4 <- c(x1[1], mid[2])
  c5 <- c(x2[1], x1[2])
  c6 <- c(x1[1], x2[2])
    
  signerrx1  <- sign(f(a, b, v, x1, c)) # old errors
  signerrx2  <- sign(f(a, b, v, x2, c))
  errmid <- f(a,b, v, mid, c)
  errc1 <- f(a, b, v, c1, c)
  errc2 <- f(a, b, v, c2, c)
  errc3 <- f(a, b, v, c3, c)
  errc4 <- f(a, b, v, c4, c)
  errc5 <- f(a, b, v, c5, c)
  errc6 <- f(a, b, v, c6, c)

  points <- list(mid, c1, c2, c3, c4, c5, c6)
  errors <- list(errmid, errc1, errc2, errc3, errc4, errc5, errc6)

  pointIndex <- sameSign(errors)
  pointMin <- closestToZero(errors, pointIndex)

  if (pointMin == 0) {
    return(NULL)
  }
  
  p <- points[[pointMin]]
  e <- errors[[pointMin]]
  
  if (sign(e[1]) == signerrx1[1]) {
    return(list(p, x2))
  } else {
    return(list(x1, p))
  }      
  
  return(NULL)
}


getStartingPosition <- function(a, b, v, c) {

  # how many mismatches are there?
  percentError <- (length(which(b == -1)) / length(b)) * 100
  alphaSampleRegion <- (-0.55/42) * percentError + 0.55
  aSR <- alphaSampleRegion
  loops <- 0
  x1 <- c()
  x2 <- c()
  signX1 <- c(-1, 1)
  signX2 <- c(1, -1)
  while ((signX1[1] != signX1[2]) ||
         (signX2[1] != signX2[2]) ||
         (signX1[1] == signX2[1]) ||
         (signX1[2] == signX2[2])) {
    if (loops < 100) {
      x1[1] <- sample(seq(from=aSR, to=aSR+0.2, by=0.001), size=1)
      x1[2] <- sample(seq(from=0.2, to=0.4, by=0.001), size=1)
      x2[1] <- sample(seq(from=aSR-0.2, to=aSR, by=0.001), size=1)
      x2[2] <- sample(seq(from=0, to=0.2, by=0.001), size=1)
    } else {
      x1 <- sample(seq(from=0, to=2, by=0.001), size=2)
      x2 <- sample(seq(from=-2, to=0, by=0.001), size=2)
    }
    loops <- loops+1
    print(loops)
    signX1 <- sign(f(a, b, v, x1, c))
    signX2 <- sign(f(a, b, v, x2, c))
  }
  return(list(x1, x2))
}

#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
solvede <- function(r, s, h, y, c) {
  
  a <- r + s;
  b <- h * y;
  v <- list(b, -1)
  
  solved <- 0
  loopCounter <- c
  tries <- 0
  
  while (solved == 0) {
    points <- getStartingPosition (a, b, v, c)
    x1 <- points[[1]];  x2 <- points[[2]]

    while (x2[1] != x1[1] || x2[2] != x1[2]) {

      newPoints <- decideOnNewPoint(x1, x2, a, b, v, c)
    
      if (is.null(newPoints)) {
        err <- f(a, b, v, x1, c) 
        if (err[1] < 1e-10 && err[2] < 1e-10 && x1[2] > 0) {
          #cat("x1: ", x1, "err: ", f(a, b, v, x1, c), "\n")
          print(tries)
          return(x1)
        } else {
          break
        }
      } else if (loopCounter > 100) {
        loopCounter <- 0
        break
      }
      
      loopCounter <- loopCounter + 1
      x1 <- newPoints[[1]];  x2 <- newPoints[[2]];      
      #cat("loop: ", loopCounter, "\n")
      #cat("x1: ", x1, "err: ", f(a, b, v, x1, c), "\n")
      #cat("x2: ", x2, "err: ", f(a, b, v, x2, c),"\n")
    }
    tries <- tries + 1
  }

  print("Out of tries!")
  return(x1)
}

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
  print(f(r+s, h*y, list(h*y, -1), result, c))
  print(result)
}


runSolverTest2 <- function() {
  epsilon <- 0.005                              # Estimated error of data set
  c <- 4                                        # c can be defined from epsilon
  s <- c                                        # step
  r <- rep(0, times=500)                        # This is the margin r(x)
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  y <- h                                        # This is the true class
  i = sample(x=seq(1, 500), replace=F, size=50)
  y[i] <- -1 * y[i]                             # Give 10% error to y(x)
  result <- solvede(r, s, h, y, c)
  print(f(r+s, h*y, list(h*y, -1), result, c))
  print(result)
}


runSolverTestN <- function(M,N) {
  epsilon <- 0.005                              # Estimated error of data set
  c <- 4                                        # c can be defined from epsilon
  s <- c                                        # step
  r <- rep(0, times=500)                        # This is the margin r(x)
  h <- sample(x=c(-1,1), replace=T, size=500)   # This is the hypothesis
  listOfResults <- list()
  for (i in seq(M, N, by=5)) {
    
    y <- h                                        # This is the true class
    i = sample(x=seq(1, 500), replace=F, size=i)
    y[i] <- -1 * y[i]                             # Give 10% error to y(x)
    result <- solvede(r, s, h, y, c)
    listOfResults <- append(listOfResults, result)
  }
  return(unlist(listOfResults))
}



functionTest1 <- function() {
  epsilon <- 0.005                              # Estimated error of data set
  c <- 3                                      # c can be defined from epsilon
  s <- 3                                        # step
  r <- rep(0.3, times=50)
  s <- rep(0.01, times=50)                       # This is the margin r(x)
  h <- c(1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,1,1,1,1)
  y <- c(1,1,1,1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,-1,-1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,1,1,1,1,1)
  a <- r + s;
  b <- h * y;
  v <- list(b, -1)
  z <- c(0.028331,  0.597402)
  result <- f(a, b, v, z, c)
  print(result)

  z <- c(-0.285680, -0.319881)
  result <- f(a, b, v, z, c)
  print(result)
  
}
