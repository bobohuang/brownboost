dyn.load("brownboost/src/bisection.so")

solvede <- function(r, s, h, y, c) {
  try (
    result <- .C("solvede", as.double(r), as.double(s), as.double(h),
                 as.double(y), as.double(c), as.integer(length(r)),
                 output=numeric(2))
       )
    return(result$output)
}

# dot product of one vector and one list
dotl <- function (v, l) {
  #print("dotl")
  #print(v[1]*l[[1]])
  #print( v[2]*l[[2]])
  #print(v[1]*l[[1]] + v[2]*l[[2]])
    return(v[1]*l[[1]] + v[2]*l[[2]])
}

myerf <- function(a) {
  return(2*pnorm(a*sqrt(2)) - 1)
}

# when the solution is found, this result should be c(0,0)
# returns true if each function returns zero
# v = [b, -1] and z = [a, t]
f <- function (a, b, v, z, c) {
  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  f2 <- sum(myerf((a + dotl(z, v))/sqrt(c)) - myerf(a/sqrt(c)))
  return(c(f1,f2))
}


createDataSet <- function(x=10) {
  dataset <- data.frame()
  
  for (i in 1:x) {
  
    exampleNum <- sample(seq(0, 1000, 1), 1)
    errorSize <- sample(seq(0.001, 0.5, 0.01), 1)
    c <- sample(seq(0.01, 10, by=0.01), 1) 
    s <- sample(seq(0.01, c, by=0.01), 1)                     # step
    rmean <- sample(seq(-3, 3, 0.001), 1)
    r <- rnorm(mean=rmean, sd=0.75, n=exampleNum)
    hone <- sample(1:100, 1)
    hseq <- c(rep(1, hone), rep(-1, 100-hone))
    h <- sample(x=hseq, replace=T, size=exampleNum)   # This is the hypothesis
    y <- h                                        # This is the true class
    i = sample(x=seq(1, exampleNum), replace=F, size=exampleNum * errorSize)
    y[i] <- -1 * y[i]
    a <- r + s;
    b <- h * y;
    v <- list(b, -1)
    result <- solvede(r, s, h, y, c)
    print(result)
    fcheck <- f(a, b, v, result, c)
    print(fcheck)
    rsummary <- summary(r)
    d <- list(exampleNum=exampleNum, errorSize=errorSize, C=c, S=s,
              sumA=sum(a), sumB=sum(b), sumH=sum(h), sumY=sum(y), sumR=sum(r),
              minR=rsummary["Min."], maxR=rsummary["Max."], medR=rsummary["Median"],
              meanR=rsummary["Mean"], percent1=(length(which(h==1))/length(h)),
              alpha=result[1], tee=result[2], check1=fcheck[1], check2=fcheck[2])
    dataset <- rbind(dataset, d)
  }
  return(dataset)
}
