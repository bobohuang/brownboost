
# dot product of one vector and one list
dotl <- function (v, l) {
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
  return(data.frame(f1=f1,f2=f2))
}

decideOnNewPoint <- function (x1, errx1, x2, errx2, mid, errmid) {

if (


}



#  a == r(x_j, y_j) + s_i :  The margin + the step s
#  and b == h(x_i) * y_i  :  The hypothesis * the prediction
solvede <- function(r, s, h, y, c) {
  x1 <- c(1, 1)
  x2 <- c(-1, -1)
  a <- r + s;
  b <- h * y;
  v <- list(b, -1)

  print(f(a, b, v, z1, c))
  print(f(a, b, v, z2, c))

  mid <- (x1+x2)/2

  while (mid != x1 && mid != x2) {

    errx1  <- f(a, b, v, x1, c)
    errmid <- f(a, b, v, mid, c)
    errx2  <- f(a, b, v, x2, c)
        
    listOfPoints <- decideOnNewPoint(x1, errx1, x2, errx2, mid, errmid)

    x1 <- list[[1]];  x2 <- list[[2]];      
    mid <- (x1+x2)/2
  }

  return(mid)
}

#  // Assumption: One of f(a) and f(b) is ≥ 0 and the other is ≤ 0
#  if f(a) <= 0 then
#      lo := a; hi := b
#  else
#      lo := b; hi := a
#  endif
# 
#  mid := lo + (hi-lo)/2
#  while (mid ≠ lo) and (mid ≠ hi) do
#     if f(mid) ≤ 0 then
#        lo := mid
#     else
#        hi := mid
#     endif
#     mid := lo + (hi-lo)/2
#  endwhile
# 
#  return mid

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
  print(f(r+s, h*y, list(h*y, -1), result, c))
  print(result)
}


