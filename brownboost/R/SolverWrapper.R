dyn.load("brownboost/src/bisection.so")

solvede <- function(r, s, h, y, c) {
  try (
    result <- .C("solvede", as.double(r), as.double(s), as.double(h),
                 as.double(y), as.double(c), as.integer(length(r)),
                 output=numeric(2))
       )
    return(result$output)
}

