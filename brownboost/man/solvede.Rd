\name{solvede}
\alias{solvede}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
solvede(alpha, tee, a, b, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{tee}{
%%     ~~Describe \code{tee} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(alpha, tee, a, b, c) {
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
