\name{get_num_scale}
\alias{get_num_scale}
\title{Non-Metric Numerical Scale}
\usage{
  get_num_scale(X)
}
\arguments{
  \item{X}{a matrix to be scaled}
}
\value{
  scaled matrix
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_num_scale} is called by
  \code{plspm}.
}
\note{
  scales a matrix X in such a way that mean(X[,j])=0 and
  varpop(X[,j])=1 this means that sum(X[,j]^2) = n if MD,
  sum(X[,j]^2, na.rm=T) = number of available elements
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

