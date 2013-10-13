\name{check_path}
\alias{check_path}
\title{Check path matrix}
\usage{
  check_path(path_matrix)
}
\arguments{
  \item{path_matrix}{square (lower triangular) boolean
  matrix for inner model}
}
\value{
  validated path matrix
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_path} is called by
  \code{check_args}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

