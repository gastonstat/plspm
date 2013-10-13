\name{check_model}
\alias{check_model}
\title{Check congruence between inner and outer models}
\usage{
  check_model(path_matrix, blocks)
}
\arguments{
  \item{path_matrix}{matrix specifying the path
  connections}

  \item{blocks}{list defining the blocks of manifest
  variables}
}
\value{
  an error if there is something wrong
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_model} is called by
  \code{check_args}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

