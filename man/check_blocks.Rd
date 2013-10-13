\name{check_blocks}
\alias{check_blocks}
\title{Check well defined blocks}
\usage{
  check_blocks(blocks, Data)
}
\arguments{
  \item{blocks}{list defining the blocks of manifest
  variables. Elements in \code{blocks} must have the same
  mode: all "numeric" or all "character".}

  \item{Data}{matrix or data frame from where to extract
  manifest variables.}
}
\value{
  validated blocks (output in numeric format)
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_blocks} is called by
  \code{check_args}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

