\name{get_dummy}
\alias{get_dummy}
\title{Non-Metric Dummy}
\usage{
  get_dummy(x)
}
\arguments{
  \item{x}{a vector whose elements are all the natural
  numbers from 1 to p}
}
\value{
  the dummy matrix
}
\description{
  Internal function. \code{get_dummy} is called by
  \code{get_dummies}. Transforms a vector of natural
  numbers from 1 to p into the corresponding p x p dummy
  matrix
}
\keyword{internal}

