\name{is_missing}
\alias{is_missing}
\title{Presence of missing values}
\usage{
  is_missing(x)
}
\arguments{
  \item{x}{a list}
}
\value{
  vector indicating whether each element in x contains NAs
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{is_missing} is called by
  \code{get_weights_nonmetric}
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

