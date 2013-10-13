\name{test_manifest_scaling}
\alias{test_manifest_scaling}
\title{Test scaling of selected manifest variables}
\usage{
  test_manifest_scaling(MV, scaling)
}
\arguments{
  \item{MV}{matrix or data frame with selected manifest
  variables}

  \item{scaling}{list with specified scaling for each
  manifest variable}
}
\value{
  gives an error when factors in MV have bad scaling
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{test_manifest_scaling} is called
  by \code{plspm}. Basically checks that unordered factors
  have nominal scaling.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

