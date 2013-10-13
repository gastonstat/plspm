\name{get_manifests}
\alias{get_manifests}
\title{Building data matrix with manifest variables}
\usage{
  get_manifests(Data, blocks)
}
\arguments{
  \item{Data}{matrix or data frame from which manifest
  variables are extracted}

  \item{blocks}{list (wtih numeric elements) indicating the
  set of manifest variables that form each block}
}
\value{
  matrix or data frame of selected manifest variables
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_manifests} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

