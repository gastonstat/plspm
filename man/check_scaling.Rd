\name{check_scaling}
\alias{check_scaling}
\title{Check types of measurement scales and metric}
\usage{
  check_scaling(scaling, scaled, blocks)
}
\arguments{
  \item{scaling}{list with type of scaling for each
  manifest variable}

  \item{scaled}{should manifest variables be standardized}

  \item{blocks}{list defining the blocks of manifest
  variables}
}
\value{
  list with validated scaling and scaled
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_scaling} is called by
  \code{check_specs}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

