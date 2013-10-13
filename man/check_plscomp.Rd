\name{check_plscomp}
\alias{check_plscomp}
\title{Check vector of PLS components (for non-metric plspm)}
\usage{
  check_plscomp(plscomp, scaling, modes)
}
\arguments{
  \item{plscomp}{optional vector indicating the number of
  PLS components (for each block) to be used when handling
  non-metric data (only used if \code{scaling} is
  provided)}

  \item{scaling}{list with type of scaling for each
  manifest variable}
}
\value{
  validated plscomp
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_plscomp} is called by
  \code{check_specs}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

