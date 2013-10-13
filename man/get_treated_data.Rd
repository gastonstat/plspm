\name{get_treated_data}
\alias{get_treated_data}
\title{Apply corresponding treatment to MV}
\usage{
  get_treated_data(MV, specs)
}
\arguments{
  \item{MV}{matrix or data frame from with manifest
  variables}

  \item{specs}{list with algorithm specifications}
}
\value{
  matrix or data frame of treated MVs (scaling)
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_treated_data} is called by
  \code{plspm}.
}
\note{
  In non-metric case, all mvs are starndardized, and those
  ordinal or nominal are rankified
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

