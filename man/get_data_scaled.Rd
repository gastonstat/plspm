\name{get_data_scaled}
\alias{get_data_scaled}
\title{Scaling data outside \code{plspm}}
\usage{
  get_data_scaled(DM, scaled)
}
\arguments{
  \item{DM}{matrix or data frame from with manifest
  variables}

  \item{scaled}{whether to scale latent variables}
}
\value{
  matrix or data frame of (un)scaled MVs
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_scaled_data} is called by
  \code{plspm.groups} and others.
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

