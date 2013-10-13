\name{get_metric}
\alias{get_metric}
\title{Type of metric based on scaling measurement}
\usage{
  get_metric(scaling)
}
\arguments{
  \item{scaling}{list with measurement scale of each
  manifest variable}
}
\value{
  metric type of metric (metric data = TRUE / non-metric
  data = FALSE)
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_metric} is called by
  \code{plspm}. It is used to decide what function apply
  for the iterative process of the PLS-PM algorithm.
}
\note{
  If scaling is NULL, apply get_weights(), otherwise apply
  get_weights_nonmetric()
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

