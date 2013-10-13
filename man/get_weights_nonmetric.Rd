\name{get_weights_nonmetric}
\alias{get_weights_nonmetric}
\title{Outer Weights Non-Metric Data}
\usage{
  get_weights_nonmetric(X, path_matrix, blocks, specs)
}
\arguments{
  \item{X}{scaled data}

  \item{path_matrix}{matrix with path connections}

  \item{blocks}{list with variables in each block}

  \item{specs}{list with algorithm specifications}
}
\value{
  list of outer weights, ODM, scores, QQ, iter
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_weights} is called by
  \code{plspm}
}
\note{
  Calculate outer weights for non-metric data (under
  Lohmoller's algorithm)
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

