\name{test_null_weights}
\alias{test_null_weights}
\title{Test outer weights convergence within specified maxiter}
\usage{
  test_null_weights(wgs, specs)
}
\arguments{
  \item{wgs}{output from 'get_weights()'}

  \item{specs}{list of specifications}
}
\value{
  gives an error when wgs is NULL (i.e. non-convergence)
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{test_weights_output} is called
  by \code{plspm}. Basically checks that weights are not
  NULL.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

