\name{test_dataset}
\alias{test_dataset}
\title{Test Data Set Availibility}
\usage{
  test_dataset(Dataset, pls_data, num_obs)
}
\arguments{
  \item{Dataset}{optional data set (with manifest
  variables) provided by user}

  \item{pls_data}{Data matrix containing the manifest
  variables used in the model. Only available when
  \code{dataset=TRUE} inside \code{plspm()}}

  \item{num_obs}{number of rows in PLS-PM Scores}
}
\value{
  TRUE if dataset is available, otherwise FALSE
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{test_dataset} checks whether a
  dataset is available for plspm.groups, rescale, it.reb,
  rebus.pls, rebus.test
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

