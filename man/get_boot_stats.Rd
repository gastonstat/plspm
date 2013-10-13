\name{get_boot_stats}
\alias{get_boot_stats}
\title{Get data frame with bootstrap statistics}
\usage{
  get_boot_stats(MATRIX, original)
}
\arguments{
  \item{MATRIX}{Matrix with bootstrapped results}

  \item{original}{vector with original values}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_boot_stats} is called by
  \code{get_boots}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

