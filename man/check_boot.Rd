\name{check_boot}
\alias{check_boot}
\title{Check bootstrap options}
\usage{
  check_boot(boot.val, br)
}
\arguments{
  \item{boot.val}{logical indicating whether to perform
  bootstrapping}

  \item{br}{number of bootstrap resamples}
}
\value{
  validated bootstrap options
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_boot} is called by
  \code{check_args}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

