\name{get_unidim}
\alias{get_unidim}
\title{Unidimensionality of reflective blocks}
\usage{
  get_unidim(DM, blocks, modes)
}
\arguments{
  \item{DM}{Data Matrix}

  \item{blocks}{vector with numbers of variables per block}

  \item{modes}{vector of modes}
}
\value{
  A data frame with the following columns:

  \item{Type.measure}{Measurement mode}

  \item{MVs}{number of manifest variables in each block}

  \item{C.alpha}{Cronbach's alpha}

  \item{DG.rho}{Dillon-Goldstein rho}

  \item{eig.1st}{First eigenvalue}

  \item{eig.2nd}{Second eigenvalue}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_unidim} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

