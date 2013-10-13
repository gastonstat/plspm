\name{unidim}
\alias{unidim}
\title{Unidimensionality of blocks}
\usage{
  unidim(Data, blocks = NULL)
}
\arguments{
  \item{Data}{matrix or data frame with variables}

  \item{blocks}{optional list with vectors indicating the
  variables in each block}
}
\value{
  A data frame with the following columns:

  \item{Block}{name of block}

  \item{MVs}{number of manifest variables in each block}

  \item{C.alpha}{Cronbach's alpha}

  \item{DG.rho}{Dillon-Goldstein rho}

  \item{eig.1st}{First eigenvalue}

  \item{eig.2nd}{Second eigenvalue}
}
\description{
  Compute unidimensionality indices (a.k.a. Composite
  Reliability indices)
}
\examples{
\dontrun{
 # load dataset satisfaction
 data(satisfaction)

 # blocks Image and Expectations
 ima_expe = list(Image=1:5, Expec=6:10)

 # compute unidimensionality indices
 unidim(satisfaction, ima_expe)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{alpha}}, \code{\link{rho}}
}

