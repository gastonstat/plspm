\name{rho}
\alias{rho}
\title{Dillon-Goldstein's rho}
\usage{
  rho(X)
}
\arguments{
  \item{X}{matrix representing one block of manifest
  variables}
}
\value{
  Dillon-Goldstein's rho
}
\description{
  Dillon-Goldstein's rho of a single block of variables
}
\examples{
\dontrun{
 # load dataset satisfaction
 data(satisfaction)

 # block Image (first 5 columns of satisfaction)
 Image = satisfaction[,1:5]

 # compute Dillon-Goldstein's rho for Image block
 rho(Image)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{alpha}}, \code{\link{unidim}}
}

