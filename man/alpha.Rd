\name{alpha}
\alias{alpha}
\title{Cronbach's alpha}
\usage{
  alpha(X)
}
\arguments{
  \item{X}{matrix representing one block of manifest
  variables}
}
\value{
  Cronbach's alpha
}
\description{
  Cronbach's alpha of a single block of variables
}
\examples{
\dontrun{
 # load dataset satisfaction
 data(satisfaction)

 # block Image (first 5 columns of satisfaction)
 Image = satisfaction[,1:5]

 # compute Cronbach's alpha for Image block
 alpha(Image)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{rho}}, \code{\link{unidim}}
}

