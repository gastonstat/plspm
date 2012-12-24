\name{rescale}
\alias{rescale}
\title{Rescale Latent Variable Scores}
\usage{
  rescale(pls, Y = NULL)
}
\arguments{
  \item{pls}{object of class \code{"plspm"}}

  \item{Y}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}
}
\value{
  A matrix with the rescaled latent variable scores
}
\description{
  Rescale standardized latent variable scores to original
  scale of manifest variables
}
\details{
  \code{rescale} requires all outer weights to be positive
}
\examples{
\dontrun{
 ## example with customer satisfaction analysis

 # load data satisfaction
 data(satisfaction)

 # define inner model matrix
 IMAG = c(0,0,0,0,0,0)
 EXPE = c(1,0,0,0,0,0)
 QUAL = c(0,1,0,0,0,0)
 VAL = c(0,1,1,0,0,0)
 SAT = c(1,1,1,1,0,0)
 LOY = c(1,0,0,0,1,0)
 sat.inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

 # define outer model list
 sat.outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

 # define vector of reflective modes
 sat.mod = rep("A", 6)

 # apply plspm
 my_pls = plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="factor",
              scaled=FALSE)

 # rescaling standardized scores of latent variables
 scores = rescale(my_pls)

 # compare standardized LVs against rescaled LVs
 summary(my_pls$latents)
 summary(scores)
 }
}
\author{
  Gaston Sanchez
}
\seealso{
  \code{\link{plspm}}
}

