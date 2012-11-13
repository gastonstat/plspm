\name{plspm.groups}
\alias{plspm.groups}
\title{Two Groups Comparison in PLS-PM}
\usage{
  plspm.groups(pls, group, Y = NULL, method = "bootstrap",
    reps = NULL)
}
\arguments{
  \item{pls}{object of class \code{"plspm"}}

  \item{group}{factor with 2 levels indicating the groups
  to be compared}

  \item{Y}{optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}

  \item{method}{method to be used in the test. Possible
  values are \code{"bootstrap"} or \code{"permutation"}}

  \item{reps}{integer indicating the number of either
  bootstrap resamples or number of permutations. If
  \code{NULL} then \code{reps}=100}
}
\value{
  An object of class \code{"plspm.groups"}

  \item{test}{Table with the results of the applied test.
  Includes: path coefficients of the global model, path
  coeffs of group1, path coeffs of group2, (absolute)
  difference of path coeffs between groups, and the test
  results with the p-value.}

  \item{global}{List with results of the inner model for
  the global model}

  \item{group1}{List with results of the inner model for
  group1}

  \item{group2}{List with results of the inner model for
  group2}
}
\description{
  Performs a group comparison test for comparing path
  coefficients between two groups. The null and alternative
  hypotheses to be tested are: H0: path coefficients are
  not significantly different; H1: path coefficients are
  significantly different
}
\details{
  \code{plspm.groups} performs a two groups comparison test
  in PLS-PM for comparing path coefficients between two
  groups. Only two methods are available: 1) bootstrap, and
  2) permutation. The bootstrap test is an adapted t-test
  based on bootstrap standard errors. The permutation test
  is a randomization test which provides a non-parametric
  option.

  When the object \code{pls} does not contain a data matrix
  (i.e. \code{pls$data=NULL}), the user must provide the
  data matrix or data frame in \code{Y}.
}
\examples{
\dontrun{
 ## example with customer satisfaction analysis
 ## group comparison based on the segmentation variable "gender"

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

 # permutation test with 100 permutations
 group_comp = plspm.groups(my_pls, satisfaction$gender,
                           method="permutation", reps=100)
 group_comp
 }
}
\author{
  Gaston Sanchez
}
\references{
  Chin, W.W. (2003) A permutation procedure for multi-group
  comparison of PLS models. In: Vilares M., Tenenhaus M.,
  Coelho P., Esposito Vinzi V., Morineau A. (Eds.)
  \emph{PLS and Related Methods - Proceedings of the
  International Symposium PLS03.} Decisia, pp. 33-43.

  Chin, W.W. (2000) Frequently Asked Questions, Partial
  Least Squares PLS-Graph. Available from:
  \url{http://disc-nt.cba.uh.edu/chin/plsfaq/multigroup.htm}
}
\seealso{
  \code{\link{plspm}}
}

