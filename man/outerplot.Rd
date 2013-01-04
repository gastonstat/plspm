\name{outerplot}
\alias{outerplot}
\title{Plot outer model}
\usage{
  outerplot(x, what = "loadings", colpos = "#6890c4BB",
    colneg = "#f9675dBB", box.prop = 0.55, box.size = 0.08,
    box.cex = 1, box.col = "gray95", lcol = "gray95",
    box.lwd = 2, txt.col = "gray40", shadow.size = 0,
    curve = 0, lwd = 2, arr.pos = 0.5, arr.width = 0.15,
    cex.txt = 0.9, ...)
}
\arguments{
  \item{x}{An object of class \code{"plspm"}.}

  \item{what}{What to plot: \code{"loadings"} or
  \code{"weights"}.}

  \item{colpos}{Color of arrows for positive path
  coefficients.}

  \item{colneg}{Color of arrows for negative path
  coefficients.}

  \item{box.prop}{Length/width ratio of ellipses and
  rectangles.}

  \item{box.size}{Size of ellipses and rectangles.}

  \item{box.cex}{Relative size of text in ellipses and
  rectangles.}

  \item{box.col}{fill color of ellipses and rectangles.}

  \item{lcol}{border color of ellipses and rectangles.}

  \item{box.lwd}{line width of the box.}

  \item{txt.col}{color of text in ellipses and rectangles.}

  \item{shadow.size}{Relative size of shadow of label box.}

  \item{curve}{arrow curvature.}

  \item{lwd}{line width of arrow.}

  \item{arr.pos}{Relative position of arrowheads on
  arrows.}

  \item{arr.width}{arrow width.}

  \item{cex.txt}{Relative size of text on arrows.}

  \item{\dots}{Further arguments passed on to
  \code{\link{plotmat}}.}
}
\description{
  Plot either outer weights or loadings in the outer model
  for objects of class \code{"plspm"}
}
\note{
  \code{outerplot} uses the function \code{\link{plotmat}}
  of package \code{diagram}. \cr
  \url{http://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf}
}
\examples{
\dontrun{
 ## typical example of PLS-PM in customer satisfaction analysis
 ## model with six LVs and reflective indicators
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
 satpls = plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="centroid",
               scaled=FALSE)

 # plot loadings
 outerplot(satpls, what="loadings")

 # plot outer weights
 outerplot(satpls, what="weights")
 }
}
\seealso{
  \code{\link{innerplot}}, \code{\link{plot.plspm}},
  \code{\link{plspm}}, \code{\link{plot.plspm}}
}

