#' @title Plots for PLS Path Models
#' 
#' @description
#' Plot method for objects of class \code{"plspm"}. This function plots either
#' the inner (i.e. structural) model with the estimated path coefficients, or
#' the outer (i.e. measurement) model with loadings or weights.
#' 
#' @details
#' \code{plot.plspm} is just a wraper of \code{\link{innerplot}} 
#' and \code{\link{outerplot}}.
#' 
#' @param x An object of class \code{"plspm"}.
#' @param what What to plot: \code{"inner"}, \code{"loadings"}, \code{"weights"}.
#' @param colpos Color of arrows for positive path coefficients.
#' @param colneg Color of arrows for negative path coefficients.
#' @param box.prop Length/width ratio of ellipses and rectangles.
#' @param box.size Size of ellipses and rectangles.
#' @param box.cex Relative size of text in ellipses and rectangles.
#' @param box.col fill color of ellipses and rectangles.
#' @param lcol border color of ellipses and rectangles.
#' @param txt.col color of text in ellipses and rectangles.
#' @param arr.pos Relative position of arrowheads on arrows.
#' @param cex.txt Relative size of text on arrows.
#' @param \dots Further arguments passed on to \code{\link{plotmat}}.
#' @note Function \code{plot.plspm} is based on the function
#' \code{\link{plotmat}} of package \code{diagram}. \cr
#' \url{http://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf}
#' @seealso \code{\link{innerplot}}, \code{\link{outerplot}}, 
#' \code{\link{plspm}}
#' @method plot plspm
#' @S3method plot plspm
#' @examples
#' 
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
#'  # load data satisfaction
#'  data(satisfaction)
#'  
#'  # define inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat.inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # define outer model list
#'  sat.outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'  
#'  # define vector of reflective modes
#'  sat.mod = rep("A", 6)
#'  
#'  # apply plspm
#'  satpls = plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="centroid", 
#'                scaled=FALSE)
#'                
#'  # plot path coefficients
#'  plot(satpls, what="inner")
#'  
#'  # plot loadings
#'  plot(satpls, what="loadings")
#'  
#'  # plot outer weights
#'  plot(satpls, what="weights")
#'  }
#'
plot.plspm <-
function(x, what="inner", colpos = "#6890c4BB", colneg = "#f9675dBB",
         box.prop=.55, box.size=0.08, box.cex=1, box.col="gray95", lcol="gray95",
         txt.col = "gray40", arr.pos=0.5, cex.txt=0.9, ...)
{
  if (what == "weights" || what == "loadings") {
    outerplot(x, what, colpos = colpos, colneg = colneg, box.prop = box.prop,
              box.size = box.size, box.cex = box.cex, box.col = box.col, 
              lcol=lcol, txt.col = txt.col, arr.pos = arr.pos, 
              cex.txt = cex.txt, ...)
    
  } else {
    innerplot(x, colpos = colpos, colneg = colneg, box.prop = box.prop,
              box.size = box.size, box.cex = box.cex, box.col = box.col, 
              lcol=lcol, txt.col = txt.col, arr.pos = arr.pos, 
              cex.txt = cex.txt, ...)
  }
}
