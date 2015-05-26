#' @title Plot inner model
#' 
#' @description 
#' Plot the inner (structural) model for objects of class \code{"plspm"}, 
#' as well as path matrices
#' 
#' @param x Either a matrix defining an inner model or an 
#' object of class \code{"plspm"}.
#' @param colpos Color of arrows for positive path coefficients.
#' @param colneg Color of arrows for negative path coefficients.
#' @param box.prop Length/width ratio of ellipses.
#' @param box.size Size of ellipses.
#' @param box.cex Relative size of text in ellipses.
#' @param box.col fill color of ellipses,
#' @param lcol border color of ellipses.
#' @param box.lwd line width of the box.
#' @param txt.col color of text in ellipses.
#' @param shadow.size Relative size of shadow of label box.
#' @param curve arrow curvature.
#' @param lwd line width of arrow.
#' @param arr.pos Relative position of arrowheads on arrows.
#' @param arr.width arrow width.
#' @param arr.lwd line width of arrow, connecting two different points, 
#' (one value, or a matrix with same dimensions as \code{x}).
#' @param cex.txt Relative size of text on arrows.
#' @param show.values should values be shown when \code{x} is a matrix.
#' @param \dots Further arguments passed on to \code{\link{plotmat}}.
#' @note \code{innerplot} uses the function
#' \code{\link{plotmat}} in package \code{diagram}. \cr
#' \url{http://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf}
#' @seealso \code{\link{plot.plspm}}, \code{\link{outerplot}}
#' @export 
innerplot <-
function(x, colpos = "#6890c4BB", colneg = "#f9675dBB",
           box.prop = 0.55, box.size = 0.08, box.cex = 1, box.col = "gray95", 
           lcol="gray95", box.lwd = 2, txt.col = "gray50", shadow.size = 0,
           curve = 0, lwd = 3, arr.pos = 0.5, arr.width = 0.2, arr.lwd = 3,
           cex.txt = 0.9, show.values = FALSE,
           ...)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (!inherits(x, "plspm") && !inherits(x, "matrix"))
    stop(paste("\nSorry, can't work with an object of class "), class(x))
  # if x is "plspm"
  if (inherits(x, "plspm"))
  {
    # get ingredients
    IDM = x$model$IDM
    lvs = nrow(IDM)
    # matrix of path coefficients
    MPC = x$path_coefs[lvs:1,]
    MPC = MPC[,lvs:1]
    names = rownames(MPC)
    # arrow matrix colors
    AM.col = MPC
    AM.col[MPC < 0] = colneg # negative path coeffs in red
    AM.col[MPC >= 0] = colpos # positive path coeffs in blue      
  } else {
    if (nrow(x) != ncol(x))
      stop("\nSorry, the provided matrix is not a square matrix")
    lvs = nrow(x)
    if (is.null(rownames(x))) {
      if (is.null(colnames(x))) {
        rownames(x) = as.character(1:lvs)
      } else {
        rownames(x) = colnames(x)
      }        
    } else {
      colnames(x) = rownames(x)
    }
    MPC = x[lvs:1,]
    MPC = MPC[,lvs:1]
    names = rownames(MPC)
    AM.col = MPC
    AM.col[MPC < 0] = colneg # negative path coeffs in red
    AM.col[MPC >= 0] = colpos # positive path coeffs in blue
    if (!show.values) cex.txt = 0
  }
  # check arr.lwd
  if (is.matrix(arr.lwd))
  {
    # Arrow Line Width
    ALW = arr.lwd[lvs:1,]
    ALW = ALW[,lvs:1]
    arr.lwd = ALW
  }

  # plot of inner model (adapted function from plotmat)
  diagram::plotmat(round(MPC, 4),     # square coefficient matrix
          name = names,               # names of elements
          box.type = "ellipse",       # shape of label box
          box.size = box.size,        # size of label box
          box.prop = box.prop,        # length/width ratio of label box
          box.col = box.col,          # fill color of label box
          lcol = lcol,                # color of box line
          box.lwd = box.lwd,          # line width of the box
          box.cex = box.cex,          # relative size of text in boxes
          txt.col = txt.col,          # color of text in boxes
          curve = curve,              # arrow curvature
          lwd = lwd,                  # line width of arrow
          cex.txt = cex.txt,          # relative size of arrow text
          arr.type = "triangle",      # type of arrowhead
          arr.pos = arr.pos,          # relative pos of arrowhead on arrow line
          arr.lwd = arr.lwd,          # width of arrow, connecting two points
          shadow.size = shadow.size,  # relative size of shadow of label box
          prefix = "",                # added in front of non-zero arrow labels
          arr.lcol = AM.col,          # color of arrow line
          arr.col = AM.col,           # color of arrowhead
          arr.width = arr.width,      # arrow width
          self.arrpos = pi/2,         # position of the self-arrow
          ...) 
}
