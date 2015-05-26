#' @title Plot outer model
#' 
#' @description Plot either outer weights or loadings in the 
#' outer model for objects of class \code{"plspm"}
#' 
#' @param x An object of class \code{"plspm"}.
#' @param what What to plot: \code{"loadings"} or \code{"weights"}.
#' @param colpos Color of arrows for positive path coefficients.
#' @param colneg Color of arrows for negative path coefficients.
#' @param box.prop Length/width ratio of ellipses and rectangles.
#' @param box.size Size of ellipses and rectangles.
#' @param box.cex Relative size of text in ellipses and rectangles.
#' @param box.col fill color of ellipses and rectangles.
#' @param lcol border color of ellipses and rectangles.
#' @param box.lwd line width of the box.
#' @param txt.col color of text in ellipses and rectangles.
#' @param shadow.size Relative size of shadow of label box.
#' @param curve arrow curvature.
#' @param lwd line width of arrow.
#' @param arr.pos Relative position of arrowheads on arrows.
#' @param arr.width arrow width.
#' @param cex.txt Relative size of text on arrows.
#' @param \dots Further arguments passed on to \code{\link{plotmat}}.
#' @note \code{outerplot} uses the function
#' \code{\link{plotmat}} of package \code{diagram}. \cr
#' \url{http://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf}
#' @seealso \code{\link{innerplot}}, \code{\link{plot.plspm}}, 
#' \code{\link{plspm}}
#' @export 
outerplot <-
function(x, what = "loadings", colpos = "#6890c4BB", colneg = "#f9675dBB",
         box.prop = 0.55, box.size = 0.08, box.cex = 1, box.col = "gray95", 
         lcol = "gray95", box.lwd = 2, txt.col = "gray40", shadow.size = 0,
         curve = 0, lwd = 2, arr.pos = 0.5, arr.width = 0.15, cex.txt = 0.9,
         ...)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (!inherits(x, "plspm"))
    stop("\nSorry, an object of class 'plspm' is required")
  check_what = what %in% c("loadings", "weights")
  if (!check_what) what = "loadings"
  
  # get ingredients
  IDM = x$model$IDM
  blocks = x$model$blocks    
  modes = x$model$specs$modes
  loadings = x$outer_model$loading
  out.weights = x$outer_model$weight
  lvs = nrow(IDM)
  mvs_names = as.character(x$outer_model$name)
  # rows and columns
  rs = c(1,1,1,2,2,2,2,2,3,2,3,3)
  cs = c(1,2,3,2,3,3,4,4,3,5,4,4)
  index.mat = cbind(1:12, rs, cs)
  # auxiliary indices
  ini_end = from_to(lengths(blocks))
  ini.vec = ini_end$from
  end.vec = ini_end$to 
  
  # =======================================================
  # Plotting
  # =======================================================
  # set graphical parameters
  op = par(mfrow = index.mat[lvs,2:3], mar=c(0, 3, 2.5, 2)) 
  # for each block
  for (k in 1:lvs)
  {
    num.mvs = length(blocks[[k]])
    names.mvs = mvs_names[ini.vec[k]:end.vec[k]]
    names.mvs = c(names.mvs, rownames(IDM)[k])
    box.types = c(rep("rect", num.mvs), "ellipse")
    # matrix with either loadings or weights
    MAT = matrix(0, num.mvs+1, num.mvs+1)
    MAT.col = MAT
    if (what == "loadings") {
      # loadings
      MAT[num.mvs+1,] = c(loadings[ini.vec[k]:end.vec[k]], 0)
      if (modes[k] == "A")  MAT = t(MAT) 
    } else {
      # outer weights
      MAT[num.mvs+1,] = c(out.weights[ini.vec[k]:end.vec[k]], 0)
    }
    # set colors
    MAT.col[MAT < 0] = colneg  # negative
    MAT.col[MAT >= 0] = colpos  # positive

    # call plotmat
    diagram::plotmat(round(MAT, 4),     # square matrix with values
            name = names.mvs,           # names of manifest variables
            box.type = box.types,       # shape of label box
            box.size = box.size,        # size of label box
            box.prop = box.prop,        # ength/width ratio of label box
            box.col = box.col,          # fill color of label box
            lcol = lcol,                # color of box line
            box.lwd = box.lwd,          # line width of the box
            box.cex = box.cex,          # relative size of text in boxes
            txt.col = txt.col,          # color of text in boxes
            shadow.size = shadow.size,  # relative size of shadow of label box
            curve = curve,              # arrow curvature
            lwd = lwd,                  # line width of arrow
            cex.txt = cex.txt,          # relative size of arrow text
            arr.type = "triangle",      # type of arrowhead
            arr.pos = arr.pos,          # relative pos of arrowhead on arrow line
            prefix = "",                # added in front of non-zero arrow labels
            arr.lcol = MAT.col,         # color of arrow line
            arr.col = MAT.col,          # color of arrowhead
            arr.width = arr.width,      # arrow width
            main = c(paste(rownames(IDM)[k]), what),
            self.arrpos = pi/2,         # position of the self-arrow
            ...)
  } # end for
  # reset graphical parameters
  par(op)
}
