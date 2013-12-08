#' @title Quantification Plot
#' 
#' @description
#' Quantification Plots for Non-Metric PLS-PM
#' 
#' @details
#' If both \code{lv} and \code{mv} are specified, only the value of \code{lv}
#' will be taken into account. \cr
#' If the given \code{lv} have more than 15 variables, only the first 15
#' are plotted.
#' 
#' @param pls a non-metric \code{"plspm"} object
#' @param lv number or name of latent variable
#' @param mv number or name of manifest variable
#' @param pch Either an integer specifying a symbol or a single character to be
#' used as the default in plotting points
#' @param col color
#' @param lty type of line
#' @param \dots Further arguments passed on to \code{\link{plot}}.
#' @export
quantiplot <- 
function(pls, lv = NULL, mv = NULL, pch = 16, col = "darkblue", lty = 2, ...) 
{
  if (class(pls) != "plspm" && is.null(pls$model$specs$scaling))  
    stop("\n'quantiplot()' requires a non-metric plspm object")
  # verify lv and mv
  check_lv_mv(lv, mv, pls$model$gens$lvs_names, pls$model$gens$mvs_names) 
  
  # plots of a block (latent variable)
  if (!is.null(lv)) {
    quantiplot_lv(pls, lv, pch=pch, col=col, lty=lty, ...)
  } else {
    # plot of a manifest variable
    if (!is.null(mv)) {
      quantiplot_mv(pls, mv, pch=pch, col=col, lty=lty, ...)
    }     
  }
  
}

check_lv_mv <- function(lv, mv, lvs_names, mvs_names) 
{
  if (is.null(lv) & is.null(mv))
    stop("\n'quantiplot()' requires non-null arguments 'lv' or 'mv'")
  
  # check lv value
  if (!is.null(lv)) {
    # lv as character
    if (is.character(lv)) {
      if (is.na(match(lv, lvs_names)))
        stop("'lv' name not recognized")
    } else {
      # lv as numeric
      if (is.numeric(lv)) {
        if(lv > length(lvs_names))
          stop("'lv' value out of limits")
      }  
    }
  }
  
  # check mv value
  if (!is.null(mv)) {
    # mv as character
    if (is.character(mv)) {
      if (is.na(match(mv, mvs_names)))
        stop("'mv' name not recognized")
    } else {
      # mv as numeric
      if (is.numeric(mv)) {
        if(mv > length(mvs_names))
          stop("'mv' value out of limits")
      }  
    }
  }
  # output
  TRUE
}

quantiplot_lv <- 
function(pls, lv, pch, col, lty, ...)
{
#  # fixed rows and columns vectors to set a layout window 
#  rs = c(1, 1, 1, 2, 2, 2, 2, 2, 3, 2, 3, 3)
#  cs = c(1, 2, 3, 2, 3, 3, 4, 4, 3, 5, 4, 4)
#  # matrix for graphical layout 
#  layout_mat = cbind(1:12, rs, cs)
  
  # fixed rows and columns vectors to set a layout window 
  rs = c(1, 1, 1, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3)
  cs = c(1, 2, 3, 2, 3, 3, 4, 4, 3, 5, 4, 4, 5, 5, 5)
  # matrix for layout 
  layout_mat = cbind(1:15, rs, cs)
  
  # get original and quantified mvs in 'lv' block
  mvs = pls$model$blocks[[lv]]
  
  # No more than 15 mvs can be plotted
  if (length(mvs) > 15) {
    warning("'quantiplot()' can show only 15 indicators")
    mvs = mvs[1:15]
  }
  
  mvs_original = pls$data[,mvs]
  mvs_quantified = pls$manifests[,mvs]
  # if only one mv then convert to matrix
  if (length(mvs) == 1) {
    mvs_original = as.matrix(mvs_original)
    mvs_quantified = as.matrix(mvs_quantified)
  }
  
  # graphical parameters
  op = par(mfrow = layout_mat[length(mvs), 2:3], mar = c(4.5, 4, 3, 2))
  # plot indicators
  for (q in seq_along(mvs)) 
  {
    plot(mvs_original[,q], mvs_quantified[,q], 
         xlab = "raw values", ylab = "scaling values", 
         pch = pch, col = col, ...)
    # don't add lines if mv is what_scale == "nom"
    what_scale = unlist(pls$model$specs$scaling)[mvs[q]]
    if (what_scale != "nom") {
      lines(sort(mvs_original[,q]), sort(mvs_quantified[,q]), 
            col = col, lty = lty)    
    }    
    title(main = colnames(pls$data)[mvs[q]], ...)   
  }
  # reset graphical parameters
  par(op)
}

quantiplot_mv <- 
function(pls, mv, pch, col, lty, ...)
{
  # column index of 'mv'
  mv_num = mv
  if (is.character(mv)) {
    mv_num = match(mv, pls$model$gens$mvs_names)
  }
  # extract original and quantified mv
  mv_original = pls$data[,mv_num]
  mv_quantified = pls$manifests[,mv_num]
  
  # for MV
  what_scale = unlist(pls$model$specs$scaling)[mv_num]
  
  # plot values
  plot(mv_original, mv_quantified, 
       xlab = "raw values", ylab = "scaling values", 
       pch = pch, col = col, ...)
  # don't add lines if mv is what_scale == "nom"
  if (what_scale != "nom") {
    lines(sort(mv_original), sort(mv_quantified), 
          col = col, lty = lty)    
  }
  # add title
  title(main = colnames(pls$data[mv]), ...)
}
