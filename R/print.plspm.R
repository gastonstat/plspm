#'@S3method print plspm
print.plspm <-
function(x, ...)
{
  cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n")
  cat(rep("-", 45), sep="")
  cat("\n$outer.mod     ", "outer model")
  cat("\n$inner.mod     ", "inner model")
  cat("\n$latents       ", "scaled latent variables (LVs)")
  cat("\n$scores        ", "LVs for scaled=FALSE")
  cat("\n$out.weights   ", "outer weights")
  cat("\n$loadings      ", "loadings")
  cat("\n$path.coefs    ", "path coefficients matrix")
  cat("\n$r.sqr         ", "R-squared")  
  if (!inherits(x, "plspm.fit"))
  {
    cat("\n$outer.cor     ", "outer correlations")
    cat("\n$inner.sum     ", "summary inner model")
    cat("\n$effects       ", "total effects")
    cat("\n$unidim        ", "unidimensionality")
    cat("\n$gof           ", "goodnes-of-fit")
    cat("\n$boot          ", "bootstrap results")
    cat("\n$data          ", "data matrix")
  }
  cat("\n")
  cat(rep("-", 45), sep="")
  cat("\nYou can also use the function 'summary'", "\n\n")    
  invisible(x)
}
