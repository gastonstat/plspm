#' @S3method print plspm.groups
print.plspm.groups <-
function(x,...)
{
    cat("GROUP COMPARISON IN PLS-PM FOR PATH COEFFICIENTS", "\n\n")
    cat("Scale of Data:        ", x$settings[[1]], "\n")
    cat("Weighting Scheme:     ", x$settings[[2]], "\n")
    cat("Selected method:      ", x$settings[[3]], "\n")
    cat("Number of replicates: ", x$reps, "\n\n")
    cat("$test", "\n")
    print(x$test, print.gap=2)
    cat("\n")
    cat("Inner models in the following objects:", "\n")
    cat("$global ", "\n")
    cat("$group1 ", "\n")
    cat("$group2 ", "\n")
}

