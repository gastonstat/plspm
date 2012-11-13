#'@S3method print rebus.test
print.rebus.test <-
function(x, ...)
{
    cat("Validation for REBUS-PLS by permutation tests", "\n")
    cat("---------------------------------------------", "\n")
    cat("List with the results for path coefficients, ", "\n") 
    cat("loadings, and GoF in the following objects:", "\n\n")
    cat("  Name    ", "  Classes", "\t", "   Elements in each object", "\n")
    for (k in 1:length(x))
    {
        cat(paste("$",names(x)[k],sep=""), "  ",  
            substr(names(x)[k], 6, 6), "and", substr(names(x)[k], 8, 8), 
            "  ", "...$paths,", "...$loadings,", "...$gof", "\n")
    }
}

