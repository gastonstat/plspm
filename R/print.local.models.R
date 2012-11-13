#'@S3method print local.models
print.local.models <-
function(x, ...)
{
    ## x$model = list(IDM, blocks, scheme, modes, scaled, boot.val, plsr, obs, br, tol, iter) 
    model <- x$glob.model$model
    boot.sam <- if(is.null(model$br)) "NULL" else model$br
    if (model[[5]]) Scale="Standardized Data" else Scale="Raw Data"
    n.clus <- length(x) - 1
    Name <- c("$glob.model", paste(rep("$loc.model",n.clus), 1:n.clus, sep="."))
    Description <- c("global model", paste(rep("local model class",n.clus), 1:n.clus, sep=" "))
    res1 <- cbind(Name, Description)    
    cat("This function calculates PLS-PM for global and local models", "\n")
    cat("-----------------------------------------------------------", "\n")    
    cat("  Number of classes", "\t", length(x)-1,"\n")
    cat("  Scale of Data", "\t", Scale, "\n")
    cat("  Weighting Scheme", "\t", model$scheme, "\n")
    cat("  Tolerance Crit", "\t", model$tol, "\n")
    cat("  Max Num of Iters", "\t", model$iter, "\n")
    cat("\n")
    cat("PLS-PM models in the following objects", "\n")
    rownames(res1) <- 1:nrow(res1)
    print(res1, print.gap=3, justify="right")
}

