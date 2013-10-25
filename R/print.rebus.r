#'@S3method print rebus
print.rebus <- function(x, ...)
{
  lvs <- x$aux[[1]]
  n.end <- x$aux[[2]]
  segs <- table(x$segments)
  segmen <- rbind(segs, 100*segs/sum(segs))
  dimnames(segmen) <- list(c("number.units","proportions(%)"), 
                           colnames(x$path.coef))
  cat("\n")
  cat("RESPONSE-BASED UNIT SEGMENTATION (REBUS)", "\n")
  cat("IN PARTIAL LEAST SQUARES PATH MODELING", "\n")
  cat("----------------------------------------------", "\n\n")    
  cat("Parameters Specification", "\n")
  cat("  Number of segments:   ", length(segs), "\n")
  cat("  Stop criterion:       ", x$aux[[4]], "\n")
  cat("  Max number of iter:   ", x$aux[[5]], "\n\n")
  cat("REBUS solution (on standardized data)", "\n")
  cat("  Number of iterations: ", x$aux[[6]], "\n")
  cat("  Rate of unit change:  ", x$aux[[3]], "\n")
  cat("  Group Quality Index:  ", x$aux[[7]], "\n\n")
  cat("REBUS Segments", "\n")
  print(segmen, print.gap=3, digits=1)
  cat("\n")
  cat("----------------------------------------------", "\n")   
  cat("$path.coef", "\n")
  print(round(x$path.coef,4), print.gap=3, justify="right")
  cat("\n")
  cat("----------------------------------------------", "\n")    
  cat("$loadings", "\n")
  print(round(x$loadings,4), print.gap=3, justify="right")
  cat("\n")
  cat("----------------------------------------------", "\n")    
  qual.index <- NULL
  qual.labs <- NULL
  ini.ind <- c(1,lvs+1,lvs+n.end+1,lvs+2*n.end+1)
  end.ind <- c(lvs,lvs+n.end,lvs+2*n.end,lvs+2*n.end+1)
  qual.nams <- c("Aver.Com","Aver.Redu","R2","GoF")
  for (i in 1:4) {        
    q <- ini.ind[i]:end.ind[i]
    qual.index <- rbind(qual.index, rep(NA,ncol(x$quality)), x$quality[q,])
    qual.labs <- c(qual.labs, qual.nams[i],
                   paste(rep(" ", length(rownames(x$quality)[q])), 
                         rownames(x$quality)[q]))
  }
  rownames(qual.index) <- qual.labs
  cat("$quality", "\n")
  print(qual.index, na.print="", print.gap=2)    
}
