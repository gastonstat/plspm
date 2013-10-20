#' @S3method summary plspm
summary.plspm <- function(object, ...)
{
  ## Reminder of model in objects "plspm"
  # y$model = list(IDM=path_matrix, blocks=blocks, specs=specs, 
  #                iter=weights$iter, boot.val=boot.val, br=br, gens=gens)
  ## Reminder of model in objects "plspm.fit"
  # y.fit$model = list(IDM=path_matrix, blocks=blocks, specs=specs, 
  #                    iter=weights$iter, gens=gens)
  
  # =======================================================
  # inputs setting
  # =======================================================  
  y = object
  IDM = y$model$IDM
  exo_endo = rep("Exogenous", nrow(IDM))
  exo_endo[rowSums(IDM) != 0] = "Endogenous"

  inputs = data.frame(Block = rownames(IDM), 
                      Type = exo_endo, 
                      Size = lengths(y$model$blocks), 
                      Mode = y$model$specs$modes,
                      row.names = 1:length(exo_endo))
  
  # =======================================================
  # results
  # =======================================================  
  if (length(y$model) == 5) 
  {
    # results for object "plspm.fit"
    res = list(inputs = inputs, 
               outer_model = y$outer_model, 
               inner_model = y$inner_model, 
               model = y$model)
  } else {
    # results for object "plspm"
    correlations = round(cor(y$scores), 4)
    res = list(inputs = inputs, 
               unidim = y$unidim, 
               outer_model = y$outer_model, 
               crossloadings = y$crossloadings, 
               inner_model = y$inner_model, 
               correlations = correlations,
               inner_summary = y$inner_summary, 
               gof = y$gof, 
               effects = y$effects, 
               boot = y$boot, 
               model = y$model)    
  }
  class(res) = "summary.plspm"
  res
}


#' @S3method print summary.plspm
print.summary.plspm <- function(x, ...)
{
  ## Reminder of model in objects "plspm"
  # x$model <- list(IDM, blocks, scheme, modes, scaled, boot.val, 
  #                 plsr, obs, br, tol, iter, n.iter, outer)
  ## Reminder of model in objects "plspm.fit"
  # x.fit$model <- list(IDM, blocks, scheme, modes, scaled, 
  #                     obs, tol, iter, n.iter, outer)
  
  # =======================================================
  # inputs setting
  # =======================================================  
  if (x$model$specs$scaled) Scale="Standardized Data" else Scale="Raw Data"
  
  cat("PARTIAL LEAST SQUARES PATH MODELING (PLS-PM)", "\n\n")
  cat("----------------------------------------------------------", "\n")   
  cat("MODEL SPECIFICATION", "\n")
  cat("1   Number of Cases     ", x$model$gens$obs, "\n")
  cat("2   Latent Variables    ", nrow(x$model$IDM), "\n")
  cat("3   Manifest Variables  ", x$model$gens$mvs, "\n")
  cat("4   Scale of Data       ", Scale, "\n")
  cat("5   Non-Metric PLS      ", !get_metric(x$model$specs$scaling), "\n")
  cat("6   Weighting Scheme    ", x$model$specs$scheme, "\n")
  cat("7   Tolerance Crit      ", x$model$specs$tol, "\n")
  cat("8   Max Num Iters       ", x$model$specs$maxiter, "\n")
  cat("9   Convergence Iters   ", x$model$iter, "\n")
  if (length(x$model) > 5)
  {
    boot_sam <- if(is.null(x$model$br)) "NULL" else x$model$br
    cat("10  Bootstrapping       ", x$model$boot.val, "\n")
    cat("11  Bootstrap samples   ", boot_sam, "\n")
  }
  cat("\n")
  cat("----------------------------------------------------------", "\n")    
  cat("BLOCKS DEFINITION", "\n")
  print(x$inputs, print.gap = 3)
  cat("\n")
  cat("----------------------------------------------------------", "\n") 
  if (length(x$model)>5) 
  {   
    cat("BLOCKS UNIDIMENSIONALITY","\n")
    print(x$unidim, print.gap = 2, digits = 3)
    cat("\n")
    cat("----------------------------------------------------------", "\n")    
  }
  cat("OUTER MODEL","\n")
  out_mod_mat = NULL
  out_mod_labs = NULL
  tmp = indexify(x$model$blocks)
  for (i in 1:length(x$model$blocks)) {
    out_mod_mat = rbind(out_mod_mat, rep(NA, 4), 
                        as.matrix(x$outer_model[tmp==i,3:6]))
    out_mod_labs = c(out_mod_labs, x$model$gens$lvs_names[i],
          paste(" ", i, x$model$gens$mvs_names[tmp==i]))
  }
  rownames(out_mod_mat) = out_mod_labs
  print(out_mod_mat, na.print="", print.gap = 2, digits = 3)
#  print(x$outer_model, print.gap = 2, digits = 3)
  cat("\n")
  cat("----------------------------------------------------------", "\n")    
  if (length(x$model) > 5)
  {
    cat("CROSSLOADINGS","\n")
    cros = NULL
    for (i in 1:length(x$model$blocks)) {
      cros = rbind(cros, rep(NA, ncol(x$crossloadings[,-c(1,2)])),
                   as.matrix(x$crossloadings[tmp==i,-c(1,2)]))
    } 
    rownames(cros) = out_mod_labs
    print(cros, na.print="", print.gap=2, digits=3)
#    print(x$crossloadings, print.gap = 2, digits = 3)
    cat("\n")
    cat("----------------------------------------------------------", "\n") 
  }
  cat("INNER MODEL","\n")
  print(x$inner_model, print.gap = 3, digits = 3)
  if (length(x$model) > 5)
  {
    cat("----------------------------------------------------------", "\n") 
    cat("CORRELATIONS BETWEEN LVs","\n")
    print(x$correlations, print.gap = 2, digits = 3)
    cat("\n")
    cat("----------------------------------------------------------", "\n") 
    cat("SUMMARY INNER MODEL","\n")
    print(x$inner_summary, print.gap = 2, digits = 3)
    cat("\n")
    cat("----------------------------------------------------------", "\n") 
    cat("GOODNESS-OF-FIT","\n")
    print(x$gof, print.gap = 2, digits = 4)
    cat("\n")
    cat("----------------------------------------------------------", "\n")   
    cat("TOTAL EFFECTS","\n")
    print(x$effects, print.gap = 2, digits = 3)
    if (!is.logical(x$boot))
    {
      cat("\n")
      cat("---------------------------------------------------------", "\n")    
      cat("BOOTSTRAP VALIDATION", "\n")
      for (i in 1:length(x$boot))
      {
        cat(names(x$boot)[i], "\n")
        print(x$boot[[i]], print.gap=2, digits=3)
        cat("\n")
      }
    }      
  }
  invisible(x)
}
