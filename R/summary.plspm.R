#' @S3method summary plspm
summary.plspm <-
function(object, ...)
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
  y = object
  IDM = y$model$IDM
  blocks = y$model$blocks
  modes = y$model$modes
  Mode = modes
  Mode[modes == "A"] = "Reflective"
  Mode[modes == "B"] = "Formative"   
  exo.endo = rowSums(IDM)
  exo.endo[rowSums(IDM) == 0] = "Exogenous"
  exo.endo[rowSums(IDM) != 0] = "Endogenous"
  blocklist = as.list(1:sum(blocks))
  for (j in 1:length(blocks))
    blocklist[[j]] = rep(j,blocks[j])
  blocklist = unlist(blocklist)
  inputs = data.frame(Block = rownames(IDM), 
                      Type = exo.endo, 
                      NMVs = y$model$blocks, 
                      Mode = Mode)
  rownames(inputs) = 1:length(exo.endo) 

  # =======================================================
  # results
  # =======================================================  
  if (length(y$model) == 10) 
  {
    # results for object "plspm.fit"
    res = list(inputs = inputs, 
               outer.mod = y$outer.mod, 
               inner.mod = y$inner.mod, 
               xxx = y$model)
  } else {
    # results for object "plspm"
    lat.cor = round(cor(y$latents), 4)
    res = list(inputs = inputs, 
               unidim = y$unidim, 
               outer.mod = y$outer.mod, 
               outer.cor = y$outer.cor, 
               inner.mod = y$inner.mod, 
               latent.cor = lat.cor,
               inner.sum = y$inner.sum, 
               gof = y$gof, 
               effects = y$effects, 
               boot = y$boot, 
               xxx = y$model)    
  }
  class(res) = "summary.plspm"
  res
}

