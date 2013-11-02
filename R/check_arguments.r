#' @title Check arguments for \code{plspm} and \code{plspm.fit}
#' 
#' @details
#' Internal function. \code{check_args} is called by \code{plspm}.
#'
#' @param Data numeric matrix or data frame containing the manifest variables.
#' @param path_matrix square (lower triangular) boolean matrix for inner model.
#' @param blocks List of vectors (either numeric or character) to specify
#' the blocks of manifest variables asociated to the latent variables.
#' @param scaling optional list of string vectors indicating the type of 
#' measurement scale for each manifest variable specified in \code{blocks}.
#' @param modes character indicating the type of measurement.
#' @param scheme string indicating the type of inner weighting scheme.
#' @param scaled logical indicating whether scaling data is performed.
#' @param tol decimal value indicating the tol criterion for covergence.
#' @param maxiter integer indicating the maximum number of iterations.
#' @param plscomp optional vector indicating the number of PLS components
#' (for each block) to be used when handling non-metric data 
#' (only used if \code{scaling} is provided)
#' @param boot.val whether bootstrap validation is performed.
#' @param br integer indicating the number bootstrap resamples.
#' @param dataset whether the data matrix should be retrieved.
#' @return list of validated arguments
#' @keywords internal
#' @template internals
#' @export
check_args <- 
function(Data, path_matrix, blocks, scaling, modes, scheme,
         scaled, tol, maxiter, plscomp, boot.val, br, dataset)
{
  # check definitions
  Data = check_data(Data)
  path_matrix = check_path(path_matrix)
  blocks = check_blocks(blocks, Data)
  specs = check_specs(blocks, scaling, modes, scheme, scaled, 
                      tol, maxiter, plscomp)
  boot_args = check_boot(boot.val, br)
  if (!is.logical(dataset)) dataset = TRUE
  
  # check congruence between inner model and outer model
  good_model = check_model(path_matrix, blocks)
  
  # list with verified arguments
  list(Data = Data,
       path_matrix = path_matrix,
       blocks = blocks,
       specs = specs,
       boot.val = boot_args$boot.val,
       br = boot_args$br, 
       dataset = dataset)
}


#' @title Check Data
#' 
#' @details
#' Internal function. \code{check_data} is called by \code{check_args}.
#'
#' @param Data numeric matrix or data frame containing the manifest variables.
#' @return validated Data 
#' @keywords internal
#' @template internals
#' @export
check_data <- function(Data)
{
  if (is_not_tabular(Data))
    stop("\nInvalid 'Data'. Must be a matrix or data frame.")
  
  if (is.matrix(Data) && !is.numeric(Data))
    stop("\nInvalid 'Data' matrix. Must be a numeric matrix.")
  
  if (nrow(Data) == 1)
    stop("\nCannot work with only one row in 'Data'")
  
  if (ncol(Data) == 1)
    stop("\nCannot work with only one column in 'Data'")
  
  if (lacks_rownames(Data))
    rownames(Data) = 1:nrow(Data)
  
  if (lacks_colnames(Data)) 
    colnames(Data) = paste("MV", 1:ncol(Data), sep="")
  
  # return
  Data
}


#' @title Check path matrix
#' 
#' @details
#' Internal function. \code{check_path} is called by \code{check_args}.
#'
#' @param path_matrix square (lower triangular) boolean matrix for inner model
#' @return validated path matrix 
#' @keywords internal
#' @template internals
#' @export
check_path <- function(path_matrix)
{
  if (is_not_matrix(path_matrix))
    stop("\n'path_matrix' must be a matrix.")
  
  if (!is_square_matrix(path_matrix))
    stop("\n'path_matrix' must be a square matrix.")
  
  if (nrow(path_matrix) == 1)
    stop("\n'path_matrix' must have more than one row")
  
  if (!is_lower_triangular(path_matrix))
    stop("\n'path_matrix' must be a lower triangular matrix")
  
  
  for (j in 1:ncol(path_matrix)) 
  {
    for (i in 1:nrow(path_matrix)) 
    {
      if (length(intersect(path_matrix[i,j], c(1,0))) == 0)
        stop("\nElements in 'path_matrix' must be '1' or '0'")
    }      
  }
  
  if (lacks_dimnames(path_matrix)) {
    LV_names = paste("LV", 1:ncol(path_matrix), sep = "")
    dimnames(path_matrix) = list(LV_names, LV_names)
  }
  if (has_rownames(path_matrix) && lacks_colnames(path_matrix)) {
    colnames(path_matrix) = rownames(path_matrix)
  }
  if (has_colnames(path_matrix) && lacks_rownames(path_matrix)) {
    rownames(path_matrix) = colnames(path_matrix)
  }
  
  # return
  path_matrix
}


#' @title Check well defined blocks
#' 
#' @details
#' Internal function. \code{check_blocks} is called by \code{check_args}.
#'
#' @param blocks list defining the blocks of manifest variables. Elements in
#' \code{blocks} must have the same mode: all "numeric" or all "character".
#' @param Data matrix or data frame from where to extract manifest variables.
#' @return validated blocks (output in numeric format)
#' @keywords internal
#' @template internals
#' @export
check_blocks <- function(blocks, Data)
{
  if (!is.list(blocks))
    stop("\n'blocks' must be a list.")
  
  # no duplicated elements within each block
  mvs_duplicated = unlist(lapply(blocks, duplicated))
  if (any(mvs_duplicated))
    stop("\nWrong 'blocks'. Duplicated variables in a block are not allowed")
  
  # all elements in blocks of same mode
  mvs_mode = unique(unlist(lapply(blocks, mode)))
  if (length(mvs_mode) > 1)
    stop("\nAll elements in 'blocks' must have the same mode")
  
  # check indices inside columns range of Data
  if (mvs_mode == "numeric") {
    blocks_in_data = match(unlist(blocks), 1:ncol(Data))
    if (any(is.na(blocks_in_data)))
      stop("\nIndices in 'blocks' outside the number of columns in 'Data'")
  }
  
  # convert character blocks to numeric blocks
  if (mvs_mode == "character") {
    data_names = colnames(Data)
    matched_names = match(unlist(blocks), data_names)
    if (any(is.na(matched_names))) {
      bad_names = unlist(blocks)[is.na(matched_names)]
      stop(sprintf("\nUnrecognized name in 'blocks': '%s'", bad_names))        
    }
    blocks = lapply(blocks, function(x, y) match(x, y), data_names)
  }
  
  # output
  blocks
}


#' @title Check bootstrap options
#' 
#' @details
#' Internal function. \code{check_boot} is called by \code{check_args}.
#'
#' @param boot.val logical indicating whether to perform bootstrapping
#' @param br number of bootstrap resamples
#' @return validated bootstrap options
#' @keywords internal
#' @template internals
#' @export
check_boot <- function(boot.val, br)
{
  if (!is.logical(boot.val)) boot.val = FALSE

  if (boot.val) {
    if (!is.null(br)) {
      if(!is_positive_integer(br) || length(br) != 1L || br < 10) {
        warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
        br = 100
      } 
    } else
      br = 100
  }
  
  # return
  list(boot.val = boot.val, br = br)
}


#' @title Check congruence between inner and outer models
#' 
#' @details
#' Internal function. \code{check_model} is called by \code{check_args}.
#'
#' @param path_matrix matrix specifying the path connections
#' @param blocks list defining the blocks of manifest variables
#' @return an error if there is something wrong 
#' @keywords internal
#' @template internals
#' @export
check_model <- function(path_matrix, blocks)
{    
  # compatibility between path_matrix and blocks
  if (length(blocks) != nrow(path_matrix))
    stop("\nNumber of rows in 'path_matrix' different from length of 'blocks'.")
  
  # output
  TRUE
}
