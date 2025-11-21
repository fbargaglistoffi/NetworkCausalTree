#' @title
#' Generate all combinations of the supplied vectors, without rows characterized
#' by equal elements
#' @description
#' Create a data frame including all combinations of the supplied vectors or factors,
#' while omitting rows characterized by equal elements.
#' @param  x A vector
#' @param  y A vector
#' @param  include.equals Boolean (dafault: FALSE)
#'
#' @return A data frame with all combinations of the elements of the vectors,
#' with no rows characterized by equal elements
#'
#' @import dplyr
#'
expand.grid.unique <- function(x, y, include.equals = FALSE){
  
  x <- unique(x)
  y <- unique(y)
  
  g <- function(i){
    z <- dplyr::setdiff(y, x[seq_len(i - include.equals)])
    if (length(z)) as.data.frame(cbind(x[i], z, deparse.level = 0))  # ← ensure data.frame
  }
  
  result <- do.call(rbind, lapply(seq_along(x), g))
  
  # Ensure final output is always a data.frame
  if (!is.null(result)) result <- as.data.frame(result)
  
  return(result)
}

