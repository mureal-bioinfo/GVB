#' GVB_merging() Function
#'
#' This function merges GVB data frames of multiple samples derived from GVB() function.
#' @param x a vector that includes object names of GVB data frames as global variables.
#' @export
#' @examples
#' GVB_merging(x = GVB_DataFrames_Vector)

GVB_merging <- function(x)
{
  GVBs_list <- mget(x, envir = globalenv())
  
  merged_GVB <- Reduce(function(a,b) merge(a,b, all=TRUE), GVBs_list)
  merged_GVB[is.na(merged_GVB)] <- 1
  
  rownames(merged_GVB) <- merged_GVB$Gene
  #merged_GVB <- as.data.frame(t(subset(merged_GVB, select = -c(Gene))))
  merged_GVB <- t(subset(merged_GVB, select = -c(Gene)))
  
  
  return(merged_GVB)
}
