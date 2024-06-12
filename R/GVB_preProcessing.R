#' GVB_preProcessing() Function
#'
#' This function preprocesses the input and converts it to a list where each element is a data frame for each gene before calculating the GVB score.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @param x a data frame derived from individual variant annotated file.
#' @param indel a logical value. If TRUE, include indel variants (Default = FALSE).
#' @export
#' @examples
#' GVB_preProcessing(x, indel=FALSE)

GVB_preProcessing <- function(x, indel=FALSE)
{
  if(indel==TRUE)
  {
    x <- x 
  }
  
  else
  {
    sub_x <- subset(x, (Ref!="") & (Ref!="-") & (Alt!="-"))
    x <- sub_x %>% filter(nchar(Ref) == nchar(Alt))
  }
  
  gene_list <- split(x, x$Gene)
  
  return(gene_list)
}
