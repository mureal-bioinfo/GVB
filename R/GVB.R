#' GVB() Function
#' 
#' This function calculates the GVB score.
#' @importFrom data.table rbindlist
#' @param input a list, the output of GVB_preProcessing() function.
#' @param score a character vector for in-silico variant score (e.g.: "SIFT", "PP2", "CADD") which is a column name of the data frame as elements of 'input' (Default = "SIFT").
#' @param threshold a numeric vector specifying cut-off of the in-silico variant score. The default value for "SIFT" is 0.7, for "PP2" it's 0.85, and for "CADD" it's 15.
#' @param direction The direction of scores to be included based on the threshold. Only one of "or_less" or "or_more" is allowed. "or_less" means <= the threshold and "or_more" means >= the threshold (Default = "or_less").
#' @param homozygote logical value. If TRUE, consider homozygotes for GVB calculation. (Default = FALSE).
#' @param method allowed values are one of "gm" or "hm". "gm" means geometric mean and "hm" means harmonic mean (Default = "gm").
#' @param indelScore a numeric value to assign the in-silico variant score of indel variants when the 'score' parameter is "SIFT" or "PP2" (Default = NULL).
#' @export
#' @examples
#' VDF_list <- GVB_preProcessing(VariantDataFrame)
#' GVB(input = VDF_list, score = "SIFT", threshold = 0.7, direction = "or_less", homozygote = TRUE, method = "hm", indelScore)

GVB <- function(input, score = "SIFT", threshold = 0.7, direction = "or_less", homozygote = FALSE, method = "gm", indelScore = NULL)
{

  if(score=="PP2")
  {
    threshold <- if(missing(threshold)) 0.85 else threshold
  }

  else if(score=="CADD")
  {
    threshold <- if(missing(threshold)) 15 else threshold
  }

  GVB_Result <- as.data.frame(rbindlist(lapply(X = input, FUN = GVB_preCalculation, score = score, threshold = threshold, direction = direction, homozygote = homozygote, method = method, indelScore = indelScore)))
  
  return(GVB_Result)
}

