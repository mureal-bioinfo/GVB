#' GVB_preCalculation() Function
#' 
#' This function is for calculating the GVB score.
#' @import psych
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
#' GVB_preCalculation(input = VDF_list, score = "SIFT", threshold = 0.7, direction = "or_less", homozygote = TRUE, method = "hm", indelScore)

GVB_preCalculation <- function(input, score, threshold, direction, homozygote, method, indelScore)
{
  gene <- input[1,'Gene']
  sample_id <- input[1,'Sample_ID']
  
  input[, score] <- replace(input[, score], input[, score]==0, 10^-8)

#score
  if(score=="SIFT")
  {
    if(is.null(indelScore))
    {
      input <- input[!is.na(input[, score]), ]
    }
    
    else if(is.numeric(indelScore))
    {
      input[, score] <- replace(input[, score], is.na(input[, score]), indelScore)
    }
    
    else
    {
      stop("Invalid indelScore!")
    }   
  }
  
  else if(score=="PPT2")
  {
    if(is.null(indelScore))
    {
      input <- input[!is.na(input[, score]), ]
    }
    
    else if(is.numeric(indelScore))
    {
      input[, score] <- replace(input[, score], is.na(input[, score]), indelScore)
    }
    
    else
    {
      stop("Invalid indelScore!")
    }
  }

  else if(score=='CADD')
  {
    input <- input
  }


#direction
  if(direction=='or_less')
  {
    input_threshold <- input[input[, score] <= threshold, ]
  }
    
  else if(direction=='or_more')
  {
    input_threshold <- input[input[, score] >= threshold, ]
  }
  
  else
  {
    stop("Invalid value for direction. It must be 'or_less' or 'or_more'.")
  }
  


  n <- nrow(input_threshold)
  
  if(n==0)
  { 
    GVB_score <- 1
  }
  
  else if(n>0)
  {
  #homozygote    
    if(homozygote==TRUE)
    {
      homo_index <- grep(2, input_threshold$Genotype)
      het_index <- grep(1, input_threshold$Genotype)
      
      homo_v <- (input_threshold[score][homo_index, ])^0.5
      het_v <- input_threshold[score][het_index, ]
      
      if(method=="gm")
      {
        GVB_score <- geometric.mean(c(het_v, homo_v))
        GVB_score <- round(GVB_score,4)
      }
      
      else if(method=="hm")
      {
        GVB_score <- harmonic.mean(c(het_v, homo_v))
        GVB_score <- round(GVB_score,4)
      }

      else
      {
        stop("Invalid value for method. It must be 'gm' or 'hm'.")
      }
    
    }
    
    else if(homozygote==FALSE)
    { 
      carrier_index <- grep("1|2", input_threshold$Genotype)
      
      carrier_v <- input_threshold[score][carrier_index, ]
      
      if(method=="gm")
      {
        GVB_score <- geometric.mean(carrier_v)
        GVB_score <- round(GVB_score,4)
      }
      
      else if(method=="hm")
      {
        GVB_score <- harmonic.mean(carrier_v)
        GVB_score <- round(GVB_score,4)
      }
      
      else
      {
        stop("Invalid value for method. It must be 'gm' or 'hm'.")
      }
    }
    
  }
  
  GVB_DataFrame <- data.frame(Gene=gene, Sample_ID=GVB_score)
  colnames(GVB_DataFrame)[grep("Sample_ID", colnames(GVB_DataFrame))] <- sample_id
  
  return(GVB_DataFrame)
}

