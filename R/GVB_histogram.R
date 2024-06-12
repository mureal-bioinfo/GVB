#' GVB_histogram() Function
#' 
#' This function plots the GVB score for one specific gene as a histogram from merged GVB data, which is the output of the GVB() function.
#' @param data a matrix, the output of GVB_merging() function.
#' @param geneSymbol a character vector as gene symbol which is a column name of the 'data'.
#' @param outputDirectory a character vector as directory for saving the plot. The file name of plot is geneSymbol_histogram.jpg.
#' @export
#' @examples
#' GVB_histogram(data = merged_gvb, geneSymbol = "Gene1")

GVB_histogram <- function(data, geneSymbol, outputDirectory=NULL)
{

  data_df <- as.data.frame(data)

  MEAN <- round(mean(data_df[,geneSymbol]),2)
  
  hist_data <- hist(x = data_df[,geneSymbol], freq = TRUE, xlab = "GVB", ylab = "Count", main = geneSymbol, plot = FALSE, warn.unused = FALSE)
  MaxCount <- round(max(hist_data$counts)*1.2)

  if(!(is.null(outputDirectory)))
  {
    OutputFilePath <- paste0(outputDirectory, "/", geneSymbol, "_histogram.jpg")
    
    jpeg(OutputFilePath, quality = 100, width = 2400, height = 1600, res = 300)
    
    hist(x = data_df[,geneSymbol], freq = TRUE, breaks = seq(0, 1, by = 0.05), xlab = "GVB", ylab = "Count", main = geneSymbol, ylim = c(0, MaxCount))
    abline(v = MEAN, lty = 2, col = 'red', lwd = 2)
    
    dev.off()
    
    cat(OutputFilePath, "\n")
  }
  
  
  else
  {
    hist(x = data_df[,geneSymbol], freq = TRUE, breaks = seq(0, 1, by = 0.05), xlab = "GVB", ylab = "Count", main = geneSymbol, ylim = c(0, MaxCount))
    abline(v = MEAN, lty = 2, col = 'red', lwd = 2)
  }
  
}