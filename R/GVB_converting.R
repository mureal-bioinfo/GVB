#' GVB_converting() Function
#'
#' This function converts a VEP tab-delimited output to a GVB input.
#' @import dplyr
#' @param vep_tab a character vector for a file path of the VEP tab output.
#' @export
#' @examples
#' GVB_converting(vep_tab = "VEP-tab_Annotation.txt")

GVB_converting <- function(vep_tab)
{
  lines <- readLines(vep_tab)
  data_lines <- lines[!grepl("^##", lines)]
  
  header_index <- grep("^#", data_lines)
  header <- sub("#", "", data_lines[header_index])
  
  data_lines[header_index] <- header
  data_string <- paste(data_lines, collapse = "\n")
  
  data <- read.table(text = data_string, sep="\t", header=TRUE)
  
  
  #Remove variants w/o Gene symbol
  data <- subset(data, SYMBOL != '-')
  
  #Location to Chromsome and Start
  data$Chromosome <- sapply(strsplit(data$Location, ":"), function(x) x[[1]])
  
  pos <- sapply(strsplit(data$Location, ":"), function(x) x[[2]])
  data$Location <- sapply(strsplit(pos, "-"), function(x) x[[1]])
  colnames(data)[grep("^Location", colnames(data))] <- "Start"

  data <- data[grep("HET|HOM", data$ZYG),]
  
  #ZYG to Genotype and Sample_ID
  data$Sample_ID <- sapply(strsplit(data$ZYG, ":"), function(x) x[[1]])
  data$Genotype <- sapply(strsplit(data$ZYG, ":"), function(x) x[[2]])
  
  data$Genotype[grep("^HET$", data$Genotype)] <- 1
  data$Genotype[grep("^HOM$", data$Genotype)] <- 2
  
  
  gvb_input <- select(data, Chromosome, Start, REF_ALLELE, Allele, SYMBOL, Consequence, CADD_PHRED, SIFT, PolyPhen, Genotype, Sample_ID)
  
  column_names <- colnames(gvb_input)
  column_names[column_names == 'REF_ALLELE'] <- 'Ref'
  column_names[column_names == 'Allele'] <- 'Alt'
  column_names[column_names == 'SYMBOL'] <- 'Gene'
  column_names[column_names == 'CADD_PHRED'] <- 'CADD'
  column_names[column_names == 'PolyPhen'] <- 'PP2'
  colnames(gvb_input) <- column_names
  
  gvb_input$CADD <- gsub("-", NA, gvb_input$CADD)
  gvb_input$SIFT <- gsub("-", NA, gvb_input$SIFT)
  gvb_input$PP2 <- gsub("-", NA, gvb_input$PP2)

  gvb_input$Chromosome <- as.integer(gvb_input$Chromosome)
  gvb_input$Start <- as.integer(gvb_input$Start)
  gvb_input$CADD <- as.numeric(gvb_input$CADD)
  gvb_input$SIFT <- as.numeric(gvb_input$SIFT)
  gvb_input$PP2 <- as.numeric(gvb_input$PP2)
  gvb_input$Genotype <- as.integer(gvb_input$Genotype)
  
  N <- nrow(gvb_input)
  if(N>=1)
  {
    cat('"Variants without a gene symbol were excluded. As a result, there are ', N, ' variants remaining."\n')
  }
  
  else
  {
    cat('"Variants without a gene symbol were excluded. As a result, there is no variant."\n')
  }
  
  
  return(gvb_input)
}