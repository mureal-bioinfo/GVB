GVB
================

### 1. Introduction

#### 1-1. GVB
Current gene-level approaches including the residual variation intolerance score (RVIS) and the gene damage index (GDI) use a population dataset to find pathogenic genes by assigning a score for each gene for a given population. 
Thus, they provide population-dependent and collective evaluation scores that are not individualized. 
On the other hand, [gene-wise variant burden](https://www.tandfonline.com/doi/full/10.2217/pgs-2020-0039) (GVB) assigns a score for each gene for each individual in a population-independent and individual-specific manner, which is necessary for clinical applications.

#### 1-2. GVB R package
GVB is an integrated gene-level measure of the cumulative impact of the multitude of deleterious variants on a given gene in a population data-free manner. 
This package facilitates GVB score calculation from *in silico* deleteriousness scores and visualization through [Ensembl Variant EffectPredictor](https://link.springer.com/article/10.1186/S13059-016-0974-4) (VEP) tab-delimited output.

<br>

### 2. Installation
``` r
install.packages("devtools")
library(devtools)
install_github("mureal-bioinfo/GVB")
```

<br>

### 3. Field requirements in VEP tab-delimited output (VEP-tab)

VEP-tab contains many fields ranging from genomics locations to *in
silico* deleteriousness scores. GVB package uses following fields.

- **Mandatory fields** (VEP version 111)

  | Filed       | Description                                      |
  |:------------|:-------------------------------------------------|
  | Location    | Coordinate format (chr:start or chr:start-end)   |
  | REF_ALLELE  | Reference allele                                 |
  | Allele      | Variant allele used to calculate the consequence |
  | SYMBOL      | Gene symbol                                      |
  | Consequence | Consequence type of this variant                 |
  | ZYG         | Zygosity of individual genotype at this locus    |
  | SIFT        | SIFT score                                       |

- **Optional fields:** Other *in silico* deleteriousness score of
  variants such as CADD, Polyphen, etc.

<br>

### 4. Reading VEP-tab files and convert these to the GVB input format

<span style="background-color: #EDEDED;">GVB_converting</span> function
converts a VEP tab-delimited output to a GVB input. This excludes
variants without the gene symbol.

``` r
library(GVB)
```

``` r
#path to individual VEP-tab files
sample_1_path <- system.file("exdata", "Sample_1_VEP-tab.txt", package = "GVB")
sample_2_path <- system.file("exdata", "Sample_2_VEP-tab.txt", package = "GVB")
sample_3_path <- system.file("exdata", "Sample_3_VEP-tab.txt", package = "GVB")
```

``` r
#convert these to GVB input format
sample_1 <- GVB_converting(sample_1_path)
```

    ## "Variants without a gene symbol were excluded. As a result, there are  690  variants remaining."

``` r
sample_2 <- GVB_converting(sample_2_path)
```

    ## "Variants without a gene symbol were excluded. As a result, there are  663  variants remaining."

``` r
sample_3 <- GVB_converting(sample_3_path)
```

    ## "Variants without a gene symbol were excluded. As a result, there are  702  variants remaining."

``` r
#converted VEP-tab object (GVB input)
str(sample_1)
```

    ## 'data.frame':    690 obs. of  11 variables:
    ##  $ Chromosome : int  22 22 22 22 22 22 22 22 22 22 ...
    ##  $ Start      : int  16449050 16487825 16871080 16914208 17006833 17072483 17073066 17103717 17264565 17265124 ...
    ##  $ Ref        : chr  "G" "-" "C" "A" ...
    ##  $ Alt        : chr  "A" "T" "A" "G" ...
    ##  $ Gene       : chr  "OR11H1" "YME1L1P1" "ABCD1P4" "SLC9B1P4" ...
    ##  $ Consequence: chr  "missense_variant" "splice_region_variant" "splice_polypyrimidine_tract_variant" "splice_region_variant" ...
    ##  $ CADD       : num  NA NA NA NA NA ...
    ##  $ SIFT       : num  0 NA NA NA NA 1 NA NA 0.02 1 ...
    ##  $ PP2        : num  0.972 NA NA NA NA 0 NA NA 0.094 0 ...
    ##  $ Genotype   : int  1 1 1 1 2 2 2 1 2 2 ...
    ##  $ Sample_ID  : chr  "S1" "S1" "S1" "S1" ...

<br>

### 5. Preprocessing the GVB input

<span style="background-color: #EDEDED;">GVB_preProcessing</span>
function preprocesses the GVB input and converts it to a list where each
element is a data frame for each gene before calculating the GVB score.
In this function, insertion and deletion variant (indel) can be excluded
or included.

``` r
#preprocess the GVB input while excluding indels
sample_1_pre <- GVB_preProcessing(sample_1, indel = FALSE)
sample_2_pre <- GVB_preProcessing(sample_2, indel = FALSE)
sample_3_pre <- GVB_preProcessing(sample_3, indel = FALSE)
```

<br>

### 6. Calculating GVB score using SIFT

<span style="background-color: #EDEDED;">GVB</span> function calculates
the GVB score and returns a data frame of the GVB scores. In this
function, you can set the *in silico* deleteriousness score, its
threshold, calculation method for the GVB score, etc.

``` r
#calculate the GVB score
sample_1_gvb <- GVB(sample_1_pre, score = "SIFT", threshold = 0.7, direction = "or_less", method = "gm")
sample_2_gvb <- GVB(sample_2_pre, score = "SIFT", threshold = 0.7, direction = "or_less", method = "gm")
sample_3_gvb <- GVB(sample_3_pre, score = "SIFT", threshold = 0.7, direction = "or_less", method = "gm")
```

``` r
#GVB score object
str(sample_1_gvb)
```

    ## 'data.frame':    282 obs. of  2 variables:
    ##  $ Gene: chr  "A4GALT" "ABCD1P4" "ABHD17AP5" "AC002472.1" ...
    ##  $ S1  : num  0.26 1 1 1 0 1 1 1 1 0.57 ...

<br>

### 7. merging multiple GVB data frames

<span style="background-color: #EDEDED;">GVB_merging</span> function
merges GVB data frames of multiple samples derived from
<span style="background-color: #EDEDED;">GVB</span> function.

``` r
#merges GVB data frames of multiple samples
merged_3samples_gvb <- GVB_merging(ls(pattern = "_gvb"))
```

``` r
#merged GVB data object
str(merged_3samples_gvb)
```

    ##  num [1:3, 1:372] 0.26 0.26 0.26 1 1 1 1 1 1 1 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:3] "S1" "S2" "S3"
    ##   ..$ : chr [1:372] "A4GALT" "ABCD1P4" "ABHD17AP5" "AC002472.1" ...

<br>

### 8. Visualize merged GVB data frame

<span style="background-color: #EDEDED;">GVB_histogram</span> function
plots the GVB scores for one specific gene as a histogram from merged
GVB data, which is the output of the
<span style="background-color: #EDEDED;">GVB</span> function. The red
vertical dashed line represents the mean GVB score of that gene across
all samples.

``` r
#plot histogram for CACNA1I gene from merged GVB data for 3 samples
GVB_histogram(data = merged_3samples_gvb, geneSymbol = "CACNA1I")
```

![1](https://github.com/mureal-bioinfo/GVB/assets/42491429/ded95b72-6e4d-447e-bf97-c0b5262cbac6)
``` r
#plot histogram for CACNA1I gene from merged GVB data for 100 samples
GVB_histogram(data = merged_gvb, geneSymbol = "CACNA1I")
```

![2](https://github.com/mureal-bioinfo/GVB/assets/42491429/545bd7f3-23ae-46bb-b00a-7cf6df0e7358)

### 9. Session Info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Korean_Korea.949  LC_CTYPE=Korean_Korea.949   
    ## [3] LC_MONETARY=Korean_Korea.949 LC_NUMERIC=C                
    ## [5] LC_TIME=Korean_Korea.949    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] GVB_1.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.14   knitr_1.43        magrittr_2.0.3    tidyselect_1.2.0 
    ##  [5] mnormt_2.1.1      lattice_0.20-45   R6_2.5.1          rlang_1.1.1      
    ##  [9] fastmap_1.1.0     fansi_1.0.4       highr_0.10        dplyr_1.1.2      
    ## [13] tools_4.1.2       grid_4.1.2        parallel_4.1.2    data.table_1.14.8
    ## [17] nlme_3.1-162      xfun_0.39         utf8_1.2.3        psych_2.3.9      
    ## [21] cli_3.6.0         withr_2.5.0       htmltools_0.5.4   yaml_2.3.7       
    ## [25] digest_0.6.31     tibble_3.2.1      lifecycle_1.0.3   vctrs_0.6.2      
    ## [29] glue_1.6.2        evaluate_0.21     rmarkdown_2.21    compiler_4.1.2   
    ## [33] pillar_1.9.0      generics_0.1.3    pkgconfig_2.0.3
