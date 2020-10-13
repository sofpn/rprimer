
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> [![R build
status](https://github.com/sofpn/rprimer/workflows/R-CMD-check/badge.svg)](https://github.com/sofpn/rprimer/actions)
<!-- badges: end -->

### Introduction

The purpose of rprimer is to design (RT)-(q/dd)PCR assays from multiple
DNA sequence alignments.

In this document, I demonstrate how to use the package by designing an
RT-(q/dd)PCR assay for detection of hepatitis E virus, which is a highly
variable RNA virus.

### Package overview

The assay design workflow is based on the following functions:

  - `getAlignmentProfile()`
  - `getAlignmentProperties()`
  - `getOligos()`
  - `getAssays()`

### Installation

You can install the development version of rprimer from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sofpn/rprimer")
```

Initial setup for the code in this document:

``` r
# library(rprimer)
devtools::load_all(".")
library(magrittr) ## Required for the pipe operator 
library(Biostrings) ## Required to import alignments 
```

### To start

The first step is to import an alignment with target sequences of
interest and mask positions with high gap frequency. Previously existing
functionality from the Biostrings-package can be used for this part.

The file “example\_alignment.txt” is provided with the rprimer package
and contains 100 hepatitis E virus sequences. I use it in the example
below:

``` r
infile <- system.file('extdata', 'example_alignment.txt', package = 'rprimer')

myAlignment <- infile %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1) 
## Mask positions with at least 50 % gaps 
```

### Step 1: `getAlignmentProfile()`

`getAlignmentProfile()` takes a `Biostrings::DNAMultipleAlignment`
object as input and returns an `RprimerProfile` object, which contains a
numeric matrix that holds the proportion of each nucleotide at each
position in the alignment. The function is a wrapper around
`Biostrings::consensusMatrix()`.

``` r
myAlignmentProfile <- getAlignmentProfile(myAlignment)
rpGetData(myAlignmentProfile)[ , 1:10] ## View the first 10 bases 
#>          1    2    3    4    5    6    7    8    9   10
#> A     0.00 0.00 0.00 0.71 0.00 0.71 0.00 0.00 0.75 0.02
#> C     0.00 0.00 0.71 0.00 0.00 0.00 0.72 0.76 0.01 0.75
#> G     0.59 0.71 0.00 0.00 0.70 0.00 0.00 0.00 0.00 0.00
#> T     0.00 0.00 0.00 0.00 0.01 0.00 0.00 0.00 0.00 0.03
#> -     0.41 0.29 0.29 0.29 0.29 0.29 0.28 0.24 0.24 0.20
#> other 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
## rpGetData is a getter method for all classes within the package 
```

The object can be visualized using `rpPlot()`. The `rc` option regulates
whether it should be displayed as a reverse complement or not.

``` r
rpPlot(myAlignmentProfile[, 1:30], rc = FALSE) ## Plot the first 30 bases 
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" style="display: block; margin: auto;" />

### Step 2: `getAlignmentProperties()`

`getAlignmentProperties()` takes an `RprimerProfile`-object as input and
returns an `RprimerProperties`-object with information about majority
and IUPAC consensus base at each position, together with gap frequency,
nucleotide identity and Shannon entropy. All nucleotides that occurs
with a frequency higher than the `iupacThreshold` will be included in
the IUPAC consensus base. Hence, for downstream oligo design, a low
`iupacThreshold` will generate oligos with more degenerate bases than a
high `iupacThreshold`.

``` r
myAlignmentProperties <- getAlignmentProperties(myAlignmentProfile, iupacThreshold = 0.05)
head(myAlignmentProperties)
#> # A tibble: 6 x 6
#>   Position Majority IUPAC  Gaps Identity Entropy
#>      <int> <chr>    <chr> <dbl>    <dbl>   <dbl>
#> 1        1 G        G     0.41      1       0   
#> 2        2 G        G     0.290     1       0   
#> 3        3 C        C     0.290     1       0   
#> 4        4 A        A     0.290     1       0   
#> 5        5 G        G     0.290     0.99    0.11
#> 6        6 A        A     0.290     1       0
```

The object can be visualised with `rpPlot()`.

### Step 3: `getOligos()`

`getOligos()` takes an `RprimerProperties` object as input and searches
for oligos based on the following constraints:

  - `maxGapFrequency` Maximum gap frequency, defaults to `0.1`.
  - `length` Oligo length, defaults to `18:22`.
  - `maxDegeneracy` Maximum number of degenerate variants of each oligo,
    defaults to `4`.
  - `gcClamp` If oligos must have a GC-clamp to be considered as valid
    (recommended for primers), defaults to `TRUE`.
  - `avoid3endRuns` If oligos with more than two runs of the same
    nucleotide at the terminal 3’ end should be excluded (recommended
    for primers), defaults to `TRUE`.
  - `avoid5endG` If oligos with a G at the terminal 5’ end should be
    avoided (recommended for probes), defaults to `FALSE`.
  - `minEndIdentity` Optional. Minimum allowed identity at the 3’ end
    (i.e. the last five bases). E.g., if set to `1`, only oligos with
    complete target conservation at the 3’ end will be considered, and
    hence no wobble bases will be present at the 3’ ends.
  - `gcRange` GC-content-range, defaults to `c(0.45, 0.55)`.
  - `tmRange` Melting temperature (Tm) range, defaults to `c(50, 65)`.
    Tm is calculated using the nearest-neighbor method. See
    `?rprimer::getOligos` for a detailed description and references.
  - `concOligo` Oligo concentration (for Tm calculation), defaults to
    `5e-07` M (500 nM)
  - `concNa` Sodium ion concentration (for Tm calculation), defaults to
    `0.05` M (50 mM).
  - `showAllVariants` If sequence, GC-content and Tm should be presented
    for all variants of each oligo (in case of degenerate bases).`TRUE`
    (slower) or `FALSE` (faster), defaults to `TRUE`.

In addition, `get_oligos()` avoids:

  - Oligos with more than than three consecutive runs of the same
    dinucleotide (e.g. “TATATATA”)
  - Oligos with more than four consecutive runs of the same nucleotide
    (e.g. “AAAAA”)
  - Oligos that are duplicated (to prevent binding at several places on
    the genome)

Below, I want to design both primers and probes. I use somewhat
different settings for the two oligo types.

An error message will return if no oligos are found.

``` r
myPrimers <- getOligos(myAlignmentProperties,
                       length = 18:22,
                       maxGapFrequency = 0.05,
                       maxDegeneracy = 4,
                       gcClamp = TRUE,
                       avoid3EndRuns = TRUE,
                       avoid5EndG = FALSE,
                       minEndIdentity = 0.99,
                       gcRange = c(0.40, 0.60),
                       tmRange = c(50, 65),
                       showAllVariants = TRUE)

myProbes <-  getOligos(myAlignmentProperties,
                       length = 16:22,
                       maxGapFrequency = 0.05,
                       maxDegeneracy = 4,
                       gcClamp = FALSE,
                       avoid3EndRuns = FALSE,
                       avoid5EndG = TRUE,
                       minEndIdentity = NULL,
                       gcRange = c(0.40, 0.60),
                       tmRange = c(50, 75),
                       showAllVariants = TRUE)
```

### Step 4: `getAssays()`

`getAssays()` finds pairs of forward and reverse primers and combines
them with probes (if selected). It takes `RprimerOligo` objects as input
and returns an `RprimerAssay` object.

Assays are designed from the following constraints:

  - `length` Amplicon length, defaults to `65:120`.
  - `maxTmDifferencePrimers` Maximum Tm difference between the two
    primers (absolute value), defaults to `2`.
  - `tmDifferenceProbes` Acceptable Tm difference between the primers
    (average Tm of the primer pair) and probe, defaults to `c(0, 20)`.

Candidate assays are displayed in a tibble. An error message will return
if no assays are found.

``` r
myAssays <- getAssays(primers = myPrimers, 
                      probes = myProbes)
```

### Further notes

degeneracy…

### Session info

``` r
sessionInfo()
#> R version 4.0.2 (2020-06-22)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 18362)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Swedish_Sweden.1252  LC_CTYPE=Swedish_Sweden.1252   
#> [3] LC_MONETARY=Swedish_Sweden.1252 LC_NUMERIC=C                   
#> [5] LC_TIME=Swedish_Sweden.1252    
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#> [1] Biostrings_2.57.2   XVector_0.29.3      IRanges_2.23.10    
#> [4] S4Vectors_0.27.12   BiocGenerics_0.35.4 magrittr_1.5       
#> [7] rprimer_0.99.0      testthat_2.3.2     
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.5                  lattice_0.20-41            
#>  [3] prettyunits_1.1.1           ps_1.3.4                   
#>  [5] utf8_1.1.4                  assertthat_0.2.1           
#>  [7] rprojroot_1.3-2             digest_0.6.25              
#>  [9] R6_2.4.1                    GenomeInfoDb_1.25.11       
#> [11] plyr_1.8.6                  backports_1.1.9            
#> [13] evaluate_0.14               ggplot2_3.3.2              
#> [15] pillar_1.4.6                zlibbioc_1.35.0            
#> [17] rlang_0.4.7                 rstudioapi_0.11            
#> [19] callr_3.4.4                 Matrix_1.2-18              
#> [21] rmarkdown_2.3               labeling_0.3               
#> [23] desc_1.2.0                  devtools_2.3.1             
#> [25] stringr_1.4.0               RCurl_1.98-1.2             
#> [27] munsell_0.5.0               DelayedArray_0.15.8        
#> [29] compiler_4.0.2              xfun_0.17                  
#> [31] pkgconfig_2.0.3             pkgbuild_1.1.0             
#> [33] htmltools_0.5.0             tidyselect_1.1.0           
#> [35] SummarizedExperiment_1.19.6 tibble_3.0.3               
#> [37] GenomeInfoDbData_1.2.3      matrixStats_0.56.0         
#> [39] fansi_0.4.1                 crayon_1.3.4               
#> [41] dplyr_1.0.2                 withr_2.2.0                
#> [43] bitops_1.0-6                grid_4.0.2                 
#> [45] gtable_0.3.0                lifecycle_0.2.0            
#> [47] scales_1.1.1                cli_2.0.2                  
#> [49] stringi_1.5.3               farver_2.0.3               
#> [51] reshape2_1.4.4              fs_1.5.0                   
#> [53] remotes_2.2.0               ellipsis_0.3.1             
#> [55] generics_0.0.2              vctrs_0.3.4                
#> [57] tools_4.0.2                 Biobase_2.49.1             
#> [59] glue_1.4.2                  purrr_0.3.4                
#> [61] processx_3.4.4              pkgload_1.1.0              
#> [63] yaml_2.2.1                  colorspace_1.4-1           
#> [65] GenomicRanges_1.41.6        sessioninfo_1.1.1          
#> [67] memoise_1.1.0               knitr_1.29                 
#> [69] usethis_1.6.1
```
