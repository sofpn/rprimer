
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> [![R build
status](https://github.com/sofpn/rprimer/workflows/R-CMD-check/badge.svg)](https://github.com/sofpn/rprimer/actions)
<!-- badges: end -->

### Package overview

rprimer provides functions for designing (RT)-(q/dd)PCR assays from
multiple DNA sequence alignments. The design process is built on three
functions:

  - `getConsensusProfile()`
  - `getOligos()`
  - `getAssays()`

### Installation

You can install rprimer from [GitHub](https://github.com/) with:

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

### Demonstration

Below, I demonstrate how to use rprimer by designing an RT-(q/d)PCR
assay for detection of hepatitis E virus, which is a highly variable RNA
virus.

### To start

The first step is to import an alignment with target sequences of
interest and, if preferred, mask positions with high gap frequency.
`readDNAMultipleAlignment()` and `maskGaps()` from Biostrings can be
used for this part.

The file “example\_alignment.txt” is provided with the rprimer package
and contains 100 hepatitis E virus sequences.

``` r
infile <- system.file('extdata', 'example_alignment.txt', package = 'rprimer')

myAlignment <- infile %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1) 
## Mask positions with at least 50 % gaps 
```

### Step 1: `getConsensusProfile()`

`getConsensusProfile()` takes a `Biostrings::DNAMultipleAlignment`
object as input and retrieves the information needed for the design
process.

All nucleotides with a frequency higher than the `iupacThreshold` will
be included in the IUPAC consensus base. Hence, oligos generated from an
object with a low `iupacThreshold` will contain more degenerate bases
than those generated with a high.

``` r
 myConsensusProfile <- getConsensusProfile(myAlignment, iupacThreshold = 0.05)
```

The ougput can be visualised using `plotConsensusProfile()`.

``` r
plotConsensusProfile(myConsensusProfile)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Moreover, the nucleotide distribution can be visualized using
`plotNucleotides()`. The `rc` option regulates whether the sequence
should be displayed as a reverse complement or not.

``` r
plotNucleotides(myConsensusProfile, from = 1, to = 30, rc = FALSE) ## Plot the first 30 bases 
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" style="display: block; margin: auto;" />

### Step 2: `getOligos()`

`getOligos()` searches for oligos from the following constraints:

  - `maxGapFrequency` Maximum gap frequency, defaults to `0.1`.
  - `length` Oligo length, defaults to `18:22`.
  - `maxDegeneracy` Maximum number of degenerate variants of each oligo,
    defaults to `4`.
  - `gcClamp` If oligos must have a GC-clamp (recommended for primers),
    defaults to `TRUE`.
  - `avoid3endRuns` If oligos with more than two runs of the same
    nucleotide at the terminal 3’ end should be avodied (recommended for
    primers), defaults to `TRUE`.
  - `avoid5endG` If oligos with a G at the terminal 5’ end should be
    avoided (recommended for probes), defaults to `FALSE`.
  - `minEndIdentity` Optional. Minimum allowed identity at the 3’ end
    (i.e. the last five bases). E.g., if set to `1`, only oligos with
    complete target conservation at the 3’ end will be considered.
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

An error message will return if no oligos are found.

Below, I want to design both primers and probes. I use somewhat
different settings for the two oligo types.

``` r
myPrimers <- getOligos(myConsensusProfile,
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

myProbes <-  getOligos(myConsensusProfile,
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
them with probes (if selected).

Assays are designed from the following constraints:

  - `length` Amplicon length, defaults to `65:120`.
  - `maxTmDifferencePrimers` Maximum Tm difference between the two
    primers (absolute value, calculated for majority primers), defaults
    to `2`.
  - `tmDifferenceProbes` Acceptable Tm difference between the primers
    (average Tm of the primer pair) and probe, defaults to `c(0, 20)`.
    The Tm-difference is calculated by subtracting the Tm of the probe
    with the average Tm of the (majority) primer pair. Thus, a negative
    Tm-difference means that the Tm of the probe is lower than the
    average Tm of the primer pair.

An error message will return if no assays are found.

``` r
myAssays <- getAssays(primers = myPrimers, probes = myProbes)
#> New names:
#> * majorityRev -> majorityRev...17
#> * iupacRev -> iupacRev...23
#> * majorityRev -> majorityRev...27
#> * iupacRev -> iupacRev...30
```

### Further notes

degeneracy… majority gap masking pitfall

*To do:*

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
#>  [1] tidyselect_1.1.0  xfun_0.17         remotes_2.2.0     reshape2_1.4.4   
#>  [5] purrr_0.3.4       colorspace_1.4-1  vctrs_0.3.4       generics_0.0.2   
#>  [9] usethis_1.6.1     htmltools_0.5.0   yaml_2.2.1        rlang_0.4.7      
#> [13] pkgbuild_1.1.0    pillar_1.4.6      glue_1.4.2        withr_2.2.0      
#> [17] sessioninfo_1.1.1 plyr_1.8.6        lifecycle_0.2.0   stringr_1.4.0    
#> [21] zlibbioc_1.35.0   munsell_0.5.0     gtable_0.3.0      devtools_2.3.1   
#> [25] memoise_1.1.0     evaluate_0.14     labeling_0.3      knitr_1.29       
#> [29] callr_3.4.4       ps_1.3.4          fansi_0.4.1       Rcpp_1.0.5       
#> [33] backports_1.1.9   scales_1.1.1      desc_1.2.0        pkgload_1.1.0    
#> [37] farver_2.0.3      fs_1.5.0          ggplot2_3.3.2     digest_0.6.25    
#> [41] stringi_1.5.3     processx_3.4.4    dplyr_1.0.2       rprojroot_1.3-2  
#> [45] grid_4.0.2        cli_2.0.2         tools_4.0.2       patchwork_1.0.1  
#> [49] tibble_3.0.3      crayon_1.3.4      pkgconfig_2.0.3   ellipsis_0.3.1   
#> [53] prettyunits_1.1.1 assertthat_0.2.1  rmarkdown_2.3     rstudioapi_0.11  
#> [57] R6_2.4.1          compiler_4.0.2
```
