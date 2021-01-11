
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> <!-- badges: end -->

## TODO

  - Assays=beautify, document, generate extdata
  - Oligos=describe design process
  - Plots=fix assay method
  - Classes=document, set validity
  - All=test
  - Vignette

## Installation

You can install rprimer from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("sofpn/rprimer")
```

``` r
# library(rprimer)
devtools::load_all(".")
```

## Overview

rprimer provides tools for designing broadly reactive primers, probes
and (RT)-(q/dd)PCR assays from a multiple DNA sequence alignment.

The design process is built on three functions:

  - `consensusProfile()`
  - `oligos()`
  - `assays()`

## Workflow

### Import alignment

The first step is to import an alignment with target sequences of
interest and, if preferred, mask positions with e.g. high gap frequency.
`readDNAMultipleAlignment()` and `maskGaps()` from Biostrings do the
work for this part. The file “example\_alignment.txt” contains an
alignment of 100 hepatitis E virus sequences.

``` r
infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")

myAlignment <- Biostrings::readDNAMultipleAlignment(infile, format = "fasta")
myAlignment <- Biostrings::maskGaps(myAlignment, 
                                    min.fraction = 0.5, 
                                    min.block.width = 1) 
```

### Step 1: `consensusProfile`

`consensusProfile()` takes a `Biostrings::DNAMultipleAlignment` as input
and returns all the information needed for the subsequent design
process.

``` r
myConsensusProfile <- consensusProfile(myAlignment, iupacThreshold = 0.05)
```

Results (first six rows):

| position |    a |    c |    g |    t | other | gaps | majority | identity | iupac | entropy |
| -------: | ---: | ---: | ---: | ---: | ----: | ---: | :------- | -------: | :---- | ------: |
|        1 | 0.00 | 0.00 | 0.59 | 0.00 |     0 | 0.41 | G        |     1.00 | G     |    0.00 |
|        2 | 0.00 | 0.00 | 0.71 | 0.00 |     0 | 0.29 | G        |     1.00 | G     |    0.00 |
|        3 | 0.00 | 0.71 | 0.00 | 0.00 |     0 | 0.29 | C        |     1.00 | C     |    0.00 |
|        4 | 0.71 | 0.00 | 0.00 | 0.00 |     0 | 0.29 | A        |     1.00 | A     |    0.00 |
|        5 | 0.00 | 0.00 | 0.70 | 0.01 |     0 | 0.29 | G        |     0.99 | G     |    0.11 |
|        6 | 0.71 | 0.00 | 0.00 | 0.00 |     0 | 0.29 | A        |     1.00 | A     |    0.00 |

The results can be visualized with `plotData()`. You can either plot the
entire genome:

``` r
plotData(myConsensusProfile)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

Or zoom into a specific region of interest:

``` r
roi <- myConsensusProfile[myConsensusProfile$position >= 5000 & myConsensusProfile$position <= 5050, ]
plotData(roi)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

The nucleotide distribution can be shown by specifying `type =
"nucleotide`:

``` r
plotData(roi, type = "nucleotide")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Step 2: `oligos`

The next step is to design oligos. You can either use the default
settings:

``` r
myOligos <- oligos(myConsensusProfile)
```

Or adjust the constraints (see the package vignette or
`?rprimer::oligos` for more information), for example:

``` r
# myOligos2 <- oligos(myConsensusProfile, 
#                     maxDegeneracyPrimer = 16,
#                     tmRangePrimer = c(58, 60),
#                     minEndIdentityPrimer =  0.99,
#                     gcClampPrimer = FALSE,
#                     avoidThreeEndRunsPrimer = FALSE,
#                     probe = FALSE)
```

Results (first six rows):

| type   | fwd   | rev  | start | end | length | iupacSequence        | iupacSequenceRc      | identity | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | sequence   | sequenceRc | gcContent  | tm         | roiStart | roiEnd |
| :----- | :---- | :--- | ----: | --: | -----: | :------------------- | :------------------- | -------: | ---------: | ------------: | -------------: | -----: | ------: | :--------- | :--------- | :--------- | :--------- | -------: | -----: |
| primer | FALSE | TRUE |    27 |  45 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |     0.97 |          4 |          0.53 |           0.11 |  57.64 |    5.01 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 60.12783…. |        1 |   7208 |
| primer | FALSE | TRUE |    27 |  46 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |     0.97 |          4 |          0.50 |           0.10 |  57.93 |    4.99 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 60.42434…. |        1 |   7208 |
| probe  | TRUE  | TRUE |    27 |  45 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |     0.97 |          4 |          0.53 |           0.11 |  56.58 |    5.01 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 59.07612…. |        1 |   7208 |
| probe  | TRUE  | TRUE |    27 |  46 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |     0.97 |          4 |          0.50 |           0.10 |  56.92 |    4.98 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 59.40578…. |        1 |   7208 |
| primer | FALSE | TRUE |    28 |  46 |     19 | TGGAGGCCCAYCAGTTYAT  | ATRAACTGRTGGGCCTCCA  |     0.97 |          4 |          0.53 |           0.11 |  57.55 |    5.16 | TGGAGGCC…. | ATGAACTG…. | 0.578947…. | 60.12783…. |        1 |   7208 |
| probe  | TRUE  | TRUE |    28 |  46 |     19 | TGGAGGCCCAYCAGTTYAT  | ATRAACTGRTGGGCCTCCA  |     0.97 |          4 |          0.53 |           0.11 |  56.51 |    5.15 | TGGAGGCC…. | ATGAACTG…. | 0.578947…. | 59.07612…. |        1 |   7208 |

The results can be visualized using `plotData()`.

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and combines them
with probes, if `probe = TRUE` was used in the oligo design step.

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

The assays can be visualized using `plotData()`:

``` r
#plotData(myAssays)
```

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
