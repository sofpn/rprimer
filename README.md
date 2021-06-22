
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="./man/figures/rprimer.png" width="15%" style="float:right; padding:10px" style="display: block; margin: auto 0 auto auto;" />

<!-- badges: start -->

[![R build
status](https://github.com/sofpn/rprimer/workflows/R-CMD-check/badge.svg)](https://github.com/sofpn/rprimer/actions)
[![Codecov test
coverage](https://codecov.io/gh/sofpn/rprimer/branch/master/graph/badge.svg)](https://codecov.io/gh/sofpn/rprimer?branch=master)
<!-- badges: end -->

*This package is in development*

## Installation

rprimer can be installed from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("sofpn/rprimer")
```

``` r
library(rprimer)
```

## Overview

rprimer provides tools for designing degenerate DNA oligos for sequence
variable targets.

The design workflow consists of five functions:

  - `consensusProfile()`
  - `oligos()`
  - `assays()`
  - `checkMatch()`
  - `plotData()`

## Shiny application

The design workflow can be run through a Shiny application. It is loaded
by:

  - `runRprimerApp()`

The application is also available
[online](https://sofpn.shinyapps.io/rprimer).

## Workflow

### Import alignment

The first step is to import an alignment with target sequences of
interest. This is done by using `readDNAMultipleAlignment()` from the
Biostrings package.

The file “example\_alignment.txt” contains an alignment of 200 hepatitis
E virus sequences.

``` r
infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")

myAlignment <- Biostrings::readDNAMultipleAlignment(infile, format = "fasta")
```

### Step 1: `consensusProfile`

`consensusProfile()` takes a `Biostrings::DNAMultipleAlignment` as input
and returns all the information needed for the subsequent design
process.

``` r
myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)
```

Results (row 100-110):

| position |    a |    c | g |    t | other | gaps | majority | identity | iupac | entropy | coverage |
| -------: | ---: | ---: | -: | ---: | ----: | ---: | :------- | -------: | :---- | ------: | -------: |
|      100 | 0.00 | 1.00 | 0 | 0.00 |     0 |    0 | C        |     1.00 | C     |    0.00 |     1.00 |
|      101 | 0.10 | 0.00 | 0 | 0.90 |     0 |    0 | T        |     0.90 | W     |    0.50 |     1.00 |
|      102 | 1.00 | 0.00 | 0 | 0.00 |     0 |    0 | A        |     1.00 | A     |    0.00 |     1.00 |
|      103 | 0.00 | 1.00 | 0 | 0.00 |     0 |    0 | C        |     1.00 | C     |    0.00 |     1.00 |
|      104 | 0.01 | 0.08 | 0 | 0.91 |     0 |    0 | T        |     0.91 | Y     |    0.48 |     0.99 |
|      105 | 0.00 | 0.00 | 1 | 0.00 |     0 |    0 | G        |     1.00 | G     |    0.00 |     1.00 |
|      106 | 0.00 | 1.00 | 0 | 0.00 |     0 |    0 | C        |     1.00 | C     |    0.05 |     1.00 |
|      107 | 0.00 | 0.64 | 0 | 0.35 |     0 |    0 | C        |     0.64 | Y     |    0.98 |     0.99 |
|      108 | 1.00 | 0.00 | 0 | 0.00 |     0 |    0 | A        |     1.00 | A     |    0.05 |     1.00 |
|      109 | 0.00 | 0.00 | 0 | 1.00 |     0 |    0 | T        |     1.00 | T     |    0.00 |     1.00 |
|      110 | 0.00 | 0.08 | 0 | 0.92 |     0 |    0 | T        |     0.92 | Y     |    0.40 |     1.00 |

The results can be visualized with `plotData()`. You can either plot the
entire genome:

``` r
plotData(myConsensusProfile)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Step 2: `oligos`

The next step is to design oligos. You can either use the default
settings as below, or adjust them (see the package vignette or `?oligos`
for more information).

``` r
myOligos <- oligos(myConsensusProfile)
```

Results (first six rows):

| type   | fwd   | rev   | start | end | length | iupacSequence        | iupacSequenceRc      | identity | coverage | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | deltaGMean | deltaGRange | sequence   | sequenceRc | gcContent  | tm         | deltaG      | method    | score | roiStart | roiEnd |
| :----- | :---- | :---- | ----: | --: | -----: | :------------------- | :------------------- | -------: | -------: | ---------: | ------------: | -------------: | -----: | ------: | ---------: | ----------: | :--------- | :--------- | :--------- | :--------- | :---------- | :-------- | ----: | -------: | -----: |
| probe  | TRUE  | TRUE  |   138 | 157 |     20 | TCYGCCYTGGCGAATGCTGT | ACAGCATTCGCCARGGCRGA |     0.96 |        1 |          4 |          0.60 |           0.10 |  63.17 |    4.33 |    \-21.61 |        2.08 | TCCGCCCT…. | ACAGCATT…. | 0.65, 0….. | 65.33623…. | \-22.6538…. | ambiguous |     5 |        1 |   7997 |
| probe  | FALSE | TRUE  |   141 | 160 |     20 | GCCYTGGCGAATGCTGTGGT | ACCACAGCATTCGCCARGGC |     0.98 |        1 |          2 |          0.62 |           0.05 |  63.18 |    1.84 |    \-21.73 |        0.83 | GCCCTGGC…. | ACCACAGC…. | 0.65, 0.6  | 64.10586…. | \-22.1475…. | ambiguous |     4 |        1 |   7997 |
| primer | TRUE  | FALSE |   142 | 160 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.98 |        1 |          2 |          0.61 |           0.05 |  61.48 |    1.95 |    \-19.84 |        0.83 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 62.45335…. | \-20.2570…. | ambiguous |     4 |        1 |   7997 |
| probe  | TRUE  | TRUE  |   142 | 160 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.98 |        1 |          2 |          0.61 |           0.05 |  60.45 |    1.94 |    \-19.84 |        0.83 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 61.41686…. | \-20.2570…. | ambiguous |     4 |        1 |   7997 |
| primer | TRUE  | FALSE |   143 | 160 |     18 | CYTGGCGAATGCTGTGGT   | ACCACAGCATTCGCCARG   |     0.98 |        1 |          2 |          0.58 |           0.06 |  59.19 |    2.04 |    \-18.35 |        0.83 | CCTGGCGA…. | ACCACAGC…. | 0.611111…. | 60.20693…. | \-18.7709…. | ambiguous |     4 |        1 |   7997 |
| probe  | TRUE  | TRUE  |   143 | 160 |     18 | CYTGGCGAATGCTGTGGT   | ACCACAGCATTCGCCARG   |     0.98 |        1 |          2 |          0.58 |           0.06 |  58.11 |    2.03 |    \-18.35 |        0.83 | CCTGGCGA…. | ACCACAGC…. | 0.611111…. | 59.12650…. | \-18.7709…. | ambiguous |     4 |        1 |   7997 |

The results can be visualized as a dashboard, using `plotData()`:

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and combines them
with probes, if probes are present in the input dataset. You can either
use the default settings as below, or adjust the design constraints (see
the package vignette or `?assays` for more information).

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

| start |  end | length | totalDegeneracy | score | startFwd | endFwd | lengthFwd | iupacSequenceFwd     | identityFwd | coverageFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | deltaGMeanFwd | deltaGRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | deltaGFwd   | methodFwd | startRev | endRev | lengthRev | iupacSequenceRev     | identityRev | coverageRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | deltaGMeanRev | deltaGRangeRev | sequenceRev | gcContentRev | tmRev      | deltaGRev   | methodRev | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | identityPr | coveragePr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | deltaGMeanPr | deltaGRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | deltaGPr    | methodPr  | roiStart | roiEnd |
| ----: | ---: | -----: | --------------: | ----: | -------: | -----: | --------: | :------------------- | ----------: | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | ------------: | -------------: | :---------- | :----------- | :--------- | :---------- | :-------- | -------: | -----: | --------: | :------------------- | ----------: | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | ------------: | -------------: | :---------- | :----------- | :--------- | :---------- | :-------- | :----- | :------ | ------: | ----: | -------: | :--------------------- | :--------------------- | ---------: | ---------: | -----------: | --------------: | ---------------: | -------: | --------: | -----------: | ------------: | :--------- | :----------- | :---------- | :--------- | :---------- | :-------- | -------: | -----: |
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | FALSE   |    5946 |  5963 |       18 | CMGGGTTGATTCTCAGCC     | GGCTGAGAATCAACCCKG     |       0.98 |          1 |            2 |            0.58 |             0.06 |    55.14 |      2.79 |      \-17.06 |          1.25 | CAGGGTTG…. | GGCTGAGA….   | 0.555555….  | 53.74554…. | \-16.4324…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    12 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | FALSE   |    5946 |  5964 |       19 | CMGGGTTGATTCTCAGCCC    | GGGCTGAGAATCAACCCKG    |       0.98 |          1 |            2 |            0.61 |             0.05 |    57.63 |      2.64 |      \-18.54 |          1.25 | CAGGGTTG…. | GGGCTGAG….   | 0.578947….  | 56.30713…. | \-17.9185…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | TRUE    |    5946 |  5965 |       20 | CMGGGTTGATTCTCAGCCCT   | AGGGCTGAGAATCAACCCKG   |       0.98 |          1 |            2 |            0.58 |             0.05 |    58.87 |      2.55 |      \-19.43 |          1.25 | CAGGGTTG…. | AGGGCTGA….   | 0.55, 0.6   | 57.59836…. | \-18.8035…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    10 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | TRUE    |    5946 |  5966 |       21 | CMGGGTTGATTCTCAGCCCTT  | AAGGGCTGAGAATCAACCCKG  |       0.98 |          1 |            2 |            0.55 |             0.05 |    59.21 |      2.43 |      \-20.08 |          1.25 | CAGGGTTG…. | AAGGGCTG….   | 0.523809….  | 57.99472…. | \-19.4553…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | FALSE   |    5946 |  5967 |       22 | CMGGGTTGATTCTCAGCCCTTC | GAAGGGCTGAGAATCAACCCKG |       0.98 |          1 |            2 |            0.57 |             0.05 |    59.91 |      2.28 |      \-21.11 |          1.25 | CAGGGTTG…. | GAAGGGCT….   | 0.545454….  | 58.77533…. | \-20.4881…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    10 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |       \-20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-20.2350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |       \-15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-15.0078…. | ambiguous | TRUE   | FALSE   |    5947 |  5964 |       18 | MGGGTTGATTCTCAGCCC     | GGGCTGAGAATCAACCCK     |       0.98 |          1 |            2 |            0.58 |             0.06 |    55.71 |      1.65 |      \-17.21 |          0.94 | AGGGTTGA…. | GGGCTGAG….   | 0.555555….  | 54.88806…. | \-16.7409…. | ambiguous |        1 |   7997 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

### Check match

`checkMatch()` shows the proportion and names of the target sequences in
the input alignment that match with the generated oligos or assays:

``` r
## Randomly select six oligos to illustrate an example 
selection <- sample(seq_len(nrow(myOligos)), size = 6)

matchTableOligos <- checkMatch(myOligos[selection, ], target = myAlignment)
```

Results:

| iupacSequence          | perfectMatch | idPerfectMatch | oneMismatch | idOneMismatch | twoMismatches | idTwoMismatches | threeMismatches | idThreeMismatches | fourOrMoreMismatches | idFourOrMoreMismatches |
| :--------------------- | -----------: | :------------- | ----------: | :------------ | ------------: | :-------------- | --------------: | :---------------- | -------------------: | :--------------------- |
| TTCTGGGGTGACMGGGTTGATT |         0.94 | M73218.1….     |        0.04 | AB222183….    |             0 | JQ953665.1      |               0 |                   |                 0.01 | MN614142….             |
| GATTCTCAGCCCTTCGCMMTC  |         0.94 | M73218.1….     |        0.05 | AB222183….    |             0 |                 |               0 | JF443725.1        |                 0.01 | MN614142….             |
| GTYGTCTCRGCCAATGGCGA   |         0.96 | M73218.1….     |        0.04 | AB630970….    |             0 |                 |               0 |                   |                 0.00 |                        |
| GGTGACMGGGTTGATTCT     |         0.93 | M73218.1….     |        0.06 | AB222183….    |             0 | JQ953665.1      |               0 |                   |                 0.01 | MN614142….             |
| GGGGTGACMGGGTTGATTCTC  |         0.93 | M73218.1….     |        0.06 | AB222183….    |             0 | JQ953665.1      |               0 |                   |                 0.01 | MN614142….             |
| TGGGGTGACMGGGTTGATT    |         0.94 | M73218.1….     |        0.04 | AB222183….    |             0 | JQ953665.1      |               0 |                   |                 0.01 | MN614142….             |

The match table can be visualized using `plotData()`:

``` r
plotData(matchTableOligos)
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" />

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.

## Citation

To cite this package, please use: `citation("rprimer")`.
