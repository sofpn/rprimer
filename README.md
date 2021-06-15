
<!-- README.md is generated from README.Rmd. Please edit that file -->

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

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Or zoom into a specific region of interest:

``` r
## Select position 5000 to 5500 
roi <- myConsensusProfile[myConsensusProfile$position >= 5000 & myConsensusProfile$position <= 5500, ]

plotData(roi)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

The nucleotide distribution can be shown by specifying `type =
"nucleotide`:

``` r
## Select position 150 to 170
selection <- myConsensusProfile[myConsensusProfile$position >= 150 & myConsensusProfile$position <= 170, ]

plotData(selection, type = "nucleotide")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

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
| probe  | TRUE  | TRUE  |   138 | 157 |     20 | TCYGCCYTGGCGAATGCTGT | ACAGCATTCGCCARGGCRGA |     0.96 |        1 |          4 |          0.60 |           0.10 |  63.17 |    4.33 |   \-216.14 |       20.80 | TCCGCCCT…. | ACAGCATT…. | 0.65, 0….. | 65.33623…. | \-226.538…. | ambiguous |     5 |        1 |   7997 |
| probe  | FALSE | TRUE  |   141 | 160 |     20 | GCCYTGGCGAATGCTGTGGT | ACCACAGCATTCGCCARGGC |     0.98 |        1 |          2 |          0.62 |           0.05 |  63.18 |    1.84 |   \-217.30 |        8.34 | GCCCTGGC…. | ACCACAGC…. | 0.65, 0.6  | 64.10586…. | \-221.475…. | ambiguous |     4 |        1 |   7997 |
| primer | TRUE  | FALSE |   142 | 160 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.98 |        1 |          2 |          0.61 |           0.05 |  61.48 |    1.95 |   \-198.40 |        8.34 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 62.45335…. | \-202.570…. | ambiguous |     4 |        1 |   7997 |
| probe  | TRUE  | TRUE  |   142 | 160 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.98 |        1 |          2 |          0.61 |           0.05 |  60.45 |    1.94 |   \-198.40 |        8.34 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 61.41686…. | \-202.570…. | ambiguous |     4 |        1 |   7997 |
| primer | TRUE  | FALSE |   143 | 160 |     18 | CYTGGCGAATGCTGTGGT   | ACCACAGCATTCGCCARG   |     0.98 |        1 |          2 |          0.58 |           0.06 |  59.19 |    2.04 |   \-183.54 |        8.34 | CCTGGCGA…. | ACCACAGC…. | 0.611111…. | 60.20693…. | \-187.709…. | ambiguous |     4 |        1 |   7997 |
| probe  | TRUE  | TRUE  |   143 | 160 |     18 | CYTGGCGAATGCTGTGGT   | ACCACAGCATTCGCCARG   |     0.98 |        1 |          2 |          0.58 |           0.06 |  58.11 |    2.03 |   \-183.54 |        8.34 | CCTGGCGA…. | ACCACAGC…. | 0.611111…. | 59.12650…. | \-187.709…. | ambiguous |     4 |        1 |   7997 |

The results can be visualized as a dashboard, using `plotData()`:

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

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
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | FALSE   |    5946 |  5963 |       18 | CMGGGTTGATTCTCAGCC     | GGCTGAGAATCAACCCKG     |       0.98 |          1 |            2 |            0.58 |             0.06 |    55.14 |      2.79 |     \-170.55 |         12.45 | CAGGGTTG…. | GGCTGAGA….   | 0.555555….  | 53.74554…. | \-164.324…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    12 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | FALSE   |    5946 |  5964 |       19 | CMGGGTTGATTCTCAGCCC    | GGGCTGAGAATCAACCCKG    |       0.98 |          1 |            2 |            0.61 |             0.05 |    57.63 |      2.64 |     \-185.41 |         12.45 | CAGGGTTG…. | GGGCTGAG….   | 0.578947….  | 56.30713…. | \-179.185…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | TRUE    |    5946 |  5965 |       20 | CMGGGTTGATTCTCAGCCCT   | AGGGCTGAGAATCAACCCKG   |       0.98 |          1 |            2 |            0.58 |             0.05 |    58.87 |      2.55 |     \-194.26 |         12.45 | CAGGGTTG…. | AGGGCTGA….   | 0.55, 0.6   | 57.59836…. | \-188.035…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    10 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | TRUE    |    5946 |  5966 |       21 | CMGGGTTGATTCTCAGCCCTT  | AAGGGCTGAGAATCAACCCKG  |       0.98 |          1 |            2 |            0.55 |             0.05 |    59.21 |      2.43 |     \-200.78 |         12.45 | CAGGGTTG…. | AAGGGCTG….   | 0.523809….  | 57.99472…. | \-194.553…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    11 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | FALSE   |    5946 |  5967 |       22 | CMGGGTTGATTCTCAGCCCTTC | GAAGGGCTGAGAATCAACCCKG |       0.98 |          1 |            2 |            0.57 |             0.05 |    59.91 |      2.28 |     \-211.11 |         12.45 | CAGGGTTG…. | GAAGGGCT….   | 0.545454….  | 58.77533…. | \-204.881…. | ambiguous |        1 |   7997 |
|  5926 | 5994 |     69 |               6 |    10 |     5926 |   5945 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |      \-208.58 |          12.45 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | \-202.350…. | ambiguous |     5975 |   5994 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |      \-152.68 |           5.21 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | \-150.078…. | ambiguous | TRUE   | FALSE   |    5947 |  5964 |       18 | MGGGTTGATTCTCAGCCC     | GGGCTGAGAATCAACCCK     |       0.98 |          1 |            2 |            0.58 |             0.06 |    55.71 |      1.65 |     \-172.09 |          9.37 | AGGGTTGA…. | GGGCTGAG….   | 0.555555….  | 54.88806…. | \-167.409…. | ambiguous |        1 |   7997 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

### Check match

`checkMatch()` shows the proportion and names of the target sequences in
the input alignment that match with the generated oligos or assays:

``` r
## Randomly select six oligos to illustrate an example 
selection <- sample(seq_len(nrow(myOligos)), size = 6)

matchTableOligos <- checkMatch(myOligos[selection, ], target = myAlignment)
```

Results:

| iupacSequence          | perfectMatch | idPerfectMatch | oneMismatch | idOneMismatch | twoMismatches | idTwoMismatches | threeMismatches | idThreeMismatches | fourOrMoreMismatches | idFourOrMoreMmismatches |
| :--------------------- | -----------: | :------------- | ----------: | :------------ | ------------: | :-------------- | --------------: | :---------------- | -------------------: | :---------------------- |
| CCCCTATWTTCATCCAACCAA  |         0.96 | M73218.1….     |        0.03 | AY575857….    |             0 |                 |               0 |                   |                 0.01 | MN614142….              |
| TTGATTCTCAGCCCTTCGC    |         0.93 | M73218.1….     |        0.06 | JF443725….    |             0 | AB222183.1      |               0 |                   |                 0.01 | MN614142….              |
| ATWTTCATCCAACCAACCCCTT |         0.94 | M73218.1….     |        0.06 | FJ457024….    |             0 |                 |               0 |                   |                 0.01 | MN614142….              |
| TTGATTCTCAGCCCTTCGC    |         0.93 | M73218.1….     |        0.06 | JF443725….    |             0 | AB222183.1      |               0 |                   |                 0.01 | MN614142….              |
| CCTATWTTCATCCAACCAAC   |         0.94 | M73218.1….     |        0.04 | AY575857….    |             0 |                 |               0 |                   |                 0.01 | MN614142….              |
| TGATTCTCAGCCCTTCGC     |         0.94 | M73218.1….     |        0.06 | JF443725….    |             0 |                 |               0 |                   |                 0.01 | MN614142….              |

The match table can be visualized using `plotData()`:

``` r
plotData(matchTableOligos)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.

## Citation

To cite this package, please use: `citation("rprimer")`.
