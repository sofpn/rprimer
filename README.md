
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

![R build
status](https://github.com/sofpn/rprimer/workflows/R-CMD-check/badge.svg)
[![Codecov test
coverage](https://codecov.io/gh/sofpn/rprimer/branch/master/graph/badge.svg)](https://codecov.io/gh/sofpn/rprimer?branch=master)
<!-- badges: end -->

*This package is in development, and will only work for R \>4.1.*

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

## Citation

To cite rprimer, please use: `citation("rprimer")`.

## Overview

rprimer provides tools for visualization of sequence conservation and
generating degenerate DNA oligos from a multiple sequence alignment. It
is especially developed for sequence variable viruses.

The package contains five functions:

  - `consensusProfile()`
  - `oligos()`
  - `assays()`
  - `plotData()`
  - `convertToDNAStringSet()`

## Workflow

### Import alignment

The first step is to import an alignment with target sequences of
interest and, if preferred, mask positions with high gap frequency.
Please use the `readDNAMultipleAlignment()` and `maskGaps()` from the
Biostrings for this part.

The file “example\_alignment.txt” contains an alignment of 200 hepatitis
E virus sequences.

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
myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)
```

Results (first six rows):

| position |    a |    c |    g |    t | other | gaps | majority | identity | iupac | entropy | coverage |
| -------: | ---: | ---: | ---: | ---: | ----: | ---: | :------- | -------: | :---- | ------: | -------: |
|        1 | 0.00 | 0.00 | 0.56 | 0.02 |     0 | 0.41 | G        |     0.95 | G     |    0.35 |     0.95 |
|        2 | 0.01 | 0.54 | 0.01 | 0.04 |     0 | 0.41 | C        |     0.90 | Y     |    0.60 |     0.97 |
|        3 | 0.54 | 0.02 | 0.01 | 0.02 |     0 | 0.41 | A        |     0.92 | A     |    0.54 |     0.92 |
|        4 | 0.02 | 0.01 | 0.56 | 0.01 |     0 | 0.40 | G        |     0.93 | G     |    0.50 |     0.93 |
|        5 | 0.56 | 0.00 | 0.02 | 0.03 |     0 | 0.40 | A        |     0.93 | A     |    0.46 |     0.93 |
|        6 | 0.01 | 0.58 | 0.02 | 0.01 |     0 | 0.38 | C        |     0.93 | C     |    0.49 |     0.93 |

The results can be visualized with `plotData()`. You can either plot the
entire genome:

``` r
plotData(myConsensusProfile)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Or zoom into a specific region of interest:

``` r
roi <- myConsensusProfile[myConsensusProfile$position >= 5000 & myConsensusProfile$position <= 5500, ]
plotData(roi)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

The nucleotide distribution can be shown, preferably within a short
range, by specifying `type = "nucleotide`:

``` r
roi2 <- myConsensusProfile[myConsensusProfile$position >= 150 & myConsensusProfile$position <= 170, ]
plotData(roi2, type = "nucleotide")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

### Step 2: `oligos`

The next step is to design oligos. You can either use the default
settings as below, or adjust the design constraints (see the package
vignette or `?oligos` for more information).

``` r
myOligos <- oligos(myConsensusProfile)
```

Results (first six rows):

| type   | fwd   | rev   | start | end | length | iupacSequence        | iupacSequenceRc      | identity | coverage | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | deltaGMean | deltaGRange | sequence   | sequenceRc | gcContent  | tm         | deltaG      | method    | roiStart | roiEnd |
| :----- | :---- | :---- | ----: | --: | -----: | :------------------- | :------------------- | -------: | -------: | ---------: | ------------: | -------------: | -----: | ------: | ---------: | ----------: | :--------- | :--------- | :--------- | :--------- | :---------- | :-------- | -------: | -----: |
| primer | FALSE | TRUE  |    26 |  43 |     18 | ATGGAGGCCCAYCAGTTY   | RAACTGRTGGGCCTCCAT   |     0.95 |        1 |          4 |          0.56 |           0.11 |  54.92 |    3.27 |     \-8.88 |        1.30 | ATGGAGGC…. | GAACTGGT…. | 0.611111…. | 56.54057…. | \-9.52761…. | ambiguous |        1 |   7238 |
| primer | FALSE | TRUE  |    26 |  44 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |     0.95 |        1 |          4 |          0.53 |           0.11 |  55.62 |    5.22 |     \-9.11 |        2.23 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 58.22863…. | \-10.2223…. | ambiguous |        1 |   7238 |
| primer | FALSE | TRUE  |    26 |  45 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |     0.96 |        1 |          4 |          0.50 |           0.10 |  55.91 |    4.97 |     \-9.17 |        2.23 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 58.39341…. | \-10.2845…. | ambiguous |        1 |   7238 |
| probe  | TRUE  | FALSE |    26 |  43 |     18 | ATGGAGGCCCAYCAGTTY   | RAACTGRTGGGCCTCCAT   |     0.95 |        1 |          4 |          0.56 |           0.11 |  53.82 |    3.27 |     \-8.88 |        1.30 | ATGGAGGC…. | GAACTGGT…. | 0.611111…. | 55.44747…. | \-9.52761…. | ambiguous |        1 |   7238 |
| probe  | TRUE  | TRUE  |    26 |  44 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |     0.95 |        1 |          4 |          0.53 |           0.11 |  54.57 |    5.20 |     \-9.11 |        2.23 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 57.17288…. | \-10.2223…. | ambiguous |        1 |   7238 |
| probe  | TRUE  | TRUE  |    26 |  45 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |     0.96 |        1 |          4 |          0.50 |           0.10 |  54.91 |    4.96 |     \-9.17 |        2.23 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 57.38719…. | \-10.2845…. | ambiguous |        1 |   7238 |

The results can be visualized using `plotData()`:

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and combines them
with probes, if probes are present in the dataset. You can either use
the default settings as below, or adjust the design constraints (see the
package vignette or `?assays` for more information).

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

| start |  end | ampliconLength | tmDifferencePrimer | tmDifferencePrimerProbe | deltaGDifferencePrimer | totalDegeneracy | startFwd | endFwd | lengthFwd | iupacSequenceFwd     | identityFwd | coverageFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | deltaGMeanFwd | deltaGRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | deltaGFwd   | methodFwd | startRev | endRev | lengthRev | iupacSequenceRev     | identityRev | coverageRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | deltaGMeanRev | deltaGRangeRev | sequenceRev | gcContentRev | tmRev      | deltaGRev   | methodRev | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | identityPr | coveragePr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | deltaGMeanPr | deltaGRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | deltaGPr    | methodPr  | deltaGDifferencePrimerProbe | roiStart | roiEnd |
| ----: | ---: | -------------: | -----------------: | ----------------------: | ---------------------: | --------------: | -------: | -----: | --------: | :------------------- | ----------: | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | ------------: | -------------: | :---------- | :----------- | :--------- | :---------- | :-------- | -------: | -----: | --------: | :------------------- | ----------: | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | ------------: | -------------: | :---------- | :----------- | :--------- | :---------- | :-------- | :----- | :------ | ------: | ----: | -------: | :--------------------- | :--------------------- | ---------: | ---------: | -----------: | --------------: | ---------------: | -------: | --------: | -----------: | ------------: | :--------- | :----------- | :---------- | :--------- | :---------- | :-------- | --------------------------: | -------: | -----: |
|  5286 | 5361 |             76 |               4.63 |                    0.45 |                   2.13 |               6 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | TRUE    |    5307 |  5328 |       22 | MGGGTTGATTCTCAGCCCTTCG | CGAAGGGCTGAGAATCAACCCK |       0.98 |          1 |            2 |            0.57 |             0.05 |    58.98 |      1.25 |      \-11.00 |          0.64 | AGGGTTGA…. | CGAAGGGC….   | 0.545454….  | 58.34934…. | \-10.6782…. | ambiguous |                      \-0.66 |        1 |   7238 |
|  5286 | 5361 |             76 |               4.63 |                    0.20 |                   2.13 |               8 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | FALSE   |    5312 |  5333 |       22 | TGATTCTCAGCCCTTCGCMMTC | GAKKGCGAAGGGCTGAGAATCA |       0.98 |          1 |            4 |            0.55 |             0.09 |    58.73 |      3.29 |      \-10.86 |          1.70 | TGATTCTC…. | GATTGCGA….   | 0.5, 0.5….  | 57.08416…. | \-10.0128…. | ambiguous |                      \-0.53 |        1 |   7238 |
|  5286 | 5361 |             76 |               4.63 |                    0.13 |                   2.13 |               8 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | FALSE   |    5314 |  5334 |       21 | ATTCTCAGCCCTTCGCMMTCC  | GGAKKGCGAAGGGCTGAGAAT  |       0.98 |          1 |            4 |            0.57 |             0.10 |    58.65 |      3.46 |      \-10.83 |          1.70 | ATTCTCAG…. | GGATTGCG….   | 0.523809….  | 56.92118…. | \-9.98307…. | ambiguous |                      \-0.50 |        1 |   7238 |
|  5286 | 5361 |             76 |               4.63 |                    2.14 |                   2.13 |               8 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | FALSE   |    5314 |  5335 |       22 | ATTCTCAGCCCTTCGCMMTCCC | GGGAKKGCGAAGGGCTGAGAAT |       0.98 |          1 |            4 |            0.59 |             0.09 |    60.67 |      3.33 |      \-11.86 |          1.70 | ATTCTCAG…. | GGGATTGC….   | 0.545454….  | 59.00487…. | \-11.0118…. | ambiguous |                      \-1.53 |        1 |   7238 |
|  5286 | 5361 |             76 |               4.63 |                    2.68 |                   2.13 |               8 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | TRUE    |    5318 |  5339 |       22 | TCAGCCCTTCGCMMTCCCCTAT | ATAGGGGAKKGCGAAGGGCTGA |       0.98 |          1 |            4 |            0.59 |             0.09 |    61.21 |      3.40 |      \-12.13 |          1.70 | TCAGCCCT…. | ATAGGGGA….   | 0.545454….  | 59.51325…. | \-11.2768…. | ambiguous |                      \-1.79 |        1 |   7238 |
|  5286 | 5361 |             76 |               4.63 |                    0.60 |                   2.13 |               8 |     5286 |   5305 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     60.84 |       2.51 |        \-11.4 |           1.17 | GGCAGTGG….  | 0.6, 0.65    | 59.58519…. | \-10.8196…. | ambiguous |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA |        0.98 |           1 |             2 |             0.52 |              0.05 |     56.21 |          2 |        \-9.27 |            0.9 | CGAAGGGG….  | 0.55, 0.5    | 57.20928…. | \-9.71979…. | ambiguous | TRUE   | TRUE    |    5320 |  5339 |       20 | AGCCCTTCGCMMTCCCCTAT   | ATAGGGGAKKGCGAAGGGCT   |       0.98 |          1 |            4 |            0.60 |             0.10 |    59.12 |      3.73 |      \-11.07 |          1.70 | AGCCCTTC…. | ATAGGGGA….   | 0.55, 0…..  | 57.25513…. | \-10.2183…. | ambiguous |                      \-0.73 |        1 |   7238 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

### Convert to fasta format

`convertToDNAStringSet()` converts oligos and assays to a
`Biostrings::DNAStringSet`-object, which can be exported to a fasta
file.

``` r
## Convert the first two assays in myAssays 
myAssaysConverted <- convertToDNAStringSet(myAssays[1:2, ], revAsRc = FALSE)

## Save as fasta-format
Biostrings::writeXStringSet(myAssaysConverted, file = "myAssays.txt")
```

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
