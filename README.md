
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rprimer <img src='man/figures/rprimer.png' align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/sofpn/rprimer/workflows/R-CMD-check/badge.svg)](https://github.com/sofpn/rprimer/actions)
<!-- badges: end -->

rprimer is an R package that designs degenerate oligos and PCR assays
from a multiple DNA sequence alignment of target sequences of interest.
The package is specifically designed for sequence variable viruses.

## Installation

To install rprimer from
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/rprimer.html),
start R (version 4.2 or later) and enter:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rprimer")
```

Attach the package by calling:

``` r
library(rprimer)
```

## Overview

The package contains five main functions:

- `consensusProfile()`
- `designOligos()`
- `designAssays()`
- `checkMatch()`
- `plotData()`

## Shiny application

The package can be run through a Shiny application (a graphical user
interface). To start the application, type `runRprimerApp()` from within
R upon installing and attaching the package.

The application can also be found online,
[here](https://sofpn.shinyapps.io/rprimer).

## Workflow

### Import alignment

The first step is to import an alignment with target sequences of
interest. This is done by using `readDNAMultipleAlignment()`.

The file “example_alignment.txt” contains an alignment of 50 hepatitis E
virus sequences.

``` r
infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")

myAlignment <- readDNAMultipleAlignment(infile, format = "fasta")
```

### Step 1: `consensusProfile`

`consensusProfile()` takes a `DNAMultipleAlignment` as input and returns
all the information needed for the subsequent design process.

``` r
myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)
```

Results (row 100-106):

| position |    a |    c |    g |    t | other | gaps | majority | identity | iupac | coverage |
|---------:|-----:|-----:|-----:|-----:|------:|-----:|:---------|---------:|:------|---------:|
|      100 | 0.00 | 1.00 | 0.00 | 0.00 |  0.00 |    0 | C        |     1.00 | C     |     1.00 |
|      101 | 1.00 | 0.00 | 0.00 | 0.00 |  0.00 |    0 | A        |     1.00 | A     |     1.00 |
|      102 | 0.16 | 0.00 | 0.84 | 0.00 |  0.00 |    0 | G        |     0.84 | R     |     1.00 |
|      103 | 0.00 | 0.00 | 1.00 | 0.00 |  0.00 |    0 | G        |     1.00 | G     |     1.00 |
|      104 | 0.00 | 0.98 | 0.00 | 0.00 |  0.02 |    0 | C        |     1.00 | C     |     1.00 |
|      105 | 0.20 | 0.00 | 0.02 | 0.78 |  0.00 |    0 | T        |     0.78 | W     |     0.98 |
|      106 | 0.00 | 0.00 | 1.00 | 0.00 |  0.00 |    0 | G        |     1.00 | G     |     1.00 |

The results can be visualized with `plotData()`:

``` r
plotData(myConsensusProfile)
#> Warning in ggplot2::geom_segment(color = "#93A8AC", ggplot2::aes(x = min(position), : All aesthetics have length 1, but the data has 7597 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Step 2: `designOligos`

The next step is to design oligos. You can either use the default
settings as below, or adjust them as preferred (see the package vignette
or `?designOligos` for more information). The default settings allow a
maximum degeneracy of four, which means that only the most conserved
regions of the genome will be considered as oligo binding sites.

``` r
myOligos <- designOligos(myConsensusProfile)
```

Results (first six rows):

| type   | fwd   | rev   | start | end | length | iupacSequence        | iupacSequenceRc      | identity | coverage | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | deltaGMean | deltaGRange | sequence   | sequenceRc | gcContent  | tm         | deltaG     | method    | score | roiStart | roiEnd |
|:-------|:------|:------|------:|----:|-------:|:---------------------|:---------------------|---------:|---------:|-----------:|--------------:|---------------:|-------:|--------:|-----------:|------------:|:-----------|:-----------|:-----------|:-----------|:-----------|:----------|------:|---------:|-------:|
| probe  | TRUE  | TRUE  |   124 | 143 |     20 | TCYGCCYTGGCGAATGCTGT | ACAGCATTCGCCARGGCRGA |     0.95 |     0.99 |          4 |          0.60 |           0.10 |  63.17 |    4.33 |     -21.61 |        2.08 | TCCGCCCT…. | ACAGCATT…. | 0.65, 0….. | 65.33623…. | -22.6538…. | ambiguous |     2 |        1 |   7597 |
| probe  | FALSE | TRUE  |   127 | 146 |     20 | GCCYTGGCGAATGCTGTGGT | ACCACAGCATTCGCCARGGC |     0.98 |     0.99 |          2 |          0.62 |           0.05 |  63.18 |    1.84 |     -21.73 |        0.83 | GCCCTGGC…. | ACCACAGC…. | 0.65, 0.6  | 64.10586…. | -22.1475…. | ambiguous |     3 |        1 |   7597 |
| primer | TRUE  | FALSE |   128 | 146 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.97 |     0.99 |          2 |          0.61 |           0.05 |  61.48 |    1.95 |     -19.84 |        0.83 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 62.45335…. | -20.2570…. | ambiguous |     3 |        1 |   7597 |
| primer | TRUE  | FALSE |   128 | 147 |     20 | CCYTGGCGAATGCTGTGGTR | YACCACAGCATTCGCCARGG |     0.96 |     0.99 |          4 |          0.60 |           0.10 |  61.61 |    3.37 |     -20.55 |        1.76 | CCCTGGCG…. | TACCACAG…. | 0.6, 0.5…. | 61.77166…. | -20.5089…. | ambiguous |     2 |        1 |   7597 |
| probe  | TRUE  | TRUE  |   128 | 146 |     19 | CCYTGGCGAATGCTGTGGT  | ACCACAGCATTCGCCARGG  |     0.97 |     0.99 |          2 |          0.61 |           0.05 |  60.45 |    1.94 |     -19.84 |        0.83 | CCCTGGCG…. | ACCACAGC…. | 0.631578…. | 61.41686…. | -20.2570…. | ambiguous |     3 |        1 |   7597 |
| probe  | TRUE  | TRUE  |   128 | 147 |     20 | CCYTGGCGAATGCTGTGGTR | YACCACAGCATTCGCCARGG |     0.96 |     0.99 |          4 |          0.60 |           0.10 |  60.63 |    3.37 |     -20.55 |        1.76 | CCCTGGCG…. | TACCACAG…. | 0.6, 0.5…. | 60.78676…. | -20.5089…. | ambiguous |     2 |        1 |   7597 |

The results can be visualized as a dashboard, using `plotData()`:

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

### Step 3: `designAssays`

`designAssays()` finds pairs of forward and reverse primers and combine
them with probes, if probes are present in the input dataset. You can
either use the default settings as below, or adjust the design
constraints (see the package vignette or `?designAssays` for more
information).

``` r
myAssays <- designAssays(myOligos)
```

Results (first six rows):

| start |  end | length | totalDegeneracy | score | startFwd | endFwd | lengthFwd | iupacSequenceFwd     | identityFwd | coverageFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | deltaGMeanFwd | deltaGRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | deltaGFwd  | methodFwd | startRev | endRev | lengthRev | iupacSequenceRev     | identityRev | coverageRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | deltaGMeanRev | deltaGRangeRev | sequenceRev | gcContentRev | tmRev      | deltaGRev  | methodRev | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | identityPr | coveragePr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | deltaGMeanPr | deltaGRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | deltaGPr   | methodPr  | roiStart | roiEnd |
|------:|-----:|-------:|----------------:|------:|---------:|-------:|----------:|:---------------------|------------:|------------:|--------------:|-----------------:|------------------:|----------:|-----------:|--------------:|---------------:|:------------|:-------------|:-----------|:-----------|:----------|---------:|-------:|----------:|:---------------------|------------:|------------:|--------------:|-----------------:|------------------:|----------:|-----------:|--------------:|---------------:|:------------|:-------------|:-----------|:-----------|:----------|:-------|:--------|--------:|------:|---------:|:-----------------------|:-----------------------|-----------:|-----------:|-------------:|----------------:|-----------------:|---------:|----------:|-------------:|--------------:|:-----------|:-------------|:------------|:-----------|:-----------|:----------|---------:|-------:|
|  5605 | 5673 |     69 |               6 |  2.00 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | FALSE   |    5625 |  5642 |       18 | CMGGGTTGATTCTCAGCC     | GGCTGAGAATCAACCCKG     |       0.97 |       0.99 |            2 |            0.58 |             0.06 |    55.14 |      2.79 |       -17.06 |          1.25 | CAGGGTTG…. | GGCTGAGA….   | 0.555555….  | 53.74554…. | -16.4324…. | ambiguous |        1 |   7597 |
|  5605 | 5673 |     69 |               6 |  2.33 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | FALSE   |    5625 |  5643 |       19 | CMGGGTTGATTCTCAGCCC    | GGGCTGAGAATCAACCCKG    |       0.97 |       0.99 |            2 |            0.61 |             0.05 |    57.63 |      2.64 |       -18.54 |          1.25 | CAGGGTTG…. | GGGCTGAG….   | 0.578947….  | 56.30713…. | -17.9185…. | ambiguous |        1 |   7597 |
|  5605 | 5673 |     69 |               6 |  2.00 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | TRUE    |    5625 |  5644 |       20 | CMGGGTTGATTCTCAGCCCT   | AGGGCTGAGAATCAACCCKG   |       0.97 |       1.00 |            2 |            0.58 |             0.05 |    58.87 |      2.55 |       -19.43 |          1.25 | CAGGGTTG…. | AGGGCTGA….   | 0.55, 0.6   | 57.59836…. | -18.8035…. | ambiguous |        1 |   7597 |
|  5605 | 5673 |     69 |               6 |  1.67 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | TRUE    |    5625 |  5645 |       21 | CMGGGTTGATTCTCAGCCCTT  | AAGGGCTGAGAATCAACCCKG  |       0.98 |       1.00 |            2 |            0.55 |             0.05 |    59.21 |      2.43 |       -20.08 |          1.25 | CAGGGTTG…. | AAGGGCTG….   | 0.523809….  | 57.99472…. | -19.4553…. | ambiguous |        1 |   7597 |
|  5605 | 5673 |     69 |               6 |  2.00 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | FALSE   |    5625 |  5646 |       22 | CMGGGTTGATTCTCAGCCCTTC | GAAGGGCTGAGAATCAACCCKG |       0.98 |       0.99 |            2 |            0.57 |             0.05 |    59.91 |      2.28 |       -21.11 |          1.25 | CAGGGTTG…. | GAAGGGCT….   | 0.545454….  | 58.77533…. | -20.4881…. | ambiguous |        1 |   7597 |
|  5605 | 5673 |     69 |               6 |  2.00 |     5605 |   5624 |        20 | GGCRGTGGTTTCTGGGGTGA |        0.98 |           1 |             2 |             0.62 |              0.05 |     62.84 |       2.51 |        -20.86 |           1.25 | GGCAGTGG….  | 0.6, 0.65    | 61.57995…. | -20.2350…. | ambiguous |     5654 |   5673 |        20 | GTTGGTTGGATGAASATAGG |           1 |           1 |             2 |              0.4 |                 0 |     50.71 |        1.1 |        -15.27 |           0.52 | GTTGGTTG….  | 0.4, 0.4     | 50.15469…. | -15.0078…. | ambiguous | TRUE   | FALSE   |    5626 |  5643 |       18 | MGGGTTGATTCTCAGCCC     | GGGCTGAGAATCAACCCK     |       0.97 |       0.99 |            2 |            0.58 |             0.06 |    55.71 |      1.65 |       -17.21 |          0.94 | AGGGTTGA…. | GGGCTGAG….   | 0.555555….  | 54.88806…. | -16.7409…. | ambiguous |        1 |   7597 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Additional step: `checkMatch`

`checkMatch()` shows the proportion and names of the target sequences in
the input alignment that match with the generated oligos or assays. See
the package vignette or `?checkMatch` for more information.

``` r
## Randomly select six oligos to illustrate an example 
selection <- sample(seq_len(nrow(myOligos)), size = 6)

matchTableOligos <- checkMatch(myOligos[selection, ], target = myAlignment)
```

Results:

| iupacSequence          | perfectMatch | idPerfectMatch | oneMismatch | idOneMismatch | twoMismatches | idTwoMismatches | threeMismatches | idThreeMismatches | fourOrMoreMismatches | idFourOrMoreMismatches | offTargetMatch | idOffTargetMatch |
|:-----------------------|-------------:|:---------------|------------:|:--------------|--------------:|:----------------|----------------:|:------------------|---------------------:|:-----------------------|---------------:|:-----------------|
| ACMGGGTTGATTCTCAGCCCTT |         0.90 | AB073912….     |        0.10 | AB481228….    |          0.00 |                 |               0 |                   |                    0 |                        |              0 |                  |
| TATWTTCATCCAACCAACCCC  |         0.98 | AB073912….     |        0.02 | HM439284.1    |          0.00 |                 |               0 |                   |                    0 |                        |              0 |                  |
| CRGTGGTTTCTGGGGTGAC    |         0.96 | AB073912….     |        0.04 | BD378055….    |          0.00 |                 |               0 |                   |                    0 |                        |              0 |                  |
| TTCATCCAACCAACCCCT     |         0.98 | AB073912….     |        0.02 | HM439284.1    |          0.00 |                 |               0 |                   |                    0 |                        |              0 |                  |
| GGTGACMGGGTTGATTCT     |         0.90 | AB073912….     |        0.08 | BD378055….    |          0.02 | JQ953665.1      |               0 |                   |                    0 |                        |              0 |                  |
| TGGTTTCTGGGGTGACMG     |         0.96 | AB073912….     |        0.04 | BD378055….    |          0.00 |                 |               0 |                   |                    0 |                        |              0 |                  |

The match table can be visualized using `plotData()`:

``` r
plotData(matchTableOligos)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## More information

Please see the [package
vignette](https://bioconductor.org/packages/devel/bioc/vignettes/rprimer/inst/doc/getting-started-with-rprimer.html)
for more information on how to use the package.

## Citation

Persson S., Larsson C., Simonsson M., Ellström P. (2022) rprimer: an
R/bioconductor package for design of degenerate oligos for sequence
variable viruses. [*BMC Bioinformatics*
23:239](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04781-0)
