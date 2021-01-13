
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> <!-- badges: end -->

## To do

  - Oligos=describe design process
  - Classes=set validity
  - All=test
  - Vignette
  - col/rowThreshold…

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

| type   | fwd   | rev   | start | end | length | iupacSequence          | iupacSequenceRc        | identity | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | sequence   | sequenceRc | gcContent  | tm         | roiStart | roiEnd |
| :----- | :---- | :---- | ----: | --: | -----: | :--------------------- | :--------------------- | -------: | ---------: | ------------: | -------------: | -----: | ------: | :--------- | :--------- | :--------- | :--------- | -------: | -----: |
| primer | TRUE  | FALSE |    51 |  70 |     20 | GCTCCTGGCATYACTACTGC   | GCAGTAGTRATGCCAGGAGC   |     0.99 |          2 |          0.58 |           0.05 |  58.86 |    2.40 | GCTCCTGG…. | GCAGTAGT…. | 0.6, 0.55  | 60.06040…. |        1 |   7208 |
| primer | TRUE  | FALSE |    51 |  72 |     22 | GCTCCTGGCATYACTACTGCYA | TRGCAGTAGTRATGCCAGGAGC |     0.98 |          4 |          0.55 |           0.09 |  60.74 |    4.72 | GCTCCTGG…. | TGGCAGTA…. | 0.590909…. | 63.08664…. |        1 |   7208 |
| probe  | FALSE | TRUE  |    51 |  72 |     22 | GCTCCTGGCATYACTACTGCYA | TRGCAGTAGTRATGCCAGGAGC |     0.98 |          4 |          0.55 |           0.09 |  59.84 |    4.73 | GCTCCTGG…. | TGGCAGTA…. | 0.590909…. | 62.19669…. |        1 |   7208 |
| primer | FALSE | TRUE  |    52 |  71 |     20 | CTCCTGGCATYACTACTGCY   | RGCAGTAGTRATGCCAGGAG   |     0.98 |          4 |          0.55 |           0.10 |  57.61 |    3.63 | CTCCTGGC…. | GGCAGTAG…. | 0.6, 0.5…. | 59.42049…. |        1 |   7208 |
| primer | TRUE  | TRUE  |    52 |  72 |     21 | CTCCTGGCATYACTACTGCYA  | TRGCAGTAGTRATGCCAGGAG  |     0.98 |          4 |          0.52 |           0.10 |  58.11 |    5.02 | CTCCTGGC…. | TGGCAGTA…. | 0.571428…. | 60.60395…. |        1 |   7208 |
| primer | TRUE  | TRUE  |    52 |  73 |     22 | CTCCTGGCATYACTACTGCYAT | ATRGCAGTAGTRATGCCAGGAG |     0.98 |          4 |          0.50 |           0.09 |  58.35 |    5.01 | CTCCTGGC…. | ATGGCAGT…. | 0.545454…. | 60.85270…. |        1 |   7208 |

The results can be visualized using `plotData()`.

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and combines them
with probes, if `probe = TRUE` was used in the oligo design step. You
can either use default settings (as below), or customize the amplicon
length and maximum allowed difference between the primers, and between
the primer pair and probe.

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

| start | end | ampliconLength | tmDifferencePrimer | tmDifferencePrimerProbe | totalDegeneracy | startFwd | endFwd | lengthFwd | iupacSequenceFwd       | identityFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | startRev | endRev | lengthRev | iupacSequenceRev   | identityRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | sequenceRev | gcContentRev | tmRev      | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | identityPr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | roiStart | roiEnd |
| ----: | --: | -------------: | -----------------: | ----------------------: | --------------: | -------: | -----: | --------: | :--------------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | -------: | -----: | --------: | :----------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | :----- | :------ | ------: | ----: | -------: | :--------------------- | :--------------------- | ---------: | -----------: | --------------: | ---------------: | -------: | --------: | :--------- | :----------- | :---------- | :--------- | -------: | -----: |
|    51 | 128 |             78 |               2.19 |                    2.38 |               8 |       51 |     70 |        20 | GCTCCTGGCATYACTACTGC   |        0.99 |             2 |             0.58 |              0.05 |     58.86 |       2.40 | GCTCCTGG….  | 0.6, 0.55    | 60.06040…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | FALSE   |      71 |    91 |       21 | YATTGAGCAGGCTGCTCTRGC  | GCYAGAGCAGCCTGCTCAATR  |       0.98 |            4 |            0.57 |             0.10 |    60.14 |      4.31 | CATTGAGC…. | GCTAGAGC….   | 0.571428….  | 59.47979…. |        1 |   7208 |
|    51 | 128 |             78 |               2.19 |                    2.17 |               6 |       51 |     70 |        20 | GCTCCTGGCATYACTACTGC   |        0.99 |             2 |             0.58 |              0.05 |     58.86 |       2.40 | GCTCCTGG….  | 0.6, 0.55    | 60.06040…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | FALSE   |      72 |    91 |       20 | ATTGAGCAGGCTGCTCTRGC   | GCYAGAGCAGCCTGCTCAAT   |       0.99 |            2 |            0.58 |             0.05 |    59.93 |      2.98 | ATTGAGCA…. | GCTAGAGC….   | 0.55, 0.6   | 58.44033…. |        1 |   7208 |
|    51 | 128 |             78 |               2.19 |                    3.20 |               8 |       51 |     70 |        20 | GCTCCTGGCATYACTACTGC   |        0.99 |             2 |             0.58 |              0.05 |     58.86 |       2.40 | GCTCCTGG….  | 0.6, 0.55    | 60.06040…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | TRUE    |      72 |    92 |       21 | ATTGAGCAGGCTGCTCTRGCW  | SGCYAGAGCAGCCTGCTCAAT  |       0.98 |            4 |            0.55 |             0.05 |    60.96 |      3.11 | ATTGAGCA…. | TGCTAGAG….   | 0.523809….  | 59.67288…. |        1 |   7208 |
|    51 | 128 |             78 |               2.19 |                    4.13 |               8 |       51 |     70 |        20 | GCTCCTGGCATYACTACTGC   |        0.99 |             2 |             0.58 |              0.05 |     58.86 |       2.40 | GCTCCTGG….  | 0.6, 0.55    | 60.06040…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | TRUE    |      72 |    93 |       22 | ATTGAGCAGGCTGCTCTRGCWG | CSGCYAGAGCAGCCTGCTCAAT |       0.98 |            4 |            0.57 |             0.05 |    61.89 |      2.72 | ATTGAGCA…. | CTGCTAGA….   | 0.545454….  | 60.53488…. |        1 |   7208 |
|    51 | 128 |             78 |               2.19 |                    2.94 |               8 |       51 |     70 |        20 | GCTCCTGGCATYACTACTGC   |        0.99 |             2 |             0.58 |              0.05 |     58.86 |       2.40 | GCTCCTGG….  | 0.6, 0.55    | 60.06040…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | TRUE    |      73 |    92 |       20 | TTGAGCAGGCTGCTCTRGCW   | SGCYAGAGCAGCCTGCTCAA   |       0.98 |            4 |            0.58 |             0.05 |    60.70 |      3.21 | TTGAGCAG…. | TGCTAGAG….   | 0.55, 0…..  | 59.38013…. |        1 |   7208 |
|    51 | 128 |             78 |               4.07 |                    2.00 |              10 |       51 |     72 |        22 | GCTCCTGGCATYACTACTGCYA |        0.98 |             4 |             0.55 |              0.09 |     60.74 |       4.72 | GCTCCTGG….  | 0.590909….   | 63.08664…. |      111 |    128 |        18 | AACYACCACAGCATTCGC |        0.98 |             2 |             0.53 |              0.06 |     56.67 |       3.23 | AACTACCA….  | 0.5, 0.5….   | 55.04825…. | TRUE   | TRUE    |      73 |    92 |       20 | TTGAGCAGGCTGCTCTRGCW   | SGCYAGAGCAGCCTGCTCAA   |       0.98 |            4 |            0.58 |             0.05 |    60.70 |      3.21 | TTGAGCAG…. | TGCTAGAG….   | 0.55, 0…..  | 59.38013…. |        1 |   7208 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
