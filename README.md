
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> <!-- badges: end -->

## TODO

  - Oligos=describe design process - new defaults? new defaults assays?
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
and (RT)-(q/dd)PCR assays from a multiple DNA sequence alignment. It is
especially developed for sequence variable targets, such as RNA-viruses.

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
with probes, if `probe = TRUE` was used in the oligo design step. You
can design assays with default settings, or by specifying the amplicon
length and maximum allowed difference between the primers, and the
primer pair and probe.

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

| start | end | ampliconLength | tmDifferencePrimer | tmDifferencePrimerProbe | totalDegeneracy | startFwd | endFwd | lengthFwd | iupacSequenceFwd       | identityFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | startRev | endRev | lengthRev | iupacSequenceRev     | identityRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | sequenceRev | gcContentRev | tmRev      | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | identityPr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | roiStart | roiEnd |
| ----: | --: | -------------: | -----------------: | ----------------------: | --------------: | -------: | -----: | --------: | :--------------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | -------: | -----: | --------: | :------------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | :----- | :------ | ------: | ----: | -------: | :--------------------- | :--------------------- | ---------: | -----------: | --------------: | ---------------: | -------: | --------: | :--------- | :----------- | :---------- | :--------- | -------: | -----: |
|    39 | 130 |             92 |               0.31 |                    0.19 |              10 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      111 |    130 |        20 | CKAACYACCACAGCATTCGC |        0.98 |             4 |             0.55 |               0.1 |     58.84 |       5.91 | CTAACTAC….  | 0.5, 0.5….   | 55.87255…. | TRUE   | TRUE    |      63 |    84 |       22 | ACTACTGCYATTGAGCAGGCTG | CAGCCTGCTCAATRGCAGTAGT |       0.99 |            2 |            0.52 |             0.05 |    59.19 |      2.69 | ACTACTGC…. | CAGCCTGC….   | 0.545454….  | 60.53103…. |        1 |   7208 |
|    39 | 130 |             92 |               0.31 |                    1.80 |              10 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      111 |    130 |        20 | CKAACYACCACAGCATTCGC |        0.98 |             4 |             0.55 |               0.1 |     58.84 |       5.91 | CTAACTAC….  | 0.5, 0.5….   | 55.87255…. | TRUE   | TRUE    |      65 |    86 |       22 | TACTGCYATTGAGCAGGCTGCT | AGCAGCCTGCTCAATRGCAGTA |       0.99 |            2 |            0.52 |             0.05 |    60.79 |      2.72 | TACTGCCA…. | AGCAGCCT….   | 0.545454….  | 62.15039…. |        1 |   7208 |
|    39 | 130 |             92 |               0.31 |                    2.24 |              12 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      111 |    130 |        20 | CKAACYACCACAGCATTCGC |        0.98 |             4 |             0.55 |               0.1 |     58.84 |       5.91 | CTAACTAC….  | 0.5, 0.5….   | 55.87255…. | TRUE   | TRUE    |      68 |    89 |       22 | TGCYATTGAGCAGGCTGCTCTR | YAGAGCAGCCTGCTCAATRGCA |       0.98 |            4 |            0.55 |             0.09 |    61.23 |      4.09 | TGCCATTG…. | TAGAGCAG….   | 0.545454….  | 61.89249…. |        1 |   7208 |
|    39 | 131 |             93 |               0.96 |                    0.51 |              10 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      112 |    131 |        20 | CCKAACYACCACAGCATTCG |        0.98 |             4 |             0.55 |               0.1 |     58.20 |       5.97 | CCTAACTA….  | 0.5, 0.5….   | 55.19725…. | TRUE   | TRUE    |      63 |    84 |       22 | ACTACTGCYATTGAGCAGGCTG | CAGCCTGCTCAATRGCAGTAGT |       0.99 |            2 |            0.52 |             0.05 |    59.19 |      2.69 | ACTACTGC…. | CAGCCTGC….   | 0.545454….  | 60.53103…. |        1 |   7208 |
|    39 | 131 |             93 |               0.96 |                    2.12 |              10 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      112 |    131 |        20 | CCKAACYACCACAGCATTCG |        0.98 |             4 |             0.55 |               0.1 |     58.20 |       5.97 | CCTAACTA….  | 0.5, 0.5….   | 55.19725…. | TRUE   | TRUE    |      65 |    86 |       22 | TACTGCYATTGAGCAGGCTGCT | AGCAGCCTGCTCAATRGCAGTA |       0.99 |            2 |            0.52 |             0.05 |    60.79 |      2.72 | TACTGCCA…. | AGCAGCCT….   | 0.545454….  | 62.15039…. |        1 |   7208 |
|    39 | 131 |             93 |               0.96 |                    2.56 |              12 |       39 |     60 |        22 | CAGTTYATYAAGGCTCCTGGCA |        0.98 |             4 |              0.5 |              0.09 |     59.15 |       4.42 | CAGTTCAT….  | 0.545454….   | 61.35995…. |      112 |    131 |        20 | CCKAACYACCACAGCATTCG |        0.98 |             4 |             0.55 |               0.1 |     58.20 |       5.97 | CCTAACTA….  | 0.5, 0.5….   | 55.19725…. | TRUE   | TRUE    |      68 |    89 |       22 | TGCYATTGAGCAGGCTGCTCTR | YAGAGCAGCCTGCTCAATRGCA |       0.98 |            4 |            0.55 |             0.09 |    61.23 |      4.09 | TGCCATTG…. | TAGAGCAG….   | 0.545454….  | 61.89249…. |        1 |   7208 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
