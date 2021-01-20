
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> <!-- badges: end -->

## To do

  - Oligos=describe design process
  - Classes=set validity
  - All=test
  - Vignette
  - col/rowThreshold…
  - randomly select input seq

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
and (RT)-(q/d)PCR assays from a multiple DNA sequence alignment.

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
alignment of 200 hepatitis E virus sequences.

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
|        2 | 0.01 | 0.54 | 0.01 | 0.04 |     0 | 0.41 | C        |     0.90 | C     |    0.60 |     0.97 |
|        3 | 0.54 | 0.02 | 0.01 | 0.02 |     0 | 0.41 | A        |     0.92 | A     |    0.54 |     0.92 |
|        4 | 0.02 | 0.01 | 0.56 | 0.01 |     0 | 0.40 | G        |     0.93 | G     |    0.50 |     0.93 |
|        5 | 0.56 | 0.00 | 0.02 | 0.03 |     0 | 0.40 | A        |     0.93 | A     |    0.46 |     0.93 |
|        6 | 0.01 | 0.58 | 0.02 | 0.01 |     0 | 0.38 | C        |     0.93 | C     |    0.49 |     0.93 |

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
#                     conservedThreeEndPrimer =  FALSE,
#                     gcClampPrimer = FALSE,
#                     avoidThreeEndRunsPrimer = FALSE,
#                     probe = FALSE)
```

Results (first six rows):

| type   | fwd   | rev  | start | end | length | iupacSequence        | iupacSequenceRc      | coverage | degeneracy | gcContentMean | gcContentRange | tmMean | tmRange | sequence   | sequenceRc | gcContent  | tm         | roiStart | roiEnd |
| :----- | :---- | :--- | ----: | --: | -----: | :------------------- | :------------------- | -------: | ---------: | ------------: | -------------: | -----: | ------: | :--------- | :--------- | :--------- | :--------- | -------: | -----: |
| primer | FALSE | TRUE |    26 |  44 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |        1 |          4 |          0.53 |           0.11 |  57.64 |    5.01 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 60.12783…. |        1 |   7238 |
| primer | FALSE | TRUE |    26 |  45 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |        1 |          4 |          0.50 |           0.10 |  57.93 |    4.99 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 60.42434…. |        1 |   7238 |
| probe  | TRUE  | TRUE |    26 |  44 |     19 | ATGGAGGCCCAYCAGTTYA  | TRAACTGRTGGGCCTCCAT  |        1 |          4 |          0.53 |           0.11 |  56.58 |    5.01 | ATGGAGGC…. | TGAACTGG…. | 0.578947…. | 59.07612…. |        1 |   7238 |
| probe  | TRUE  | TRUE |    26 |  45 |     20 | ATGGAGGCCCAYCAGTTYAT | ATRAACTGRTGGGCCTCCAT |        1 |          4 |          0.50 |           0.10 |  56.92 |    4.98 | ATGGAGGC…. | ATGAACTG…. | 0.55, 0….. | 59.40578…. |        1 |   7238 |
| probe  | TRUE  | TRUE |    27 |  44 |     18 | TGGAGGCCCAYCAGTTYA   | TRAACTGRTGGGCCTCCA   |        1 |          4 |          0.56 |           0.11 |  56.14 |    5.20 | TGGAGGCC…. | TGAACTGG…. | 0.611111…. | 58.72365…. |        1 |   7238 |
| probe  | TRUE  | TRUE |    27 |  45 |     19 | TGGAGGCCCAYCAGTTYAT  | ATRAACTGRTGGGCCTCCA  |        1 |          4 |          0.53 |           0.11 |  56.51 |    5.15 | TGGAGGCC…. | ATGAACTG…. | 0.578947…. | 59.07612…. |        1 |   7238 |

The results can be visualized using `plotData()`.

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and combines them
with probes, if probes are present in the dataset. You can either use
default settings (as below), or customize the amplicon length and
maximum allowed difference between the primers, and between the primer
pair and probe.

``` r
myAssays <- assays(myOligos)
```

Results (first six rows):

| start |  end | ampliconLength | tmDifferencePrimer | tmDifferencePrimerProbe | totalDegeneracy | startFwd | endFwd | lengthFwd | iupacSequenceFwd       | coverageFwd | degeneracyFwd | gcContentMeanFwd | gcContentRangeFwd | tmMeanFwd | tmRangeFwd | sequenceFwd | gcContentFwd | tmFwd      | startRev | endRev | lengthRev | iupacSequenceRev      | coverageRev | degeneracyRev | gcContentMeanRev | gcContentRangeRev | tmMeanRev | tmRangeRev | sequenceRev | gcContentRev | tmRev      | plusPr | minusPr | startPr | endPr | lengthPr | iupacSequencePr        | iupacSequenceRcPr      | coveragePr | degeneracyPr | gcContentMeanPr | gcContentRangePr | tmMeanPr | tmRangePr | sequencePr | sequenceRcPr | gcContentPr | tmPr       | roiStart | roiEnd |
| ----: | ---: | -------------: | -----------------: | ----------------------: | --------------: | -------: | -----: | --------: | :--------------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | -------: | -----: | --------: | :-------------------- | ----------: | ------------: | ---------------: | ----------------: | --------: | ---------: | :---------- | :----------- | :--------- | :----- | :------ | ------: | ----: | -------: | :--------------------- | :--------------------- | ---------: | -----------: | --------------: | ---------------: | -------: | --------: | :--------- | :----------- | :---------- | :--------- | -------: | -----: |
|  5296 | 5361 |             66 |               3.59 |                    3.01 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA  |           1 |             2 |             0.52 |              0.05 |     58.01 |       1.96 | CGAAGGGG….  | 0.55, 0.5    | 58.98779…. | TRUE   | TRUE    |    5318 |  5339 |       22 | TCAGCCCTTCGCMMTCCCCTAT | ATAGGGGAKKGCGAAGGGCTGA |          1 |            4 |            0.59 |             0.09 |    62.81 |      3.38 | TCAGCCCT…. | ATAGGGGA….   | 0.545454….  | 61.13138…. |        1 |   7238 |
|  5296 | 5361 |             66 |               3.59 |                    1.11 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5342 |   5361 |        20 | CRAAGGGGTTGGTTGGATGA  |           1 |             2 |             0.52 |              0.05 |     58.01 |       1.96 | CGAAGGGG….  | 0.55, 0.5    | 58.98779…. | TRUE   | TRUE    |    5320 |  5339 |       20 | AGCCCTTCGCMMTCCCCTAT   | ATAGGGGAKKGCGAAGGGCT   |          1 |            4 |            0.60 |             0.10 |    60.91 |      3.71 | AGCCCTTC…. | ATAGGGGA….   | 0.55, 0…..  | 59.05920…. |        1 |   7238 |
|  5296 | 5362 |             67 |               0.82 |                    1.63 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5342 |   5362 |        21 | GCRAAGGGGTTGGTTGGATGA |           1 |             2 |             0.55 |              0.05 |     60.77 |       1.83 | GCGAAGGG….  | 0.571428….   | 61.68777…. | TRUE   | TRUE    |    5318 |  5339 |       22 | TCAGCCCTTCGCMMTCCCCTAT | ATAGGGGAKKGCGAAGGGCTGA |          1 |            4 |            0.59 |             0.09 |    62.81 |      3.38 | TCAGCCCT…. | ATAGGGGA….   | 0.545454….  | 61.13138…. |        1 |   7238 |
|  5296 | 5361 |             66 |               4.69 |                    3.56 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5343 |   5361 |        19 | CRAAGGGGTTGGTTGGATG   |           1 |             2 |             0.55 |              0.05 |     56.91 |       2.08 | CGAAGGGG….  | 0.578947….   | 57.95012…. | TRUE   | TRUE    |    5318 |  5339 |       22 | TCAGCCCTTCGCMMTCCCCTAT | ATAGGGGAKKGCGAAGGGCTGA |          1 |            4 |            0.59 |             0.09 |    62.81 |      3.38 | TCAGCCCT…. | ATAGGGGA….   | 0.545454….  | 61.13138…. |        1 |   7238 |
|  5296 | 5361 |             66 |               4.69 |                    1.66 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5343 |   5361 |        19 | CRAAGGGGTTGGTTGGATG   |           1 |             2 |             0.55 |              0.05 |     56.91 |       2.08 | CGAAGGGG….  | 0.578947….   | 57.95012…. | TRUE   | TRUE    |    5320 |  5339 |       20 | AGCCCTTCGCMMTCCCCTAT   | ATAGGGGAKKGCGAAGGGCT   |          1 |            4 |            0.60 |             0.10 |    60.91 |      3.71 | AGCCCTTC…. | ATAGGGGA….   | 0.55, 0…..  | 59.05920…. |        1 |   7238 |
|  5296 | 5362 |             67 |               1.72 |                    2.08 |               8 |     5296 |   5317 |        22 | TCTGGGGTGACMGGGTTGATTC |           1 |             2 |             0.57 |              0.05 |      61.6 |       2.28 | TCTGGGGT….  | 0.545454….   | 60.45857…. |     5343 |   5362 |        20 | GCRAAGGGGTTGGTTGGATG  |           1 |             2 |             0.58 |              0.05 |     59.87 |       1.93 | GCGAAGGG….  | 0.6, 0.55    | 60.84061…. | TRUE   | TRUE    |    5318 |  5339 |       22 | TCAGCCCTTCGCMMTCCCCTAT | ATAGGGGAKKGCGAAGGGCTGA |          1 |            4 |            0.59 |             0.09 |    62.81 |      3.38 | TCAGCCCT…. | ATAGGGGA….   | 0.545454….  | 61.13138…. |        1 |   7238 |

The assays can be visualized using `plotData()`:

``` r
plotData(myAssays)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
