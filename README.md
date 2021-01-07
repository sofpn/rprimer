
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start --> <!-- badges: end -->

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

rprimer provides tools for designing broadly reactive (RT)-(q/dd)PCR
assays from a multiple DNA sequence alignment. The design process is
built on three functions:

  - `consensusProfile()`: takes a multiple DNA sequence alignment and
    returns all the information needed for the subsequent design
    process, arranged in an `RprimerProfile` object
  - `oligos()`: takes an `RprimerProfile` as input and designs primers
    and probes (if selected) based on user specified constraints on
    e.g. length, GC-content, melting temperature, maximum degeneracy
    and minimum end-conservation. Returns an `RprimerOligo` object
  - `assays()`: takes an `RprimerOligo` as input and designs assays with
    desired length and difference in melting temperature. Returns an
    `RprimerAssay` object

The `Rprimer`-classes are extensions of the `DataFrame` class from
S4Vectors, and behave in a similar manner as traditional data frames.

## Workflow

### Import data

The first step is to import an alignment with target sequences of
interest and, if preferred, mask positions with e.g. high gap frequency.
`readDNAMultipleAlignment()` and `maskGaps()` from Biostrings do the
work for this part.

``` r
infile <- system.file('extdata', 'example_alignment.txt', package = 'rprimer')

myAlignment <- Biostrings::readDNAMultipleAlignment(infile, format = "fasta")
myAlignment <- Biostrings::maskGaps(myAlignment, 
                                    min.fraction = 0.5, 
                                    min.block.width = 1) 
```

### Step 1: `consensusProfile`

`consensusProfile()` takes a `Biostrings::DNAMultipleAlignment` object
as input and returns all the information needed for the subsequent
design process.

``` r
myConsensusProfile <- consensusProfile(myAlignment, iupacThreshold = 0.05)
```

Output:

| position |    a |    c |    g |    t | other | gaps | majority | identity | iupac | entropy |
| -------: | ---: | ---: | ---: | ---: | ----: | ---: | :------- | -------: | :---- | ------: |
|        1 | 0.00 | 0.00 | 0.59 | 0.00 |     0 | 0.41 | G        |     1.00 | G     |    0.00 |
|        2 | 0.00 | 0.00 | 0.71 | 0.00 |     0 | 0.29 | G        |     1.00 | G     |    0.00 |
|        3 | 0.00 | 0.71 | 0.00 | 0.00 |     0 | 0.29 | C        |     1.00 | C     |    0.00 |
|        4 | 0.71 | 0.00 | 0.00 | 0.00 |     0 | 0.29 | A        |     1.00 | A     |    0.00 |
|        5 | 0.00 | 0.00 | 0.70 | 0.01 |     0 | 0.29 | G        |     0.99 | G     |    0.11 |
|        6 | 0.71 | 0.00 | 0.00 | 0.00 |     0 | 0.29 | A        |     1.00 | A     |    0.00 |

The output can be visualized with `plotData()`:

``` r
plotData(myConsensusProfile)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Step 2: `oligos`

`oligos()` searches for oligos from an `RprimerProfile`-object.

``` r
myOligos <- oligos(myConsensusProfile)
```

The oligo candidates can be visualized using `plotData()`:

``` r
plotData(myOligos)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Step 3: `assays`

`assays()` finds pairs of forward and reverse primers and, if selected,
combines them with probes.

``` r
#myAssays <- assays(myOligos)
```

The assay candidates can be visualized using `plotData()`:

``` r
#plotData(myAssays)
```

## More information

The package vignette contains more information. It is loaded by
`browseVignettes("rprimer")`.
