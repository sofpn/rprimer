## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- internal setup, echo=FALSE, warning=FALSE, message=FALSE----------------
library(kableExtra)
# devtools::load_all(".")
library(rprimer)

## -----------------------------------------------------------------------------
filepath <- system.file("extdata", "example_alignment.txt", package = "rprimer")

myAlignment <- Biostrings::readDNAMultipleAlignment(filepath, format = "fasta")

## -----------------------------------------------------------------------------
## Mask everything but position 3000 to 4000 and 5000 to 6000
myMaskedAlignment <- myAlignment

Biostrings::colmask(myMaskedAlignment, invert = TRUE) <- c(3000:4000, 5000:6000)

## -----------------------------------------------------------------------------
myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)

## ---- echo=FALSE--------------------------------------------------------------
myConsensusProfile[100:110, ] %>%
  kbl(., digits = 2) %>%
  kable_material(c("striped", "hover"), position = "right") %>%
  scroll_box(width = "650px", height = "300px")

## ---- fig.width=12, fig.height=6----------------------------------------------
plotData(myConsensusProfile)

## ---- fig.width=12, fig.height=6----------------------------------------------
## Select position 5000 to 5800 in the consensus profile 
selection <- myConsensusProfile[
  myConsensusProfile$position >= 5000 & myConsensusProfile$position <= 5800, 
  ]

plotData(selection)

## ---- fig.width=12, fig.height=6----------------------------------------------
## Select position 1000 to 1020 in the consensus profile 
ntSelection <- myConsensusProfile[
  myConsensusProfile$position >= 1000 & myConsensusProfile$position <= 1020, 
]

plotData(ntSelection, type = "nucleotide")

## ---- fig.width=12, fig.height=6, warn=FALSE----------------------------------
myMaskedConsensusProfile <- consensusProfile(myMaskedAlignment, ambiguityThreshold = 0.05)

plotData(myMaskedConsensusProfile)

## -----------------------------------------------------------------------------
myOligos <- oligos(myConsensusProfile)

## ---- echo=FALSE--------------------------------------------------------------
myOligos[1:10, ] %>%
  kbl(., digits = 2) %>%
  kable_material(c("striped", "hover"), position = "right") %>%
  scroll_box(width = "650px", height = "300px")

## -----------------------------------------------------------------------------
myOligos$sequence[[1]] ## All sequence variants of the first oligo (i.e., first row) 
myOligos$gcContent[[1]] ## GC-content of all variants of the first oligo 
myOligos$tm[[1]] ## Tm of all variants of the first oligo 

## ---- fig.align="center", fig.width=12, fig.height=8--------------------------
plotData(myOligos)

## -----------------------------------------------------------------------------
myOligosModified <- oligos(myConsensusProfile, 
                           maxDegeneracyPrimer = 32, 
                           tmPrimer = c(50, 70),
                           lengthPrimer = c(16, 24), 
                           maxDegeneracyProbe = 1)

## ---- fig.align="center", fig.width=12, fig.height=8--------------------------
plotData(myOligosModified)

