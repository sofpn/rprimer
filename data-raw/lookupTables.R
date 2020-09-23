allBases <- "ACGTRYSWKMBDHVN-"
dnaBases <- "ACGT-"

complementLookup <- c(
    "A" = "T",
    "T" = "A",
    "G" = "C",
    "C" = "G",
    "R" = "Y",
    "Y" = "R",
    "M" = "K",
    "K" = "M",
    "S" = "W",
    "W" = "S",
    "N" = "N",
    "H" = "D",
    "D" = "H",
    "V" = "B",
    "B" = "V",
    "-" = "-"
  )

iupacLookup <- c(
    "-,A" = "A",
    "A" = "A",
    "-,C" = "C",
    "C" = "C",
    "-,G" = "G",
    "G" = "G",
    "-,T" = "T",
    "T" = "T",
    "-,A,G" = "R",
    "A,G" = "R",
    "-,C,T" = "Y",
    "C,T" = "Y",
    "-,C,G" = "S",
    "C,G" = "S",
    "-,A,T" = "W",
    "A,T" = "W",
    "-,G,T" = "K",
    "G,T" = "K",
    "-,A,C" = "M",
    "A,C" = "M",
    "-,C,G,T" = "B",
    "C,G,T" = "B",
    "-,A,G,T" = "D",
    "A,G,T" = "D",
    "-,A,C,T" = "H",
    "A,C,T" = "H",
    "-,A,C,G" = "V",
    "A,C,G" = "V",
    "-,A,C,G,T" = "N",
    "A,C,G,T" = "N",
    "-" = "-"
  )

degenerateLookup <- c(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "R" = "A,G",
    "Y" = "C,T",
    "S" = "C,G",
    "W" = "A,T",
    "K" = "G,T",
    "M" = "A,C",
    "B" = "C,G,T",
    "D" = "A,G,T",
    "H" = "A,C,T",
    "V" = "A,C,G",
    "N" = "A,C,G,T",
    "." = ".",
    "+" = "+",
    "-" = "-"
  )

degeneracyLookup <- c(
    "A" = 1,
    "C" = 1,
    "C" = 1,
    "G" = 1,
    "T" = 1,
    "R" = 2,
    "Y" = 2,
    "S" = 2,
    "W" = 2,
    "K" = 2,
    "M" = 2,
    "B" = 3,
    "D" = 3,
    "H" = 3,
    "V" = 3,
    "N" = 4,
    "." = 0,
    "+" = 0,
    "-" = 0
  )

bases <- c(
  "AA",
  "TT",
  "AT",
  "TA",
  "CA",
  "TG",
  "GT",
  "AC",
  "CT",
  "AG",
  "GA",
  "TC",
  "CG",
  "GC",
  "GG",
  "CC",
  "."
)
dH <- c(
  7.9,
  7.9,
  7.2,
  7.2,
  8.5,
  8.5,
  8.4,
  8.4,
  7.8,
  7.8,
  8.2,
  8.2,
  10.6,
  9.8,
  8.0,
  8.0,
  0
) * -1 * 1000
dS <- c(
  22.2,
  22.2,
  20.4,
  21.3,
  22.7,
  22.7,
  22.4,
  22.4,
  21.0,
  21.0,
  22.2,
  22.2,
  27.2,
  24.4,
  19.9,
  19.9,
  0
) * -1

nnLookup <- tibble::tibble(bases, dH, dS)

gasConstant <- 1.987

usethis::use_data(
  allBases,
  dnaBases,
  complementLookup,
  iupacLookup,
  degenerateLookup,
  degeneracyLookup,
  nnLookup,
  gasConstant,
  internal = TRUE,
  overwrite = TRUE
)
