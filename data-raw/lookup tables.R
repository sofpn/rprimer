complement <- c(
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
    "-" = "-",
    "other" = "other"
)

iupac <- c(
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
    "-" = "-",
    "other" = "other"
)

degenerates <- c(
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

degeneracy <- c(
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
    "-" = 1
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
    "Initiation",
    "AT_penalty",
    "Symmetry_corr"
)

dH <- c(
    7.6,
    7.6,
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
    -0.2,
    -2.2,
    0
) * -1

dS <- c(
    21.3,
    21.3,
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
    5.7,
    -6.9,
    1.4
) * -1

nn <- data.frame(bases, dH, dS)

lookup <- list(
    "complement" = complement,
    "iupac" = iupac,
    "degenerates" = degenerates,
    "degeneracy" = degeneracy,
    "nn" = nn
)

usethis::use_data(lookup, internal = TRUE, overwrite = TRUE)
