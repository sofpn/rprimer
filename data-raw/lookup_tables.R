complement_lookup <-
  c(
    "a" = "t",
    "t" = "a",
    "g" = "c",
    "c" = "g",
    "r" = "y",
    "y" = "r",
    "m" = "k",
    "k" = "m",
    "s" = "w",
    "w" = "s",
    "n" = "n",
    "h" = "d",
    "d" = "h",
    "v" = "b",
    "b" = "v",
    "-" = "-",
    "." = "."
  )

iupac_lookup <-
  c(
    "-,a" = "a",
    "a" = "a",
    "-,c" = "c",
    "c" = "c",
    "-,g" = "g",
    "g" = "g",
    "-,t" = "t",
    "t" = "t",
    "-,a,g" = "r",
    "a,g" = "r",
    "-,c,t" = "y",
    "c,t" = "y",
    "-,c,g" = "s",
    "c,g" = "s",
    "-,a,t" = "w",
    "a,t" = "w",
    "-,g,t" = "k",
    "g,t" = "k",
    "-,a,c" = "m",
    "a,c" = "m",
    "-,c,g,t" = "b",
    "c,g,t" = "b",
    "-,a,g,t" = "d",
    "a,g,t" = "d",
    "-,a,c,t" = "h",
    "a,c,t" = "h",
    "-,a,c,g" = "v",
    "a,c,g" = "v",
    "-,a,c,g,t" = "n",
    "a,c,g,t" = "n",
    "." = ".",
    "-" = "-"
  )

degenerate_lookup <-
  c(
    "a" = "a",
    "c" = "c",
    "g" = "g",
    "t" = "t",
    "r" = "a,g",
    "y" = "c,t",
    "s" = "c,g",
    "w" = "a,t",
    "k" = "g,t",
    "m" = "a,c",
    "b" = "c,g,t",
    "d" = "a,g,t",
    "h" = "a,c,t",
    "v" = "a,c,g",
    "n" = "a,c,g,t",
    "." = ".",
    "-" = "-"
  )

degeneracy_lookup <-
  c(
    "a" = 1,
    "c" = 1,
    "c" = 1,
    "g" = 1,
    "t" = 1,
    "r" = 2,
    "y" = 2,
    "s" = 2,
    "w" = 2,
    "k" = 2,
    "m" = 2,
    "b" = 3,
    "d" = 3,
    "h" = 3,
    "v" = 3,
    "n" = 4,
    "." = 0,
    "-" = 0
  )

bases <- c(
  "aa",
  "tt",
  "at",
  "ta",
  "ca",
  "tg",
  "gt",
  "ac",
  "ct",
  "ag",
  "ga",
  "tc",
  "cg",
  "gc",
  "gg",
  "cc",
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

nearest_neighbor_lookup <- tibble::tibble(bases, dH, dS)

gas_constant <- 1.987

usethis::use_data(
  complement_lookup,
  iupac_lookup,
  degenerate_lookup,
  degeneracy_lookup,
  nearest_neighbor_lookup,
  gas_constant,
  internal = TRUE,
  overwrite = TRUE
)
