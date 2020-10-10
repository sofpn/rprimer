## Biostrings::vmatchPattern(DNAStringSet(aln)) -> irange object

#' Get oligos
#'
#' \code{getOligos} identifies oligos (primers and probes) from
#' sequence properties.
#'
#' @param x an RprimerProperties object.
#'
#' @param length
#' Oligo length. The minimum allowed
#' value is 14 and the maximum allowed value is 30.
#' It defaults to \code{18:22}.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency.
#' A number [0, 1]. It defaults to 0.1, which means that
#' positions with a gap frequency above 0.1 will not be
#' considered as an oligo region.
#'
#' @param maxDegeneracy
#' Maximum number of variants.
#' The minimum allowed value is 1 and the maximum
#' allowed value is 32. It defaults to 4.
#'
#' @param gcClamp
#' \code{TRUE} or \code{FALSE}.
#' If oligos with no GC-clamp
#' should be replaced with \code{NA}
#' (recommended for primers). It defaults to \code{TRUE}.
#'
#' @param avoid3EndRuns
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}
#' (recommended for primers).
#' It defaults to \code{TRUE}.
#'
#' @param avoid5EndG
#' \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}
#' (recommended for probes). It defaults to \code{FALSE}.
#'
#' @param gcRange
#' Accepted GC-content (proportion, not %). A numeric vector [0, 1].
#' It defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRange
#' Accepted Tm. A numeric vector [30, 90].
#' It defaults to \code{c(48, 70)}.
#'
#' @param concOligo
#' Oligo concentration in M, for Tm calculation. A numeric vector
#' [0.2e-07, 2e-06], i.e. between 20 nM and 20000 nM.
#' It defaults to 5e-07 M (500 nM).
#'
#' @param concNa
#' The sodium ion concentration in M, for Tm calculation.
#' A numeric vector [0.01, 1]. It defaults to 0.05 M (50 mM).
#'
#' @param showAllVariants
#' If sequence, GC-content and Tm should be presented for all
#' variants of each oligo (in case of degenerate bases).
#' \code{TRUE} or \code{FALSE}. It defaults to \code{TRUE}.
#'
#' @section Excluded oligos:
#' The function excludes all oligos with:
#' * More than than three consecutive runs of the same dinucleotide
#' (e.g. 'TATATATA')
#' * More than four consecutive runs of the
#' same nucleotide (e.g. 'AAAAA')
#' It also excludes oligos that are duplicated, to prevent binding to several
#' places in the target.
#' These checks are made on the majority oligo sequences.
#'
#' @section Tm:
#' The melting temperature is calculated using the nearest-neighbor method,
#' with the following assumptions:
#'
#' * Oligos are not expected to be self-complementary (hence no symmetry
#' correction is done)
#' * The oligo concentration is assumed to be much higher
#' than the target concentration
#'
#' See references for table values and equations.
#'
#' @section Note:
#' GC-content and Tm are calculated based on the majority oligos, and
#' may thus be misleading for degenerate (IUPAC) oligos.
#'
#' @return
#' A tibble (a data frame) with all oligo
#' candidates. An error message will return if no oligos are found.
#'
#' The tibble contains the following information:
#'
#' \describe{
#'   \item{Begin}{position where the oligo begins}
#'   \item{End}{position where the oligo ends}
#'   \item{Length}{length of the oligo}
#'   \item{Majority}{majority sequence}
#'   \item{IUPAC}{IUPAC sequence (i.e. with degenerate bases)}
#'   \item{Majority_RC}{majority sequence, reverse complement}
#'   \item{IUPAC_RC}{IUPAC sequence, reverse complement}
#'   \item{Degeneracy}{number of variants}
#'   \item{GC_majority}{GC content (majority sequence), proportion}
#'   \item{Tm_majority}{melting temperature}
#' }
#'
#' @examples
#' data("exampleRprimerProperties")
#'
#' getOligos <- function(
#' exampleRprimerProperties,
#' length = 18:22,
#' maxGapFrequency = 0.1,
#' maxDegeneracy = 4,
#' gcClamp = TRUE,
#' avoid3EndRuns = TRUE,
#' avoid5EndG = FALSE,
#' gcRange = c(0.45, 0.55),
#' tmRange = c(48, 65),
#' concOligo = 5e-07,
#' concNa = 0.05,
#' showAllVariants = TRUE
#' )
#'
#' @references
#' Tm-calculation:
#'
#' SantaLucia, J, et al. (1996)
#' Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability.
#' Biochemistry, 35: 3555-3562 (Formula and salt correction method
#' are from here)
#'
#' Allawi, H. & SantaLucia, J. (1997)
#' Thermodynamics and NMR of Internal G-T Mismatches in DNA.
#' Biochemistry, 36, 34: 10581–10594
#' (Duplex initiation parameters are from here)
#'
#' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @export
getOligos <- function(x,
                      length = 18:22,
                      maxGapFrequency = 0.1,
                      maxDegeneracy = 4,
                      gcClamp = TRUE,
                      avoid3EndRuns = TRUE,
                      avoid5EndG = FALSE,
                      gcRange = c(0.45, 0.55),
                      tmRange = c(48, 70),
                      concOligo = 5e-07,
                      concNa = 0.05,
                      showAllVariants = TRUE
                      ) {
  if (!is.logical(showAllVariants)) {
    stop("'showAllVariants' must be set to 'TRUE' or 'FALSE'.", call. = FALSE)
  }
  allOligos <- purrr::map_dfr(length, function(i) {
    oligos <- .generateOligos(
      x, oligoLength = i, maxGapFrequency = maxGapFrequency,
      maxDegeneracy = maxDegeneracy
    )
    oligos <- .exclude(oligos)
    oligos <- .addGcContent(oligos, gcRange = gcRange)
    oligos <- .addTm(
      oligos, concOligo = concOligo, concNa = concNa, tmRange = tmRange
    )
    oligos <- .addReverseComplement(oligos)
    oligos <- .filterOligos(
      oligos, gcClamp = gcClamp, avoid5EndG = avoid5EndG,
      avoid3EndRuns = avoid3EndRuns
    )
    oligos
  })
  if (nrow(allOligos) == 0L)
    stop("No oligos were found.", call. = FALSE)
  if (showAllVariants == TRUE) {
    allOligos <- .expandOligos(
      allOligos, concOligo = concOligo, concNa = concNa
    )
  }
 allOligos
}

# Helpers =====================================================================

#' Split sequence
#'
#' @param x A character vector of length one.
#'
#' @return A character vector of length \code{nchar(x)}.
#'
#' @keywords internal
#'
#' @noRd
.splitSequence <- function(x) {
  x <- unlist(strsplit(x, split = ""), use.names = FALSE)
  x
}

#' Divide a DNA sequence into n-sized chunks
#'
#' \code{.getNmers} divides a character vector into chunks of size \code{n}.
#'
#' @param x A character vector.
#'
#' @param n Chunk-size (a positive integer).
#'
#' @return A character vector.
#'
#' @keywords internal
#'
#' @noRd
.getNmers <- function(x, n) {
  begin <- seq_len(length(x) - n + 1)
  end <- begin + n - 1
  nmer <- purrr::map_chr(
    begin, ~ paste(x[begin[[.x]]:end[[.x]]], collapse = "")
  )
  nmer
}

#' Calculate running sums
#'
#' \code{.runningSum} calculates 'running' sums of a numeric vector. Each
#' sum is calculated in a size of \code{n}, in steps of 1 (i.e., if
#' \code{n = 20}, the sum will be calculated from element 1 to 20,
#' then from element 2 to 21, then from element 3 to 22, etc.)
#'
#' @param x A numeric vector.
#'
#' @param n The size of each sum that is to be calculated,
#'
#' @return The running sums of \code{x}
#' (a numeric vector of length \code{length(x) - n + 1}).
#'
#' @keywords internal
#'
#' @noRd
.runningSum <- function(x, n) {
  cumul <- c(0, cumsum(x))
  runsum <- cumul[seq(n + 1, length(cumul))] - cumul[seq_len(length(cumul) - n)]
  runsum
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{.countDegeneracy} returns the number of unique variants of
#' a DNA sequence with degenerate bases.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#'
#' @return The number of unique sequences of x (an integer).
#'
#' @keywords internal
#'
#' @noRd
.countDegeneracy <- function(x) {
  x <- toupper(x)
  x <- .splitSequence(x)
  nNucleotides <- degeneracyLookup[x]
  degeneracy <- prod(nNucleotides)
  degeneracy
}

#' Generate oligos of a specific length
#'
#' @inheritParams getOligos
#'
#' @return A tibble with oligos.
#'
#' @keywords internal
#'
.generateOligos <- function(x,
                            oligoLength = 20,
                            maxGapFrequency = 0.1,
                            maxDegeneracy = 4) {
  if (!(min(oligoLength) >= 14 && max(oligoLength) <= 30)) {
    stop("'oligoLength' must be between 14 and 30.", call. = FALSE)
  }
  if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
    stop("'maxGapFrequency' must be between 0 and 1.", call. = FALSE)
  }
  if (!(maxDegeneracy >= 1 && maxDegeneracy <= 32)) {
    stop("'maxDegeneracy' must be between 1 and 32.", call. = FALSE)
  }
  Majority <- .getNmers(x$Majority, n = oligoLength)
  IUPAC <- .getNmers(x$IUPAC, n = oligoLength)
  Degeneracy <- as.integer(purrr::map_dbl(IUPAC, ~ .countDegeneracy(.x)))
  Begin <- seq_along(Majority)
  End <- as.integer(seq_along(Majority) + oligoLength - 1)
  Length <- oligoLength
  gapBin <- ifelse(x$Gaps > maxGapFrequency, 1L, 0L)
  gapPenalty <- .runningSum(gapBin, n = oligoLength)
  oligos <- tibble::tibble(
    Begin, End, Length, Majority, IUPAC,
    Degeneracy, gapPenalty
  )
  uniqueOligos <- match(oligos$Majority, unique(oligos$Majority))
  oligos <- oligos[uniqueOligos, ]
  oligos <- oligos[oligos$gapPenalty == 0, ]
  oligos <- oligos[oligos$Degeneracy <= maxDegeneracy, ]
  oligos <- dplyr::select(oligos, -gapPenalty)
  oligos
}

#' Exclude non optimal oligos
#'
#' \code{.exclude} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}.
#'
#' @param x A tibble with oligos.
#'
#' @details
#' An oligo is excluded if:
#' - It has more than three runs of the same dinucleotide (e.g. "TATATATA")
#' - It has more than four runs of the same nucleotide (e.g. ("AAAAA"))
#'
#' @return
#' A tibble where non optimal oligos have been excluded.
#'
#' @keywords internal
#'
#' @noRd
.exclude <- function(x) {
  dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})"
  mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
  x <- x[!grepl(dinucleotideRepeats, x$Majority), ]
  x <- x[!grepl(mononucleotideRepeates, x$Majority), ]
  x
}

#' Get the reverse complement of a DNA sequence
#'
#' \code{.reverseComplement} finds the reverse complement of a DNA sequence.
#'
#' @param x A DNA sequence (a character vector of length one).
#'
#' @return The reverse complement. Non valid bases will return as \code{NA}.
#'
#' @keywords internal
#'
#' @noRd
.reverseComplement <- function(x) {
  x <- toupper(x)
  x <- strsplit(x, split = "")
  complement <- complementLookup[unlist(x)]
  complement <- unname(complement)
  rc <- rev(complement)
  rc <- paste(rc, collapse = "")
  rc
}

#' Add reverse complement to generated oligos
#'
#' @param x A tibble with oligos.
#'
#' @return A tibble.
#'
#' @keywords internal
#'
#' @noRd
.addReverseComplement <- function(x) {
  Majority_RC <- purrr::map_chr(x$Majority, ~.reverseComplement(.x))
  IUPAC_RC <- purrr::map_chr(x$IUPAC, ~.reverseComplement(.x))
  x <- tibble::add_column(
    x, Majority_RC, IUPAC_RC, .before = "Degeneracy"
  )
  x
}

#' Identify GC clamp
#'
#' @param x One or more oligo sequences (a character vector).
#'
#' @return
#' A character vector, where oligo sequences without GC clamp have
#' been replaced with NA.
#'
#' @keywords internal
#'
#' @noRd
.getOligosWithGcClamp <- function(x) {
  ends <- purrr::map(x, ~.splitSequence(.x))
  ends <- purrr::map(ends, ~.x[(length(.x) - 4):length(.x)])
  gc <- purrr::map_dbl(ends, ~.gcContent(.x))
  x[gc >= 4/5] <- NA
  x[gc <= 1/5] <- NA
  x
}

#' Replace unwanted oligos with NA
#'
#' @param x A tibble with oligos (with reverse complements).
#'
#' @inheritParams getOligos
#'
#' @return A tibble where unwanted oligos are replaced with NA.
#'
#' @keywords internal
#'
#' @noRd
.filterOligos <- function(x,
                          gcClamp = TRUE,
                          avoid5EndG = FALSE,
                          avoid3EndRuns = TRUE) {
  if (any(!is.logical(
    c(gcClamp, avoid5EndG, avoid3EndRuns)
  ))) {
    stop(
      "'gcClamp', 'avoid5EndG',  and 'avoid3EndRuns' must be set to \n
      'TRUE' or 'FALSE'", call. = FALSE
    )
  }
  if (gcClamp == TRUE) {
    x$Majority <- .getOligosWithGcClamp(x$Majority)
    x$Majority_RC <- .getOligosWithGcClamp(x$Majority_RC)
  }
  if (avoid5EndG == TRUE) {
    x$Majority[grepl("^G", x$Majority)] <- NA
    x$Majority_RC[grepl("^G", x$Majority_RC)] <- NA
  }
  if (avoid3EndRuns == TRUE) {
    x$Majority[grepl("([A-Z])\\1\\1$", x$Majority)] <- NA
    x$Majority_RC[grepl("([A-Z])\\1\\1$", x$Majority_RC)] <- NA
  }
  x$IUPAC[is.na(x$Majority)] <- NA
  x$IUPAC_RC[is.na(x$Majority_RC)] <- NA
  invalidOligos <- is.na(x$Majority) & is.na(x$Majority_RC)
  x <- x[!invalidOligos, ]
  x
}

#' Calculate GC content of a DNA sequence
#'
#' \code{.gcContent} finds the GC content of a DNA sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#'
#' @return The GC content of x.
#'
#' @keywords internal
.gcContent <- function(x) {
  x <- toupper(x)
  x <- .splitSequence(x)
  gcCount <- length(which(x == "C" | x == "G"))
  totalCount <- length(which(x == "A" | x == "C" | x == "G" | x == "T"))
  gc <- gcCount / totalCount
  gc
}

#' Add GC content to generated oligos
#'
#' @param x A tibble with oligos.
#'
#' @inheritParams getOligos
#'
#' @return A tibble.
#'
#' @keywords internal
#'
#' @noRd
.addGcContent <- function(x, gcRange = c(0.45, 0.65)) {
  if (!(min(gcRange) >= 0 && max(gcRange) <= 1)) {
    stop(
      "'gcRange' must be from 0 to 1, e.g. c(0.45, 0.65).", call. = FALSE
    )
  }
  GC_majority <- purrr::map_dbl(x$Majority, ~ .gcContent(.x))
  x <- tibble::add_column(
    x, GC_majority, .after = "Degeneracy"
  )
  x <- x[x$GC_majority >= min(gcRange) & x$GC_majority <= max(gcRange), ]
  x
}

#' Split a DNA sequence into nearest neighbors
#'
#' \code{.nnSplit} splits an oligo sequence into nearest neighbors
#' (for calculation of deltaG, deltaH and Tm)
#'
#' @param x a DNA sequence with at least two bases, e.g. 'CTTA'
#' (a character vector of length one).
#'
#' @return The nearest neighbors of x (a character vector).
#'
#' @keywords internal
#'
#' @noRd
.nnSplit <- function(x) {
  x <- .splitSequence(x)
  from <- (seq_along(x) - 1)[-1]
  to <- seq_along(x)[-1]
  nn <- purrr::map2_chr(from, to, function(i, j) {
    paste(x[i:j], collapse = "")
  })
  nn
}

#' Calculate dH or dS of nearest neighbors using lookup tables
#'
#' \code{.getNnTableValues} finds the corresponding dH or dS values
#' for nearest-neighbor pairs.
#'
#' @param x A matrix (of type character) with nearest-neighbor pairs
#' of DNA sequences (e.g. \code{c('CT', 'TT', 'TA')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#'
#' @return The corresponding values for dH or dS (in cal/M).
#'
#' @keywords internal
#'
#' @noRd
.getNnTableValues <- function(x, table = "dH") {
  if (table == "dH") {
    selected <- nnLookup$dH
  } else {
    selected <- nnLookup$dS
  }
  matching <- selected[match(x, nnLookup$bases)]
  if (is.null(ncol(x))) {
    matching
  } else {
    matrix(matching, ncol = ncol(x), byrow = FALSE)
  }
}

#' Initiation of DNA sequences for Tm calculation
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @return The initiaion values for x.
#'
#' @keywords internal
#'
#' @noRd
.init3End <- function(x) {
  if (grepl("(t|a)$", x)) {
    c(H = 2.3 * 1000, S = 4.1)
  } else {
    c(H = 0.1 * 1000, S = -2.8)
  }
}

#' @describeIn .init3End
.init5End <- function(x) {
  if (grepl("^(t|a)", x)) {
    c(H = 2.3 * 1000, S = 4.1)
  } else {
    c(H = 0.1 * 1000, S = -2.8)
  }
}

#' Melting temperature
#'
#' \code{.tm} calculates the melting temperature of one or
#' more perfectly matching DNA duplexes (i.e. oligo-target duplexes),
#' using the nearest neighbor method.
#'
#' @param oligos One or more DNA sequences (a character vector).
#'
#' @inheritParams getOligos
#'
#' @return The melting temperature(s) of x.
#'
#' @keywords internal
#'
#' @noRd
.tm <- function(oligos, concOligo = 5e-07, concNa = 0.05) {
  if (!is.double(concOligo) || concOligo < 2e-07 || concOligo > 2.0e-06) {
    stop("'concOligo' must be from
           0.2e-07 M (20 nM) to 2e-06 M (2000 nM).", call. = FALSE)
  }
  if (!is.double(concNa) || concNa < 0.01 || concNa > 1) {
    stop("'concNa' must be from 0.01 to 1 M.", call. = FALSE)
  }
  oligos <- toupper(oligos)
  # Find initiation values
  initH <- purrr::map_dbl(oligos, function(x) {
    .init5End(x)[["H"]] + .init3End(x)[["H"]]
  })
  initS <- purrr::map_dbl(oligos, function(x) {
    .init5End(x)[["S"]] + .init3End(x)[["S"]]
  })
  # Split to nearest neighbors
  nn <- purrr::map(oligos, .nnSplit)
  # Check oligo length
  oligoLength <- purrr::map_int(nn, length)
  # I made a matrix based tm-calculation,
  # which means that all oligos must be of the same length
  if (length(unique(oligoLength)) != 1) {
    stop("All oligos must be of equal length.", call. = FALSE)
  }
  # Find nearest neighbor values for dH and dS
  nn <- do.call("rbind", nn)
  dhResult <- .getNnTableValues(nn, "dH")
  dsResult <- .getNnTableValues(nn, "dS")
  # Sum dH and dS
  sumdH <- rowSums(dhResult) + initH
  sumdS <- rowSums(dsResult) + initS
  # Correct delta S for salt conc.
  N <- nchar(oligos[[1]]) - 1 # Number of phosphates
  sumdS <- sumdS + 0.368 * N * log(concNa)
  tm <- sumdH / (sumdS + gasConstant * log(concOligo))
  tm <- tm - 273.15
  tm
}

#' Add Tm to oligos
#'
#' @param x A tibble with oligos.
#'
#' @inheritParams getOligos
#'
#' @keywords internal
#'
#' @noRd
.addTm <- function(x,
                   concOligo = 5e-07,
                   concNa = 0.05,
                   tmRange = c(55, 75)
                   ) {
  if (!(min(tmRange) >= 20 && max(tmRange) <= 90)) {
    stop(
      "'tmRange' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE
    )
  }
  Tm_majority <- .tm(x$Majority, concOligo = concOligo, concNa = concNa)
  x <- tibble::add_column(
    x, Tm_majority, .after = "Degeneracy"
  )
  x <- x[x$Tm_majority >= min(tmRange) & x$Tm_majority <= max(tmRange), ]
  x
}

#' Find all variants of a DNA sequence
#'
#' @param x
#' A DNA sequence, i.e. a character vector of length one,
#' e.g. "ACTTTGR".
#'
#' @return
#' All variants of \code{x}. A character vector.
#'
#' @keywords internal
#'
#' @noRd
.expandDegenerates <- function(x) {
  x <- .splitSequence(x)
  expanded <- purrr::map(x, function(i) {
    allBases <- unname(degenerateLookup[[i]])
    allBases <- unlist(strsplit(allBases, split = ","))
    allBases
  })
  expanded <- expand.grid(
    expanded[seq_along(expanded)], stringsAsFactors = FALSE
  )
  expanded <- purrr::map(
    seq_len(nrow(expanded)), ~paste(expanded[.x, ], collapse = "")
  )
  expanded <- unlist(expanded, use.names = FALSE)
  expanded
}

#' Add info on sequence, GC-content and Tm of all oligo variants
#'
#' @param x A tibble with oligos.
#'
#' @inheritParams getOligos
#'
#' @return A tibble with oligos, including all variants.
#'
#' @keywords internal
#'
#' @noRd
.expandOligos <- function(x, concOligo = 5e-7, concNa = 0.05) {
  All <- purrr::map(x$IUPAC, function(x) {
    if (is.na(x)) "" else .expandDegenerates(x)
  })
  x <- tibble::add_column(x, All, .before = "Degeneracy")
  All_RC <- purrr::map(x$IUPAC_RC, function(x) {
    if (is.na(x)) "" else .expandDegenerates(x)
  })
  x <- tibble::add_column(x, All_RC, .after = "All")
  GC_all <- purrr::map(seq_len(nrow(x)), function(i) {
      toCalculate <- ifelse(!is.na(x[i, ]$Majority), x[i, ]$All, x[i, ]$All_RC)
      toCalculate <- unlist(toCalculate)
      purrr::map_dbl(toCalculate, ~.gcContent(.x))
  })
  Tm_all <- purrr::map(seq_len(nrow(x)), function(i) {
    toCalculate <- ifelse(!is.na(x[i, ]$Majority), x[i, ]$All, x[i, ]$All_RC)
    toCalculate <- unlist(toCalculate)
    purrr::map_dbl(toCalculate, ~.tm(.x))
  })
  x <- tibble::add_column(x, Tm_all, GC_all)
  x
}
