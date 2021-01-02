
#data("exampleRprimerProfile")
x <- exampleRprimerProfile
ols <- getOligos(x, maxDegeneracyPrimer = 1, maxDegeneracyProbe = 1)

# argument rbind true, false ........ # detect repeats if only one variant...
# arg 1 is not a vector order????
# if 1 degenerate only
# if no primers found?

# correct tm in .filterProbes, # return a function???????????????
# document getOligos and exdata, and write tests
# then getAssays, then plot, then vignette

# Exported =====================================================================

#' Get oligos
#'
#' \code{getOligos()} identifies oligos (primers and probes)
#' from an \code{RprimerProfile} object.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency (for both primers and probes).
#' A number [0, 1]. Defaults to 0.1.
#'
#' @param lengthPrimer
#' Primer length. A numeric vector [14, 30].
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracyPrimer
#' Maximum number of variants of each primer. A number [1, 32]. Defaults to 4.
#'
#' @param gcClampPrimer
#' \code{TRUE} or \code{FALSE}.
#' If primers with no GC-clamp
#' should be replaced with \code{NA}.
#' Defaults to \code{TRUE}. A GC-clamp
#' is identified as two to three G or
#' C:s within the last five bases (3'-end) of the oligo.
#'
#' @param avoidTreeEndRunsPrimer
#' \code{TRUE} or \code{FALSE}.
#' If primers with more than two runs
#' of the same nucleotide at the 3' end should be considered as invalid.
#'
#' @param minEndIdentityPrimer
#' Optional. If specified, a number [0, 1]. The minimum allowed identity
#' at the 3'-end of the primer (i.e. the last five bases).
#'
#' @param gcRangePrimer
#' GC-content range for primers (proportion, not percent).
#' A numeric vector [0, 1]. Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRangePrimer
#' Tm range for primers.
#' A numeric vector [30, 90]. Defaults to \code{c(55, 65)}.
#'
#' @param concPrimer
#' Primer concentration in nM, for Tm calculation. A number
#' [20, 2000] Defaults to 500 nM.
#'
#' @param probe
#' If probes should be designed as  well. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param lengthProbe
#' Probe length. A numeric vector [14, 30].
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracyProbe
#' Maximum number of variants of each probe. A number [1, 32]. Defaults to 4.
#'
#' @param avoidFiveEndGProbe
#' \code{TRUE} or \code{FALSE}. If probes with G
#' at the 5'-end should be excluded. Defaults to \code{TRUE}.
#'
#' @param gcRangeProbe
#' GC-content range for probes (proportion, not %). A numeric vector [0, 1].
#' Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRangeProbe
#' Tm range for probes.
#' A numeric vector [30, 90]. Defaults to \code{c(55, 65)}.
#'
#' @param concProbe
#' Primer concentration in nM, for Tm calculation. A numeric vector
#' [20, 2000] Defaults to 250 nM.
#'
#' @param concNa
#' The sodium ion concentration in the PCR reaction in M, for Tm calculation.
#' A numeric vector [0.01, 1]. Defaults to 0.05 M (50 mM).
#'
#' @section Excluded oligos:
#' \code{getOligos()} excludes:
#'
#' \itemize{
#' \item Majority oligos with more than than three consecutive runs of ##############################
#' the same dinucleotide (e.g. "TATATATA")
#' \item Majority oligos with more than four consecutive runs of
#' the same nucleotide  (e.g. "AAAAA")
#'
#' @section Tm-calculation:
#'
#' Melting temperatures are calculated using SantaLucia's nearest-neighbor
#' method, with the following assumptions:
#'
#' \itemize{
#'   \item Oligos are not expected to be self-complementary (i.e. no symmetry
#'   correction is done).
#'   \item The oligo concentration is assumed to be much higher
#'   than the target concentration.
#' }
#'
#' See references for table values and equations.
#' Table values can also be found by calling \code{rprimer:::lookup$nn}.
#'
#' @return
#' An \code{RprimerOligo} object.
#' An error message will return if no oligos are found.
#'
#' The object contains the following columns:
#'
#' \describe{
#'   \item{type}{Type of oligo, primer or probe.}
#'   \item{start}{Position where the oligo starts.}
#'   \item{end}{Position where the oligo ends.}
#'   \item{length}{Length of the oligo.}
 ########################################################################################
#'   \item{identity}{Average identity score of the oligo.}
#'   \item{iupac}{IUPAC sequence (i.e. with degenerate bases).}
#'   \item{iupacRc}{IUPAC sequence, reverse complement.}
#'   \item{degeneracy}{Number of variants of the degenerate oligo.}
#'   \item{sequence}{Lists with all sequence variants of the oligos.}
#'   \item{sequenceRc}{Lists with all sequence variants of the oligos, reverse
#'   complements.}
#'   \item{gcContent}{Lists with the GC content of all
#'   sequence variants of the oligos.}
#'   \item{tm}{Lists with the Tm of all sequence variants of the oligos.}
#'   \item{alignmentStart}{Start position of the input consensus profile.}
#'   \item{alingnmentEnd}{End position of the input consensus profile.}
#' }
#'
#' @examples
#' ## Load example data
#' data("exampleRprimerProfile")
#' target <- exampleRprimerProfile
#'
#' ## Design primers and probes with default values
#' getOligos(target)
#'
#' ## Design primers and probes only within a specific region of interest
#' roi <- target[target$position >= 5000 & target$position <= 6000, ]
#' getOligos(roi)
#'
#' ## Design primers only
#' getOligos(roi, probe = FALSE)
#'
#' @references
#' SantaLucia Jr, J., & Hicks, D. (2004).
#' The thermodynamics of DNA structural motifs.
#' Annu. Rev. Biophys. Biomol. Struct., 33, 415-440.
#'
#' @export
getOligos <- function(x,
                      maxGapFrequency = 0.1,
                      lengthPrimer = 18:22,
                      maxDegeneracyPrimer = 4,
                      gcClampPrimer = TRUE,
                      avoidThreeEndRunsPrimer = TRUE,
                      minEndIdentityPrimer = 0.98,
                      gcRangePrimer = c(0.45, 0.55),
                      tmRangePrimer = c(55, 65),
                      concPrimer = 500,
                      probe = TRUE,
                      lengthProbe = 18:22,
                      maxDegeneracyProbe = 4,
                      avoidFiveEndGProbe = TRUE,
                      gcRangeProbe = c(0.45, 0.55),
                      tmRangeProbe = c(55, 70),
                      concProbe = 250,
                      concNa = 0.05) {
  if (!methods::is(x, "RprimerProfile")) {
    stop("'x' must be an RprimerProfile object.")
  }
  if (nrow(x) < max(c(lengthPrimer, lengthProbe))) {
    stop(paste(
      "In order to search for oligos, the number of rows in 'x'
        must be at least", max(c(lengthPrimer, lengthProbe), ".")
    ),
    call. = FALSE
    )
  }
  if (!(min(lengthPrimer) >= 14 && max(lengthPrimer) <= 30)) {
    stop("'lengthPrimer' must be from 14 to 30.", call. = FALSE)
  }
  if (!(min(lengthProbe) >= 14 && max(lengthProbe) <= 30)) {
    stop("'lengthProbe' must be from 14 to 30.", call. = FALSE)
  }
  if (!(maxDegeneracyPrimer >= 1 && maxDegeneracyPrimer <= 32)) {
    stop("'maxDegeneracyPrimer' must be from 1 to 32.", call. = FALSE)
  }
  if (!(maxDegeneracyProbe >= 1 && maxDegeneracyProbe <= 32)) {
    stop("'maxDegeneracyProbe' must be from 1 to 32.", call. = FALSE)
  }
  if (!is.logical(probe)) {
    stop("'probe' must be set to TRUE or FALSE", call. = FALSE)
  }
  if (probe) {
    lengthOligo <- c(lengthPrimer, lengthProbe)
    lengthOligo <- unique(lengthOligo)
    lengthOligo <- lengthOligo[order(lengthOligo)]
  } else {
    lengthOligo <- lengthPrimer
  }
  maxDegeneracy <- if (probe) {
    max(c(maxDegeneracyPrimer, maxDegeneracyProbe))
  } else {
    maxDegeneracyPrimer
  }
  oligos <- .designOligos(x,
                          maxGapFrequency = maxGapFrequency,
                          lengthOligo = lengthOligo,
                          maxDegeneracy = maxDegeneracy,
                          concOligo = concPrimer,
                          concNa = concNa)
  primers <- .filterPrimers(oligos,
                            lengthPrimer = lengthPrimer,
                            maxDegeneracyPrimer = maxDegeneracyPrimer,
                            minEndIdentityPrimer = minEndIdentityPrimer,
                            gcClampPrimer = gcClampPrimer,
                            avoidThreeEndRunsPrimer = avoidThreeEndRunsPrimer,
                            gcRangePrimer = gcRangePrimer,
                            tmRangePrimer = tmRangePrimer)
  if (nrow(primers) == 0L) {
    stop("No primers were found.", call. = FALSE)
  }
  if (probe) {
    probes <- .filterProbes(oligos,
                            lengthProbe = lengthProbe,
                            maxDegeneracyProbe = maxDegeneracyProbe,
                            avoidFiveEndGProbe = avoidFiveEndGProbe,
                            gcRangeProbe = gcRangeProbe,
                            tmRangeProbe = tmRangeProbe,
                            concProbe = concProbe)
    if (nrow(probes) == 0L) {
      stop("No probes were found.", call. = FALSE)
    }
    oligos <- list(primers, probes)
  } else {
    oligos <- primers
  }
  oligos
#  oligos <- .arrangeData(oligos)
#  RprimerOligo(oligos)
}

# Helpers/internal functions ===================================================

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
  unlist(strsplit(x, split = ""), use.names = FALSE)
}

#' Arrange a vector into a matrix of "n-mers"
#'
#' \code{.getNmers()} divides a vector into a matrix with \code{n} columns.
#' It is used to generate all possible oligos of a specific length
#' from a DNA sequence.
#'
#' Helper function to \code{.generateOligos()}.
#'
#' @param x A vector.
#'
#' @param n Length of each "mer" (a positive integer).
#'
#' @return A matrix with \code{n} columns. The first row will contain
#' \code{x[1:n]}, the second \code{x[2:(n + 1)]}, and so on.
#' In this way, the matrix will contain all possible consecutive
#' sequences of \code{x} of length \code{n}.
#'
#' @examples
#' .getNmers(c("A", "G", "T", "T", "C", "G"), n = 4)
#'
#' @keywords internal
#'
#' @noRd
.getNmers <- function(x, n) {
  start <- seq_len(length(x) - n + 1)
  end <- start + n - 1
  nmers <- lapply(start, function(i) x[start[i]:end[i]])
  do.call("rbind", nmers)
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{.countDegeneracy()} finds the number of unique variants of
#' a DNA sequence with (or without) degenerate bases.
#'
#' Helper function to \code{.generateOligos()}.
#'
#' @param x A DNA sequence (a character vector).
#'
#' @return The number of sequence variants of x (an integer).
#'
#' @examples
#' .countDegeneracy(c("A", "R", "T", "T", "N", "G"))
#'
#' @keywords internal
#'
#' @noRd
.countDegeneracy <- function(x) {
  nNucleotides <- lookup$degeneracy[x]
  prod(nNucleotides)
}

#' Generate oligos of a specific length
#'
#' \code{.generateOligos()} is the first step of the oligo-design process.
#' It finds all possible oligos of a specific length from an
#' \code{RprimerProfile} object, and returns a list containing
#' start and end position,
#' length, IUPAC sequence (the oligo DNA sequence with wobble bases), degeneracy
#' (number of variants of each oligo), maximum gap frequency,
#' mean overall identity, and minimum 3'-end identity at
#' both forward and reverse direction (the 3'-end is here seen as the last five
#' bases of the oligo).
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param lengthOligo An integer [14, 30], defaults to 20.
#'
#' @return A list with all candidate oligos.
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateOligos(exampleRprimerProfile, lengthOligo = 18)
#'
#' @keywords internal
#'
#' @noRd
.generateOligos <- function(x, lengthOligo = 20) {
  iupacSequence <- .getNmers(x$iupac, lengthOligo)
  start <- seq_len(nrow(iupacSequence)) + min(x$position) - 1
  end <- seq_len(nrow(iupacSequence)) + lengthOligo - 1 + min(x$position) - 1
  length <- rep(lengthOligo, nrow(iupacSequence))
  degeneracy <- apply(iupacSequence, 1, .countDegeneracy)
  gapFrequency <- apply(.getNmers(x$gaps, lengthOligo), 1, max)
  identity <- .getNmers(x$identity, lengthOligo)
  endIdentityFwd <- apply(
    identity[, (ncol(identity) - 4):ncol(identity)], 1, min
  )
  endIdentityRev <- apply(identity[, seq_len(5)], 1, min)
  identity <- rowMeans(identity)
  alignmentStart <- rep(min(x$position, na.rm = TRUE), nrow(iupacSequence))
  alignmentEnd <- rep(max(x$position, na.rm = TRUE), nrow(iupacSequence))
  oligos <- list(
    "start" = start, "end" = end, "length" = length,
    "iupacSequence" = iupacSequence, "degeneracy" = degeneracy,
    "gapFrequency" = gapFrequency, "identity" = identity,
    "endIdentityFwd" = endIdentityFwd, "endIdentityRev" = endIdentityRev,
    "alignmentStart" = alignmentStart, "alignmentEnd" = alignmentEnd
  )
  oligos
}

#' Remove oligos with too high gap frequency and degeneracy
#'
#' \code{.filterOligos()} is the second step of the oligo-design process. It
#' removes oligos with gap frequency and degeneracy above the user-specified
#' thresholds.
#'
#' @param x An output from \code{.generateOligos()}.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency. A number [0, 1], defaults to 0.1.
#'
#' @param maxDegeneracy
#' Maximum allowed number of variants of each oligo.
#' An integer, defaults to 4.
#'
#' @return
#' A list with the same structure as \code{x}, but where invalid oligos have
#' been removed.
#'
#' @keywords internal
#'
#' @noRd
.filterOligos <- function(x, maxGapFrequency = 0.1, maxDegeneracy = 4) {
  invalid <- unique(c(
    which(x$degeneracy > maxDegeneracy),
    which(x$gapFrequency > maxGapFrequency)
  ))
  if (length(invalid > 0L)) {
    lapply(x, function(x) {
      if (is.matrix(x)) x[-invalid, , drop = FALSE] else x[-invalid]
    })
   } else {
    x
  }
}

#' Get all variants of a DNA sequence with wobble bases
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A DNA sequence, i.e. a character vector.
#'
#' @return All sequence variants of \code{x}. A matrix.
#'
#' @examples
#' .expandDegenerates(c("A", "R", "T", "T", "N", "G"))
#'
#' @keywords internal
#'
#' @noRd
.expandDegenerates <- function(x) {
  bases <- lapply(x, function(i) {
    all <- unname(lookup$degenerates[[i]])
    unlist(strsplit(all, split = ","))
  })
  allVariants <- expand.grid(bases[seq_along(bases)], stringsAsFactors = FALSE)
  allVariants <- as.matrix(allVariants)
  colnames(allVariants) <- NULL
  allVariants
}

#' Turn a list of oligos into a matrix with oligos
#'
#' \code{.makeOligoMatrix()} is part of a "workaround" to avoid loops when
#' finding e.g. the reverse complement, the presence of GC-clamp etc.
#' of many DNA sequences (oligos).
#' It takes a list of DNA sequences as input (where each element contains a
#' matrix with all variants of a specific degenerate oligo)
#' and returns a single matrix, where the belonging of each oligo sequence
#' is identified by a rowname ID.
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A list with oligos.
#'
#' @return A matrix with all sequence variants, with oligo IDs as rownames.
#'
#' @keywords internal
#'
#' @noRd
.makeOligoMatrix <- function(x) {
  degeneracy <- vapply(x, nrow, integer(1))
  id <- lapply(seq_along(degeneracy), function(x) rep(x, degeneracy[x]))
  id <- unlist(id)
  x <- do.call("rbind", x)
  rownames(x) <- id
  x
}

#' Find the reverse complement of a DNA sequence
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A vector or matrix with DNA sequence(s).
#'
#' @return The reverse complement of \code{x}, a matrix with the same
#' dimension as \code{x}, and with the same rownames.
#'
#' @examples
#' .reverseComplement(c("A", "R", "T", "T", "N", "G"))
#'
#' @keywords internal
#'
#' @noRd
.reverseComplement <- function(x) {
  if (!is.matrix(x)) x <- t(matrix(x))
  rev <- x[, rev(seq_len(ncol(x)))]
  rc <- unname(lookup$complement[rev])
  rc <- matrix(rc, ncol = ncol(x), byrow = FALSE)
  rownames(rc) <- rownames(rev)
  rc
}

#' Identify oligos with a GC-clamp
#'
#' \code{.detectGcClamp()} detects the presence of a GC-clamp on oligos
#' (good to have on primers).
#' A GC-clamp is identified as two to three G or C:s at
#' the 3'-end (= the last five bases).
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A numeric vector or matrix, where each row corresponds to a specific oligo.
#' In this matrix, 1 corresponds to G or C, and 0 corresponds to A or T.
#'
#' @param fwd
#' If the check should be done in forward direction,
#' defaults to \code{TRUE}. If \code{FALSE}, the check will be performed in
#' reverse direction.
#'
#' @return
#' A logical vector of length \code{nrow(x)}, where \code{TRUE} indicates the
#' presence of a GC-clamp.
#'
#' @examples
#' seq <- c("A", "C", "G", "G", "T", "T", "A", "A")
#' gc <- ifelse(seq == "C" | seq == "G", 1, 0)
#' .detectGcClamp(gc, fwd = TRUE)
#' .detectGcClamp(gc, fwd = FALSE)
#'
#' @keywords internal
#'
#' @noRd
.detectGcClamp <- function(x, fwd = TRUE) {
  if (!is.matrix(x)) x <- t(matrix(x))
  gc <- if (fwd) {
    end <- x[, seq((ncol(x) - 4), ncol(x)), drop = FALSE]
    ifelse(rowSums(end) >= 2 & rowSums(end) <= 3, TRUE, FALSE)
    } else {
      end <- x[, seq_len(5),  drop = FALSE]
      ifelse(5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3, TRUE, FALSE)
    }
  gc
}

#' Identify oligos with runs of the same nucleotide at the 3' end
#'
#' \code{.detectThreeEndRuns()} detects if the same nucleotide is repeated at
#' at least 3 times at the terminal 3'-end of an oligo (e.g. "AAA")
#' (bad to have on primers).
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A character vector or matrix, where each row corresponds to a specific oligo.
#'
#' @param fwd
#' If the check should be done in forward direction,
#' defaults to \code{TRUE}. If \code{FALSE}, the check will be performed in
#' reverse direction.
#'
#' @return
#' A logical vector of length \code{nrow(x)}, where \code{TRUE} indicates the
#' presence of a 3'-end run.
#'
#' @examples
#' seq <- c("A", "C", "G", "G", "T", "T", "A", "A")
#' .detectThreeEndRuns(seq, fwd = TRUE)
#'
#' @keywords internal
#'
#' @noRd
.detectThreeEndRuns <- function(x, fwd = TRUE) {
  if (!is.matrix(x)) x <- t(matrix(x))
  end <- if (fwd) {
    x[, seq((ncol(x) - 2), ncol(x)), drop = FALSE]
  } else {
    x[, seq_len(3), drop = FALSE]
  }
  apply(end, 1, function(x) {
    all(x == "A") | all(x == "C") | all(x == "T") | all(x == "G")
  })
}

#' Get all variants of oligos with degenerate bases
#'
#' \code{.getAllVariants()} is the third step of the oligo-design process.
#' It returns all sequence variants of each oligo, both in sense and anti-sense
#' (reverse complement) direction. For each sequence variant, information is
#' also provided on GC-content, the presence of a GC-clamp (good for primers),
#' terminal 3'-end run (bad for primers) and a terminal
#' five end G (good for probes), and melting temperature.
#'
#' Helper function to \code{.designOligos()},
#'
#' @param x An output from \code{.filterOligos()}.
#'
#' @param concOligo
#' Oligo concentration in nM (for tm-calculation).
#'
#' @param concNa
#' Sodium ion concentration in the PCR reaction in M (for tm-calculation).
#'
#' @return
#' A nested list with sequence, reverse complement, GC-content,
#' GC-clamp, 3' end runs, 5' end G and melting temperature of each oligo.
#'
#' @keywords internal
#'
#' @noRd
.getAllVariants <- function(x, concOligo = 500, concNa = 0.05) {
  all <- list()
  all$sequence <- apply(x$iupacSequence, 1, .expandDegenerates)
  ## some kind of workaround if there is only one variant of each oligo... #######################
  if (!is.list(all$sequence)) all$sequence <- list(t(all$sequence))
  all$sequence <- .makeOligoMatrix(all$sequence)
  all$sequenceRc <- .reverseComplement(all$sequence)
  gc <- ifelse(all$sequence == "C" | all$sequence == "G", 1, 0)
  n <- rowSums(ifelse(
    all$sequence == "A" | all$sequence == "C" |
      all$sequence == "G" | all$sequence == "T", 1, 0
  ))
  all$gcContent <- rowSums(gc)/n
  all$gcClampFwd <- .detectGcClamp(gc, TRUE)
  all$gcClampRev <- .detectGcClamp(gc, FALSE)
  all$threeEndRunsFwd <- .detectThreeEndRuns(all$sequence, TRUE)
  all$threeEndRunsRev <- .detectThreeEndRuns(all$sequence, FALSE)
  all$fiveEndGPlus <- ifelse(all$sequence[, 1] == "G", TRUE, FALSE)
  all$fiveEndGMinus <- ifelse(
    all$sequence[, ncol(all$sequence)] == "C", TRUE, FALSE
  )
  all$tm <- .tm(all$sequence, concOligo, concNa)
  all$sequence <- apply(all$sequence, 1, paste, collapse = "")
  all$sequenceRc <- apply(all$sequenceRc, 1, paste, collapse = "")
  lapply(all, function(x) unname(split(unname(x), f = as.integer(names(x)))))
}

#' Calculate mean values
#'
#' When all sequence variants of each oligo are generated with
#' \code{.getAllVariants()}, the next step is to compute mean values on
#' GC-content and melting temperature, and proportion of the sequence
#' variants holding a GC-clamp, terminal 3'-end run and terminal 5'-end G.
#'
#' Helper function to \code{.designOligos()}.
#'
#' @param x An output from \code{.getAllVariants()}.
#'
#' @return A data frame with mean values for all numerical and
#' logical variables. Each row corresponds to an oligo.
#'
#' @keywords internal
#'
#' @noRd
.getMeanValues <- function(x) {
  remove <- c("sequence", "sequenceRc")
  x <- x[!(names(x) %in% remove)]
  means <- lapply(x, function(y) vapply(y, mean, double(1)))
  means <- do.call("cbind.data.frame", means)
  names(means) <- paste0(names(means), "Mean")
  means
}

#' Turn a list with oligos into a data frame
#'
#' Helper function to \code{.designOligos()}.
#'
#' @param x An output from \code{.generateOligos()} or \code{.filterOligos()}.
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
.makeOligoDf <- function(x) {
  x <- within(x, rm(gapFrequency))
  x$iupacSequenceRc <- .reverseComplement(x$iupacSequence)
  x$iupacSequence <- apply(x$iupacSequence, 1, paste, collapse = "")
  x$iupacSequenceRc <- apply(x$iupacSequenceRc, 1, paste, collapse = "")
  x <- do.call("cbind.data.frame", x)
  x[c(
    "start", "end", "length", "iupacSequence", "iupacSequenceRc", "identity",
    "degeneracy",
    "endIdentityFwd", "endIdentityRev", "alignmentStart", "alignmentEnd"
  )]
}

#' Turn a list with all oligo sequence variants into a data frame
#'
#' Helper function to \code{.designOligos()}.
#'
#' @param x An output from \code{.getAllVariants()}.
#'
#' @return
#' A data frame with sequence, reverse complement, GC content and tm for all
#' sequence variants.
#'
#' @keywords internal
#'
#' @noRd
.makeAllVariantDf <- function(x) {
  x <- x[c("sequence", "sequenceRc", "gcContent", "tm")]
  data.frame(do.call("cbind", x))
}

#' Design oligos
#'
#' Helper function to \code{getOligos()}.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency. A number [0, 1], defaults to 0.1.
#'
#' @param maxDegeneracy
#' Maximum allowed number of variants of each oligo.
#' An integer [1, 32], defaults to 4.
#'
#' @param concOligo
#' Oligo concentration in nM. A number
#' [20, 2000] Defaults to 250.
#'
#' @param concNa
#' Sodium ion concentration in the PCR reaction in M.
#' A number [0.01, 1]. Defaults to 0.05 M (50 mM).
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
.designOligos <- function(x,
                      lengthOligo = 18:22,
                      maxGapFrequency = 0.1,
                      maxDegeneracy = 4,
                      concOligo = 500,
                      concNa = 0.05) {
  allOligos <- lapply(lengthOligo, function(i) {
    oligos <- .generateOligos(x, lengthOligo = i)
    oligos <- .filterOligos(oligos,
                            maxGapFrequency = maxGapFrequency,
                            maxDegeneracy = maxDegeneracy)
    all <- .getAllVariants(oligos, concOligo = concOligo, concNa = concNa)
    means <- .getMeanValues(all)
    all <- .makeAllVariantDf(all)
    oligos <- .makeOligoDf(oligos)
    cbind(oligos, means, all)
  })
  allOligos <- do.call("rbind", allOligos)
  allOligos <- allOligos[c(
    "start", "end", "length", "iupacSequence", "iupacSequenceRc",
    "identity", "degeneracy",
    "gcContentMean", "tmMean", "sequence", "sequenceRc", "gcContent", "tm",
    "endIdentityFwd", "endIdentityRev",
    "gcClampFwdMean", "gcClampRevMean", "threeEndRunsFwdMean",
    "threeEndRunsRevMean", "fiveEndGPlusMean", "fiveEndGMinusMean",
    "alignmentStart", "alignmentEnd"
  )]
  allOligos
}

#' Get proportion of sequence variants within range
#'
#' \code{.getProportionInRange()} is used to find the proportion of sequence
#' variants of each oligo that falls within a specified range for
#' e.g. GC-content or tm.
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @param x A list.
#'
#' @param range The specified range. A numeric vector of length two.
#'
#' Helper function to \code{.designOligos()}.
#'
#' @keywords internal
#'
#' @noRd
.getProportionInRange <- function(x, range) {
  valid <- lapply(x, function(y) {
    ifelse(y >= min(range) & y <= max(range), TRUE, FALSE)
  })
  valid <- do.call("rbind", valid)
  rowMeans(valid)
}

#' Detect mono- and dinucleotide repeats
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @param x ###############################################################
#'
#' @return A vector
#'
#' @keywords internal
#'
#' @noRd
.detectRepeats <- function(x) {
  dinuclRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})" # sry ;)
  mononuclRepeates <- "([A-Z])\\1\\1\\1\\1"
  ok <- lapply(x, function(x) {
    ifelse(grepl(dinuclRepeats, x) | grepl(mononuclRepeates, x), TRUE, FALSE)
  })
  ok <- do.call("rbind", ok)
  rowMeans(ok) ###########################################################
}

#' Find oligos that pass the criteria for being a primer
#'
#' Helper function to \code{getOligos()}
#'
#' @param x An output from \code{.generateOligos()}
#'
#' @param lengthPrimer
#' Primer length. A numeric vector [14, 30].
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracyPrimer
#' Maximum allowed number of variants of each oligo.
#' An integer [1, 32], defaults to 4.
#'
#' @param minEndIdentityPrimer
#' A number. The minimum allowed identity
#' at the 3' end of the primer (i.e. the last five bases), defaults to 0.
#'
#' @param gcClampPrimer
#' \code{TRUE} or \code{FALSE}. If primers must have a GC-clamp.
#' Defaults to \code{TRUE}. A GC-clamp
#' is identified as two to three G or
#' C:s within the last five bases (3'-end) of the primer.
#'
#' @param avoidThreeEndRunsPrimer
#' \code{TRUE} or \code{FALSE}.
#' If primers with more than two runs
#' of the same nucleotide at the terminal 3'-end should be excluded.
#' Defaults to \code{TRUE}.
#'
#' @param gcRangePrimer
#' GC-content range for primers (proportion, not percent).
#' A numeric vector. Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRangePrimer
#' Tm range for primers.
#' A numeric vector. Defaults to \code{c(55, 65)}.
#'
#' @param acceptanceThreshold
#' Propoprtion of the sequence variants of each oligo that must
#' meet the design constraints. Defaults to 0.9.
#'
#' @keywords internal
#'
#' @noRd
.filterPrimers <- function(x,
                           lengthPrimer = 18:22,
                           maxDegeneracyPrimer = 4,
                           minEndIdentityPrimer = 0,
                           gcClampPrimer = TRUE,
                           avoidThreeEndRunsPrimer = TRUE,
                           gcRangePrimer = c(0.45, 0.55),
                           tmRangePrimer = c(55, 65),
                           acceptanceThreshold = 0.9) {

if (!(min(gcRangePrimer) >= 0 && max(gcRangePrimer) <= 1)) {
  stop(
    "'gcRangePrimer' must be from 0 to 1, e.g. c(0.45, 0.65).",
    call. = FALSE
  )
}
if (!(min(tmRangePrimer) >= 20 && max(tmRangePrimer) <= 90)) {
  stop("'tmRangePrimer' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
}
if (is.null(minEndIdentityPrimer)) minEndIdentityPrimer <- 0
if (minEndIdentityPrimer < 0 || minEndIdentityPrimer > 1) {
  stop(
    "'minEndIdentityPrimer' must be either NULL or from 0 to 1.", call. = FALSE
  )
}
if (any(!is.logical(c(gcClampPrimer, avoidThreeEndRunsPrimer)))) {
  stop(
    "'gcClampPrimer' and 'avoidThreeEndRunsPrimer'
    must be set to TRUE or FALSE",
    call. = FALSE
  )
}
if (!(min(gcRangePrimer) >= 0 && max(gcRangePrimer) <= 1)) {
  stop(
    "'gcRangePrimer' must be from 0 to 1, e.g. c(0.45, 0.65).",
    call. = FALSE
  )
}
if (!(min(tmRangePrimer) >= 20 && max(tmRangePrimer) <= 90)) {
  stop("'tmRangePrimer' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
}
  if (!gcClampPrimer) {
    x$gcClampFwdMean <- 1
    x$gcClampRevMean <- 1
  }
  if (!avoidThreeEndRunsPrimer) {
    x$threeEndRunsFwdMean <- 1
    x$threeEndRunsRevMean <- 1
  }
  x <- x[x$length >= min(lengthPrimer) & x$length <= max(lengthPrimer), ]
  x <- x[x$degeneracy <= maxDegeneracyPrimer, ]
  validFwd <- ifelse(
    x$endIdentityFwd >= minEndIdentityPrimer &
      x$gcClampFwdMean >= acceptanceThreshold &
      1 - x$threeEndRunsFwdMean >= acceptanceThreshold, TRUE, FALSE
  )
  validRev <- ifelse(
    x$endIdentityRev >= minEndIdentityPrimer &
      x$gcClampRevMean >= acceptanceThreshold &
      1 - x$threeEndRunsRevMean >= acceptanceThreshold, TRUE, FALSE
  )
  x <- cbind(x, validFwd, validRev)
  x <- x[x$validFwd | x$validRev, ]
  repeats <- .detectRepeats(x$sequence)
  valid <- ifelse(1 - repeats >= acceptanceThreshold, TRUE, FALSE)
  x <- x[valid, ]
  okForGc <- .getProportionInRange(x$gcContent, gcRangePrimer)
  x <- x[okForGc >= acceptanceThreshold, ]
  okForTm <- .getProportionInRange(x$tm, tmRangePrimer)
  x <- x[okForTm >= acceptanceThreshold, ]
  type <- rep("primer", nrow(x))
  cbind(type, x)
}

#' Find oligos that pass the criteria for being a probe
#'
#' Helper function to \code{getOligos()}
#'
#' @param x An output from \code{.generateOligos()}
#'
#' @param lengthProbe
#' Probe length. A numeric vector.
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracyProbe
#' Maximum allowed number of variants of each oligo.
#' An integer, defaults to 4.
#'
#' @param avoidFiveEndGProbe
#' If probes with a G at the terminal 5'-end should be avoided.
#' Defaults to TRUE.
#'
#' @param gcRangeProbe
#' GC-content range for probes (proportion, not percent).
#' A numeric vector. Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRangeProbe
#' Tm range for probes.
#' A numeric vector. Defaults to \code{c(55, 65)}.
#'
#' @param concProbe ######################################################
#'
#' @param acceptanceThreshold
#' Proportion of the sequence variants of each oligo that must
#' meet the design constraints. Defaults to 0.9.
#'
#' @keywords internal
#'
#' @noRd
.filterProbes <- function(x,
                          lengthProbe = 18:22,
                          maxDegeneracyProbe = 4,
                          avoidFiveEndGProbe = TRUE,
                          gcRangeProbe = c(0.45, 0.55),
                          tmRangeProbe = c(55, 65),
                          concProbe = 250,
                          acceptanceThreshold = 0.9) {
if (!is.logical(avoidFiveEndGProbe)) {
  stop("'avoidFiveEndGProbe' must be set to TRUE or FALSE", call. = FALSE)
}
if (!(min(gcRangeProbe) >= 0 && max(gcRangeProbe) <= 1)) {
  stop(
    "'gcRangeProbe' must be from 0 to 1, e.g. c(0.45, 0.65).",
    call. = FALSE
  )
}
if (!(min(tmRangeProbe) >= 20 && max(tmRangeProbe) <= 90)) {
  stop("'tmRangeProbe' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
}
  if (!avoidFiveEndGProbe) {
    x$fiveEndGPlusMean <- 1
    x$fiveEndGMinusMean <- 1
  }
  x <- x[x$length >= min(lengthProbe) & x$length <= max(lengthProbe), ]
  x <- x[x$degeneracy <= maxDegeneracyProbe, ]
  validFwd <- ifelse(
      1 - x$fiveEndGPlusMean >= acceptanceThreshold, TRUE, FALSE
  )
  validRev <- ifelse(
    1 - x$fiveEndGMinusMean >= acceptanceThreshold, TRUE, FALSE
  )
  x <- cbind(x, validFwd, validRev)
  x <- x[x$validFwd | x$validRev, ]
  repeats <- .detectRepeats(x$sequence)
  valid <- ifelse(1 - repeats >= acceptanceThreshold, TRUE, FALSE)
  x <- x[valid, ]
  okForGc <- .getProportionInRange(x$gcContent, gcRangeProbe)
  x <- x[okForGc >= acceptanceThreshold, ]
  # correct tm for probe conc .....######################################################################
  okForTm <- .getProportionInRange(x$tm, tmRangeProbe)
  x <- x[okForTm >= acceptanceThreshold, ]
  type <- rep("probe", nrow(x))
  x <- cbind(type, x)
  x
}

#' Arrange oligo data
#'
#' \code{.arrangeOligos()} drops unnecessary columns and sorts oligos based
#' on their start position.
#'
#' Helper function to \code{getOligos()}.
#'
#' @param x A data frame with oligos.
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
.arrangeData <- function(x) {
  keep <- c(
    "type", "validFwd", "validRev", "start", "end", "length", "iupacSequence",
    "iupacSequenceRc",
    "identity", "degeneracy", "gcContentMean", "tmMean", "sequence",
    "sequenceRc", "gcContent", "tm", "alignmentStart",
    "alignmentEnd"
  )
  x <- x[keep]
  x <- x[order(x$start), ]
  rownames(x) <- NULL
  x
}
