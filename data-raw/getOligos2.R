
data("exampleRprimerProfile") #
x <- exampleRprimerProfile


# exclude in .filterPrimers and .filterProbes
# correct tm
# document (also exdata) and write tests!!!!!!!!
# then getAssays, then plot, then vignette

#' Get oligos
#'
getOligos <- function(x) {
  if (!is.logical(avoidFiveEndG)) {
    stop("'avoidFiveEndG' must be set to TRUE or FALSE", call. = FALSE)
  }
  if (!(min(gcRange) >= 0 && max(gcRange) <= 1)) {
    stop(
      "'gcRange' must be from 0 to 1, e.g. c(0.45, 0.65).",
      call. = FALSE
    )
  }
  if (!(min(tmRange) >= 20 && max(tmRange) <= 90)) {
    stop("'tmRange' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
  }
  if (is.null(minEndIdentity)) minEndIdentity <- 0
  if (minEndIdentity < 0 || minEndIdentity > 1) {
    stop(
      "'minEndIdentity' must be either NULL or from 0 to 1.", call. = FALSE
    )
  }
  if (any(!is.logical(c(gcClamp, avoidThreeEndRuns)))) {
    stop(
      "'gcClamp'  and 'avoidThreeEndRuns' must be set to TRUE or FALSE",
      call. = FALSE
    )
  }
  if (!(min(gcRange) >= 0 && max(gcRange) <= 1)) {
    stop(
      "'gcRange' must be from 0 to 1, e.g. c(0.45, 0.65).",
      call. = FALSE
    )
  }
  if (!(min(tmRange) >= 20 && max(tmRange) <= 90)) {
    stop("'tmRange' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
  }
  oligoLength <- if (probe) c(primerLength, probeLength) else primerLength
  maxDegeneracy <- if (probe) {
    max(maxDegeneracyPrimer, maxDegeneracyProbe)
  } else {
    maxDegeneracyPrimer
  }
  oligos <- .designOligos(x,
                          maxGapFrequency = maxGapFrequency,
                          oligoLength = oligoLength,
                          maxDegeneracy = maxDegeneracy,
                          concOligo = concPrimer,
                          concNa = concNa)
  primers <- .filterPrimers(oligos,
                            length = primerLength,
                            maxDegeneracy = maxDegeneracyPrimer,
                            minEndIdentity = minEndIdentityPrimer,
                            gcClamp = gcClampPrimer,
                            avoidThreeEndG = avoidThreeEndGPrimer,
                            gcRange = gcRangePrimer,
                            tmRange = tmRangePrimer)
  if (probe) {
    probes <- .filterProbes(oligos,
                            length = probeLength,
                            maxDegeneracy = maxDegeneracyProbe,
                            avoidFiveEndG = avoidFiveEndGProbe,
                            gcRange = gcRangeProbe,
                            tmRange = tmRangeProbe,
                            concProbe = concProbe)
    oligos <- rbind(primers, probes)
  } else {
    oligos <- primers
  }
  oligos <- .arrangeData(oligos)
  oligos
}

# Helpers ======================================================================

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
#' \code{.countDegeneracy()} finds the number of variants of
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
#' \code{RprimerProfile} object, and returns a list with start and end position,
#' length, IUPAC sequence (i.e. the DNA sequence with wobble bases), degeneracy
#' (number of variants of each oligo), maximum gap frequency,
#' mean overall identity and minimum 3'-end identity at
#' both forward and reverse direction (the 3'-end is here seen as the last five
#' bases of the oligo).
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param oligoLength An integer [14, 30], defaults to 20.
#'
#' @return A list with all candidate oligos.
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateOligos(exampleRprimerProfile, oligoLength = 18)
#'
#' @keywords internal
#'
#' @noRd
.generateOligos <- function(x, oligoLength = 20) {
  if (!(min(oligoLength) >= 14 && max(oligoLength) <= 30)) {
    stop("'oligoLength' must be from 14 to 30.", call. = FALSE)
  }
  iupacSequence <- .getNmers(x$iupac, oligoLength)
  start <- seq_len(nrow(iupacSequence)) + min(x$position) - 1
  end <- seq_len(nrow(iupacSequence)) + oligoLength - 1 + min(x$position) - 1
  length <- rep(oligoLength, nrow(iupacSequence))
  degeneracy <- apply(iupacSequence, 1, .countDegeneracy)
  gapFrequency <- apply(.getNmers(x$gaps, oligoLength), 1, max)
  identity <- .getNmers(x$identity, oligoLength)
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
#' removes oligos with too high gap frequency and degeneracy.
#'
#' @param x An output from \code{.generateOligos()}.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency. A number [0, 1], defaults to 0.1.
#'
#' @param maxDegeneracy
#' Maximum allowed number of variants of each oligo.
#' An integer [1, 32], defaults to 4.
#'
#' @return
#' A list with the same structure as \code{x}, but where invalid oligos have
#' been removed.
#'
#' @keywords internal
#'
#' @noRd
.filterOligos <- function(x, maxGapFrequency = 0.1, maxDegeneracy = 4) {
  if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
    stop("'maxGapFrequency' must be from 0 to 1.", call. = FALSE)
  }
  if (!(maxDegeneracy >= 1 && maxDegeneracy <= 32)) {
    stop("'maxDegeneracy' must be from 1 to 32.", call. = FALSE)
  }
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
#' \code{makeOligoMatrix()} is part of a "workaround" to avoid loops when
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
#' A GC-clamp is here identified as the presence of two to three G or C:s at
#' the 3'-end (i.e. the last five bases) of an oligo.
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A numeric vector or matrix, where each row corresponds to a specific oligo.
#' In this matrix, 1 corresponds to G or C, and 0 corresponds to A or T.
#'
#' @param fwd
#' If the check should be done in forward direction (i.e. on forward oligos),
#' defaults to \code{TRUE}. If \code{FALSE}, the check will be performed in
#' reverse direction (on reverse oligos).
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
  if (fwd) {
    end <- x[, seq(ncol(x) - 4, ncol(x)), drop = FALSE]
    ifelse(rowSums(end) >= 2 & rowSums(end) <= 3, TRUE, FALSE)
  } else {
    end <- x[, seq_len(5),  drop = FALSE]
    ifelse(5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3, TRUE, FALSE)
  }
}

#' Identify oligos with runs of the same nucleotide at the 3' end
#'
#' \code{.detectThreeEndRuns()} detects if the same nucleotide is repeated at
#' at least 3 times at the terminal 3'-end of an oligo (e.g. "AAA")
#' (bad to have on primers and probes).
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A character vector or matrix, where each row corresponds to a specific oligo.
#'
#' @param fwd
#' If the check should be done in forward direction (i.e. on forward oligos),
#' defaults to \code{TRUE}. If \code{FALSE}, the check will be performed in
#' reverse direction (on reverse oligos).
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
    x[, seq(ncol(x) - 2, ncol(x)), drop = FALSE]
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
#' terminal 3'-end run (bad for primers and probes) and a terminal
#' five end G (good for probes), and melting temperature.
#'
#' Helper function to \code{.designOligos()},
#'
#' @param x An output from \code{.filterOligos()}.
#'
#' @param concOligo
#' Oligo concentration in nM (for tm-calculation). A number
#' [20, 2000] Defaults to 250 nM.
#'
#' @param concNa
#' Sodium ion concentration in the PCR reaction in M (for tm-calculation).
#' A number [0.01, 1]. Defaults to 0.05 M (50 mM).
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
                      oligoLength = 18:22,
                      maxGapFrequency = 0.1,
                      maxDegeneracy = 4,
                      concOligo = 500,
                      concNa = 0.05) {
  allOligos <- lapply(oligoLength, function(i) {
    oligos <- .generateOligos(x, oligoLength = i)
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

#' Exclude
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @param x
#'
#' @return
#'
#' @keywords internal
#'
#' @noRd
.exclude <- function(x) {
  dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})" # sry ;)
  mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
  x <- x[!grepl(dinucleotideRepeats, x)]
  x <- x[!grepl(mononucleotideRepeates, x)]
  x
}

#' Find oligos that pass the criteria for being a primer
#'
#' Helper function to \code{getOligos()}
#'
#' @param x An output from \code{.generateOligos()}
#'
#' @param length
#' Primer length. A numeric vector [14, 30].
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracy
#' Maximum allowed number of variants of each oligo.
#' An integer [1, 32], defaults to 4.
#'
#' @param minEndIdentity
#' A number. The minimum allowed identity
#' at the 3' end of the primer (i.e. the last five bases), defaults to 0.
#'
#' @param gcClamp
#' \code{TRUE} or \code{FALSE}. If primers must have a GC-clamp.
#' Defaults to \code{TRUE}. A GC-clamp
#' is identified as two to three G or
#' C:s within the last five bases (3'-end) of the primer.
#'
#' @param avoidThreeEndRuns
#' \code{TRUE} or \code{FALSE}.
#' If primers with more than two runs
#' of the same nucleotide at the terminal 3'-end should be excluded.
#' Defaults to \code{TRUE}.
#'
#' @param gcRange
#' GC-content range for primers (proportion, not percent).
#' A numeric vector. Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRange
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
                           length = 18:22,
                           maxDegeneracy = 4,
                           minEndIdentity = 0,
                           gcClamp = TRUE,
                           avoidThreeEndRuns = TRUE,
                           gcRange = c(0.45, 0.55),
                           tmRange = c(55, 65),
                           acceptanceThreshold = 0.9) {
  if (!gcClamp) {
    x$gcClampFwdMean <- 1
    x$gcClampRevMean <- 1
  }
  if (!avoidThreeEndRuns) {
    x$threeEndRunsFwdMean <- 1
    x$threeEndRunsRevMean <- 1
  }
  x <- x[x$length >= min(length) & x$length <= max(length), ]
  x <- x[x$degeneracy <= maxDegeneracy, ]
  validFwd <- ifelse(
    x$endIdentityFwd >= minEndIdentity &
      x$gcClampFwdMean >= acceptanceThreshold &
      1 - x$threeEndRunsFwdMean >= acceptanceThreshold, TRUE, FALSE
  )
  validRev <- ifelse(
    x$endIdentityRev >= minEndIdentity &
      x$gcClampRevMean >= acceptanceThreshold &
      1 - x$threeEndRunsRevMean >= acceptanceThreshold, TRUE, FALSE
  )
  x <- cbind(x, validFwd, validRev)
  x <- x[x$validFwd | x$validRev, ]
  okForGc <- .getProportionInRange(x$gcContent, gcRange)
  x <- x[okForGc >= acceptanceThreshold, ]
  okForTm <- .getProportionInRange(x$tm, tmRange)
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
#' @param length
#' Primer length. A numeric vector.
#' Defaults to \code{18:22}.
#'
#' @param maxDegeneracy
#' Maximum allowed number of variants of each oligo.
#' An integer, defaults to 4.
#'
#' @param avoidFiveEndG
#' If probes with a G at the terminal 5'-end should be avoided.
#' Defaults to TRUE.
#'
#' @param gcRange
#' GC-content range for probes (proportion, not percent).
#' A numeric vector. Defaults to \code{c(0.45, 0.55)}.
#'
#' @param tmRange
#' Tm range for probes.
#' A numeric vector. Defaults to \code{c(55, 65)}.
#'
#' @param acceptanceThreshold
#' Proportion of the sequence variants of each oligo that must
#' meet the design constraints. Defaults to 0.9.
#'
#' @keywords internal
#'
#' @noRd
.filterProbes <- function(x,
                          length = 18:22,
                          maxDegeneracy = 4,
                          avoidFiveEndG = TRUE,
                          gcRange = c(0.45, 0.55),
                          tmRange = c(55, 65),
                          concProbe = 250,
                          acceptanceThreshold = 0.9) {

  if (!avoidFiveEndG) {
    x$fiveEndGPlusMean <- 1
    x$fiveEndGMinusMean <- 1
  }
  x <- x[x$length >= min(length) & x$length <= max(length), ]
  x <- x[x$degeneracy <= maxDegeneracy, ]
  validFwd <- ifelse(
      1 - x$fiveEndGPlusMean >= acceptanceThreshold, TRUE, FALSE
  )
  validRev <- ifelse(
    1 - x$fiveEndGMinusMean >= acceptanceThreshold, TRUE, FALSE
  )
  x <- cbind(x, validFwd, validRev)
  x <- x[x$validFwd | x$validRev, ]
  okForGc <- .getProportionInRange(x$gcContent, gcRange)
  x <- x[okForGc >= acceptanceThreshold, ]
  # correct tm for probe conc .....######################################################################
  okForTm <- .getProportionInRange(x$tm, tmRange)
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
    "identity", "degeneracy", "gcContentMean", "tmMean", "alignmentStart",
    "alignmentEnd"
  )
  x <- x[keep]
  x <- x[order(x$start), ]
  rownames(x) <- NULL
  x
}
