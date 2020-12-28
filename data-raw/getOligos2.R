
#===============================================================================

#' Get oligos
#'
#' \code{getOligos()} identifies oligos (primers and probes)
#' from an \code{RprimerProfile} object.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param lengthPrimer
#' Primer length. A numeric vector [14, 30].
#' Defaults to \code{18:22}.
#'
#' @param maxGapFrequencyPrimer
#' Maximum allowed gap frequency for primers.
#' A number [0, 1]. Defaults to 0.1.
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
#' C:s within the last five bases (3' end) of the oligo.
#'
#' @param avoid3EndRunsPrimer
#' \code{TRUE} or \code{FALSE}.
#' If primers with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}.
#' Defaults to \code{TRUE}.
#'
#' @param minEndIdentityPrimer
#' Optional. If specified, a number [0, 1]. The minimum allowed identity
#' at the 3' end of the primer (i.e. the last five bases).
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
#' @param maxGapFrequencyProbe
#' Maximum allowed gap frequency for probes.
#' A number [0, 1]. Defaults to 0.1.
#'
#' @param maxDegeneracyProbe
#' Maximum number of variants of each probe. A number [1, 32]. Defaults to 4.
#'
#' @param avoid5EndGProbe
#' \code{TRUE} or \code{FALSE}. If probes with G
#' at the 5' end should be excluded. Defaults to \code{TRUE}.
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
#' \item Majority oligos with more than than three consecutive runs of
#' the same dinucleotide (e.g. "TATATATA")
#' \item Majority oligos with more than four consecutive runs of
#' the same nucleotide  (e.g. "AAAAA")
#' \item Majority oligos that are duplicated
#' (to prevent binding at several places on the genome)
#' }
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
#' See references for table values and equations.  INSERT FORMULA
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
#'   \item{majority}{Majority sequence.}
#'   \item{majorityRc}{Majority sequence, reverse complement.}
#'   \item{gcMajority}{GC content (majority sequence), proportion.}
#'   \item{tmMajority}{Melting temperature.}
#'   \item{identity}{Average identity score of the oligo.}
#'   \item{iupac}{IUPAC sequence (i.e. with degenerate bases).}
#'   \item{iupacRc}{IUPAC sequence, reverse complement.}
#'   \item{degeneracy}{Number of variants of the degenerate oligo.}
#'   \item{all}{Lists with all sequence variants of the oligos.}
#'   \item{allRc}{Lists with all sequence variants of the oligos, reverse
#'   complements.}
#'   \item{gcAll}{Lists with the GC content of all
#'   sequence variants of the oligos.}
#'   \item{tmAll}{Lists with the Tm of all sequence variants of the oligos.}
#'   \item{alignmentStart}{Start position of the input consensus profile.}
#'   \item{alingnmentEnd}{End position of the input consensus profile.}
#' }
#'
#' @examples
#' data("exampleRprimerProfile")
#' ## Design primers and probes with default values
#' getOligos(exampleRprimerProfile)
#'
#' @references
#' SantaLucia Jr, J., & Hicks, D. (2004).
#' The thermodynamics of DNA structural motifs.
#' Annu. Rev. Biophys. Biomol. Struct., 33, 415-440.
#'
#' @export
getOligos <- function(x,
                      lengthPrimer = 18:22,
                      maxGapFrequencyPrimer = 0.1,
                      maxDegeneracyPrimer = 4,
                      gcClampPrimer = TRUE,
                      avoid3EndRunsPrimer = TRUE,
                      minEndIdentityPrimer = 0.98,
                      gcRangePrimer = c(0.45, 0.55),
                      tmRangePrimer = c(55, 65),
                      concPrimer = 500,
                      probe = TRUE,
                      lengthProbe = 18:22,
                      maxGapFrequencyProbe = 0.1,
                      maxDegeneracyProbe = 4,
                      avoid5EndGProbe = TRUE,
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
  x <- as.data.frame(x)
  allOligos <- .getPrimers(x,
                           lengthPrimer = lengthPrimer,
                           maxGapFrequencyPrimer = maxGapFrequencyPrimer,
                           maxDegeneracyPrimer = maxDegeneracyPrimer,
                           gcClampPrimer = gcClampPrimer,
                           avoid3EndRunsPrimer = avoid3EndRunsPrimer,
                           minEndIdentityPrimer = minEndIdentityPrimer,
                           gcRangePrimer = gcRangePrimer,
                           tmRangePrimer = tmRangePrimer,
                           concPrimer = concPrimer,
                           concNa = concNa
  )
  if (nrow(allOligos) == 0L) {
    stop("No primers were found.", call. = FALSE)
  }
  if (probe) {
    allProbes <- .getProbes(x,
                            lengthProbe = lengthProbe,
                            maxGapFrequencyProbe = maxGapFrequencyProbe,
                            maxDegeneracyProbe = maxDegeneracyProbe,
                            avoid5EndGProbe = avoid5EndGProbe,
                            gcRangeProbe = gcRangeProbe,
                            tmRangeProbe = tmRangeProbe,
                            concProbe = concProbe,
                            concNa = concNa
    )
    if (nrow(allProbes) == 0L) {
      stop("No probes were found.", call. = FALSE)
    }
    allOligos <- dplyr::bind_rows(allOligos, allProbes)
  }
  drop <- c("identity3End", "identity3EndRc")########################
  allOligos <- allOligos[!names(allOligos) %in% drop] ###################
  allOligos <- allOligos[order(allOligos$start), ]
  RprimerOligo(allOligos)
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
  unlist(strsplit(x, split = ""), use.names = FALSE)
}

#' "Split" a vector into n-sized chunks, and make a matrix
#'
#' \code{.getNmers} divides a vector into chunks of size \code{n}.
#'
#' @param x A vector.
#'
#' @param n Chunk-size (a positive integer).
#'
#' @return A matrix with \code{n} columns.
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
#' \code{.countDegeneracy()} returns the number of unique variants of
#' a DNA sequence with degenerate bases.
#'
#' @param x
#' A DNA sequence (a character vector).
#'
#' @return The number of unique sequences of x (an integer).
#'
#' @keywords internal
#'
#' @noRd
.countDegeneracy <- function(x) {
  x <- toupper(x)
  nNucleotides <- lookup$degeneracy[x]
  prod(nNucleotides)
}

# addOligoLength at the beginning
.generateOligos <- function(x, oligoLength = 20) {
  if (!(min(oligoLength) >= 14 && max(oligoLength) <= 30)) {
    stop("'oligoLength' must be from 14 to 30.", call. = FALSE)
  }
  gapFrequency <- apply(.getNmers(x$gaps, oligoLength), 1, max)
  identity <- .getNmers(x$identity, oligoLength)
  endIdentityFwd <- apply(
    identity[, (ncol(identity) - 4):ncol(identity)], 1, min
  )
  endIdentityRev <- apply(identity[, seq_len(5)], 1, min)
  majority <- .getNmers(x$majority, oligoLength)
  iupac <- .getNmers(x$iupac, oligoLength)
  degeneracy <- apply(iupac, 1, .countDegeneracy)
  start <- seq_len(nrow(majority)) + min(x$position) - 1
  end <- seq_len(nrow(majority)) + oligoLength - 1 + min(x$position) - 1
  alignmentStart <- rep(min(x$position, na.rm = TRUE), nrow(majority))
  alignmentEnd <- rep(max(x$position, na.rm = TRUE), nrow(majority))
  oligos <- list(
    "gapFrequency" = gapFrequency, "identity" = identity,
    "endIdentityFwd" = endIdentityFwd, "endIdentityRev" = endIdentityRev,
    "majority" = majority, "iupac" = iupac, "degeneracy" = degeneracy,
    "start" = start, "end" = end, "alignmentStart" = alignmentStart,
    "alignmentEnd" = alignmentEnd
  )
  oligos
}

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
#' @param x A DNA sequence, i.e. a character vector.
#'
#' @return All sequence variants of \code{x}. A matrix.
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

.makeOligoMatrix <- function(x) {
  degeneracy <- vapply(x, nrow, integer(1))
  id <- lapply(seq_along(degeneracy), function(x) rep(x, degeneracy[x]))
  id <- unlist(id)
  x <- do.call("rbind", x)
  rownames(x) <- id
  x
}

.reverseComplement <- function(x) {
  rev <- x[, rev(seq_len(ncol(x)))]
  rc <- unname(lookup$complement[rev])
  rc <- matrix(rc, ncol = ncol(x), byrow = FALSE)
  rownames(rc) <- rownames(rev)
  rc
}

.detectGcClamp <- function(x, fwd = TRUE) {
  if (fwd) {
    end <- x[, (ncol(x) - 4):ncol(x)]
    ifelse(rowSums(end) >= 2 & rowSums(end) <= 3, TRUE, FALSE)
  } else {
    end <- x[, seq_len(5)]
    ifelse(5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3, TRUE, FALSE)
  }
}

.detectThreeEndRuns <- function(x, fwd = TRUE) {
  end <- if (fwd) x[, (ncol(x) - 3):ncol(x)] else x[, seq_len(4)]
  apply(end, 1, function(x) {
    all(x == "A") | all(x == "C") | all(x == "T") | all(x == "G")
  })
}

.getAllVariants <- function(x, concOligo = 500, concNa = 0.05) {
  all <- list()
  all$sequence <- apply(x$iupac, 1, .expandDegenerates)
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
    all$sequence[, ncol(all$sequence)] == "G", TRUE, FALSE
  )
  all$tm <- .tm(all$sequence, concOligo, concNa)
  all$sequence <- apply(all$sequence, 1, paste, collapse = "")
  all$sequenceRc <- apply(all$sequenceRc, 1, paste, collapse = "")
  lapply(all, function(x) unname(split(unname(x), f = as.integer(names(x)))))
}

.getMeanValues <- function(x) {
  x <- x[-c(1, 2)]
  means <- lapply(x, function(y) vapply(y, mean, double(1)))
  do.call("cbind.data.frame", means)
}

.makeOligoDf <- function(x) {
  x <- within(x, rm(identity, gapFrequency, majority))
  x$iupacRc <- .reverseComplement(x$iupac)
  x$iupac <- apply(x$iupac, 1, paste, collapse = "")
  x$iupacRc <- apply(x$iupacRc, 1, paste, collapse = "")
  x <- do.call("cbind.data.frame", x)
  x[ , c(5:6, 3, 9, 1:2, 4, 7:8)]
}

.makeAllVariantDf <- function(x) {
  x <- x[-c(4, 5, 6, 7, 8, 9)]
  data.frame(do.call("cbind", x))
}

data("exampleRprimerProfile")

# function loop ol len
x <- .generateOligos(exampleRprimerProfile, oligoLength = 14)
x <- .filterOligos(x, maxDegeneracy = 4, minEndIdentity = NULL)
all <- .getAllVariants(x)
means <- .getMeanValues(all)
all <- .makeAllVariantDf(all)
x <- .makeOligoDf(x)
x <- cbind(x, means, all)
# .exclude - on all
# do call rbind


.exclude <- function(x) {
  dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})"
  mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
  x <- x[!grepl(dinucleotideRepeats, x$majority), ]
  x <- x[!grepl(mononucleotideRepeates, x$majority), ]
  x
}

# rearrange cols
# unique() unlist(s) if length(unique(s)) != length(s) ....?
## end

#.filterPrimers <- function(x) {}
#.filterProbes
# e.g. if 80% of variants are ok

# gcRange = c(0.45, 0.65)) {
#  if (!(min(gcRange) >= 0 && max(gcRange) <= 1)) {
#    stop(
#      "'gcRange' must be from 0 to 1, e.g. c(0.45, 0.65).",
#      call. = FALSE
#    )
#  }


#concOligo = 500, concNa = 0.05, tmRange = c(55, 65)) {
#  if (!(min(tmRange) >= 20 && max(tmRange) <= 90)) {
#    stop("'tmRange' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
#  }

#if (any(!is.logical(c(gcClamp, avoid5EndG, avoid3EndRuns)))) {
#  stop(
#    "'gcClamp', 'avoid5EndG',  and 'avoid3EndRuns' must be set to
#        TRUE or FALSE",
#    call. = FALSE
#  )
#}

#if (is.null(minEndIdentity)) minEndIdentity <- 0
#if (minEndIdentity < 0 || minEndIdentity > 1) {
#  stop(
#    "'minEndIdentity' must be either NULL or from 0 to 1.",
#    call. = FALSE
#  )
#}

# then document and write tests!!!!!!!!!!!!!!!!!!!!!!!1

# then mtrx based deltaG

######################################################################
