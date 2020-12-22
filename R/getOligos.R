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
#' Tm is calculated using the nearest-neighbor method,
#' with the following assumptions:
#' 1) Oligos are not expected to be self-complementary (i.e. no symmetry
#' correction is done);
#' 2) The oligo concentration is assumed to be much higher
#' than the target concentration. See references for table values and equations.
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
    drop <- c("identity3End", "identity3EndRc")
    allOligos <- allOligos[!names(allOligos) %in% drop]
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
    start <- seq_len(length(x) - n + 1)
    end <- start + n - 1
    nmer <- purrr::map_chr(start, function(i) {
        paste(x[start[[i]]:end[[i]]], collapse = "")
    })
    nmer
}

#' Calculate running sums
#'
#' \code{.runningSum()} calculates 'running' sums of a numeric vector. Each
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

#' Calculate mean identity for 3' ends of oligos.
#'
#' @param x A numeric vector of nucleotide identities.
#'
#' @param n Oligo length, an integer.
#'
#' @return
#' Min identity score for the last five bases of the oligo(s).
#' A numeric matrix with two columns, the first represent plus strand oligos
#' and the second represent minus strand oligos (reverse complements).
#'
#' @keywords internal
#'
#' @noRd
.countEndIdentity <- function(x, n) {
    start <- seq_len(length(x) - n + 1)
    end <- start + n - 1
    frame <- purrr::map(start, ~ x[start[[.x]]:end[[.x]]])
    endScore <- purrr::map(frame, function(x) {
        lastFive <- min(x[(length(x) - 4):length(x)])
        firstFive <- min(x[seq_len(5)])
        c(lastFive, firstFive)
    })
    endScore <- do.call("rbind", endScore)
    colnames(endScore) <- c("pos", "neg")
    endScore
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{.countDegeneracy()} returns the number of unique variants of
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
    nNucleotides <- rprimerGlobals$degeneracyLookup[x]
    degeneracy <- prod(nNucleotides)
    degeneracy
}

#' Generate oligos of a specific length
#'
#' @inheritParams getOligos
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @noRd
.generateOligos <- function(x,
                            oligoLength = 20,
                            maxGapFrequency = 0.1,
                            maxDegeneracy = 4) {
    if (!(min(oligoLength) >= 14 && max(oligoLength) <= 30)) {
        stop("'oligoLength' must be from 14 to 30.", call. = FALSE)
    }
    if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
        stop("'maxGapFrequency' must be from 0 to 1.", call. = FALSE)
    }
    if (!(maxDegeneracy >= 1 && maxDegeneracy <= 32)) {
        stop("'maxDegeneracy' must be from 1 to 32.", call. = FALSE)
    }
    majority <- .getNmers(x$majority, n = oligoLength)
    iupac <- .getNmers(x$iupac, n = oligoLength)
    degeneracy <- purrr::map_dbl(iupac, ~ .countDegeneracy(.x))
    start <- seq_along(majority) + min(x$position) - 1
    end <- seq_along(majority) + oligoLength - 1 + min(x$position) - 1
    length <- oligoLength
    identity <- .runningSum(x$identity, n = oligoLength) / oligoLength
    endIdentity <- .countEndIdentity(x$identity, n = oligoLength)
    identity3End <- endIdentity[, "pos"]
    identity3EndRc <- endIdentity[, "neg"]
    gapBin <- ifelse(x$gaps > maxGapFrequency, 1L, 0L)
    gapPenalty <- .runningSum(gapBin, n = oligoLength)
    alignmentStart <- min(x$position)
    alignmentEnd <- max(x$position)
    oligos <- data.frame(
        start, end, length, majority, identity, identity3End, identity3EndRc,
        iupac, degeneracy, gapPenalty, alignmentStart, alignmentEnd
    )
    uniqueOligos <- match(oligos$majority, unique(oligos$majority))
    oligos <- oligos[uniqueOligos, ]
    oligos <- oligos[oligos$gapPenalty == 0, ]
    oligos <- oligos[oligos$degeneracy <= maxDegeneracy, ]
    oligos <- dplyr::select(oligos, -gapPenalty)
    oligos
}

#' Exclude non optimal oligos
#'
#' \code{.exclude()} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}.
#'
#' @param x A data frame with oligos.
#'
#' @details
#' An oligo is excluded if:
#' - It has more than three runs of the same dinucleotide (e.g. "TATATATA")
#' - It has more than four runs of the same nucleotide (e.g. "AAAAA")
#'
#' @return
#' A data frame where non optimal oligos have been excluded.
#'
#' @keywords internal
#'
#' @noRd
.exclude <- function(x) {
    dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})"
    mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
    x <- x[!grepl(dinucleotideRepeats, x$majority), ]
    x <- x[!grepl(mononucleotideRepeates, x$majority), ]
    x
}

#' Get the reverse complement of a DNA sequence
#'
#' \code{.reverseComplement()} finds the reverse complement of a DNA sequence.
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
    complement <- rprimerGlobals$complementLookup[unlist(x)]
    complement <- unname(complement)
    rc <- rev(complement)
    rc <- paste(rc, collapse = "")
    rc
}

#' Add reverse complement to generated oligos
#'
#' @param x A data frame with oligos.
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @noRd
.addReverseComplement <- function(x) {
    majorityRc <- purrr::map_chr(x$majority, ~ .reverseComplement(.x))
    iupacRc <- purrr::map_chr(x$iupac, ~ .reverseComplement(.x))
    x <- tibble::add_column(x, majorityRc, .after = "majority")
    x <- tibble::add_column(x, iupacRc, .after = "iupac")
    x
}

#' Calculate GC content of a DNA sequence
#'
#' \code{.gcContent()} finds the GC content of a DNA sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#'
#' @return The GC content of x.
#'
#' @keywords internal
#'
#' @noRd
.gcContent <- function(x) {
    x <- toupper(x)
    x <- .splitSequence(x)
    gcCount <- length(which(x == "C" | x == "G"))
    totalCount <- length(which(x == "A" | x == "C" | x == "G" | x == "T"))
    gcCount / totalCount
}

#' Add GC content to generated oligos
#'
#' @param x A data frame with oligos.
#'
#' @inheritParams getOligos
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @noRd
.addGcContent <- function(x, gcRange = c(0.45, 0.65)) {
    if (!(min(gcRange) >= 0 && max(gcRange) <= 1)) {
        stop(
            "'gcRange' must be from 0 to 1, e.g. c(0.45, 0.65).",
            call. = FALSE
        )
    }
    gcMajority <- purrr::map_dbl(x$majority, ~ .gcContent(.x))
    x <- tibble::add_column(x, gcMajority, .before = "identity")
    x <- x[x$gcMajority >= min(gcRange) & x$gcMajority <= max(gcRange), ]
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
    ends <- purrr::map(x, ~ .splitSequence(.x))
    ends <- purrr::map(ends, ~ .x[(length(.x) - 4):length(.x)])
    gc <- purrr::map_dbl(ends, ~ .gcContent(.x))
    x[gc > 3 / 5] <- NA
    x[gc < 2 / 5] <- NA
    x
}

#' Replace unwanted oligos with NA
#'
#' @param x A data frame with oligos (with reverse complements).
#'
#' @inheritParams getOligos
#'
#' @return A data frame where unwanted oligos are replaced with NA.
#'
#' @keywords internal
#'
#' @noRd
.filterOligos <- function(x,
                          gcClamp = TRUE,
                          avoid5EndG = FALSE,
                          avoid3EndRuns = TRUE,
                          minEndIdentity = NULL) {
    if (any(!is.logical(c(gcClamp, avoid5EndG, avoid3EndRuns)))) {
        stop(
            "'gcClamp', 'avoid5EndG',  and 'avoid3EndRuns' must be set to
        TRUE or FALSE",
            call. = FALSE
        )
    }
    if (is.null(minEndIdentity)) minEndIdentity <- 0
    if (minEndIdentity < 0 || minEndIdentity > 1) {
        stop(
            "'minEndIdentity' must be either NULL or from 0 to 1.",
            call. = FALSE
        )
    }
    if (gcClamp) {
        x$majority <- .getOligosWithGcClamp(x$majority)
        x$majorityRc <- .getOligosWithGcClamp(x$majorityRc)
    }
    if (avoid5EndG) {
        x$majority[grepl("^G", x$majority)] <- NA
        x$majorityRc[grepl("^G", x$majorityRc)] <- NA
    }
    if (avoid3EndRuns) {
        x$majority[grepl("([A-Z])\\1\\1$", x$majority)] <- NA
        x$majorityRc[grepl("([A-Z])\\1\\1$", x$majorityRc)] <- NA
    }
    x$majority[x$identity3End < minEndIdentity] <- NA
    x$majorityRc[x$identity3EndRc < minEndIdentity] <- NA
    x$iupac[is.na(x$majority)] <- NA
    x$iupacRc[is.na(x$majorityRc)] <- NA
    invalidOligos <- is.na(x$majority) & is.na(x$majorityRc)
    x <- x[!invalidOligos, ]
    x
}

#' Add Tm to oligos
#'
#' @param x A data frame with oligos.
#'
#' @inheritParams getOligos
#'
#' @keywords internal
#'
#' @noRd
.addTm <- function(x, concOligo = 500, concNa = 0.05, tmRange = c(55, 65)) {
    if (!(min(tmRange) >= 20 && max(tmRange) <= 90)) {
        stop("'tmRange' must be from 20 to 90, e.g. c(55, 60).", call. = FALSE)
    }
    if (nrow(x) >= 1) {
        tmMajority <- .tm(x$majority, concOligo = concOligo, concNa = concNa)
    } else {
        tmMajority <- numeric(0)
    }
    x <- tibble::add_column(x, tmMajority, .before = "identity")
    x <- x[x$tmMajority >= min(tmRange) & x$tmMajority <= max(tmRange), ]
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
        allBases <- unname(rprimerGlobals$degenerateLookup[[i]])
        allBases <- unlist(strsplit(allBases, split = ","))
        allBases
    })
    expanded <- expand.grid(
        expanded[seq_along(expanded)],
        stringsAsFactors = FALSE
    )
    expanded <- purrr::map(
        seq_len(nrow(expanded)), ~ paste(expanded[.x, ], collapse = "")
    )
    expanded <- unlist(expanded, use.names = FALSE)
    expanded
}

#' Add info on sequence, GC-content and Tm of all oligo variants
#'
#' @param x A data frame with oligos.
#'
#' @inheritParams getOligos
#'
#' @return A data frame with oligos, including all variants.
#'
#' @keywords internal
#'
#' @noRd
.expandOligos <- function(x, concOligo = 5e-7, concNa = 0.05) {
    all <- purrr::map(x$iupac, function(x) {
        if (is.na(x)) "" else .expandDegenerates(x)
    })
    x <- tibble::add_column(x, all, .before = "alignmentStart")
    allRc <- purrr::map(x$iupacRc, function(x) {
        if (is.na(x)) "" else .expandDegenerates(x)
    })
    x <- tibble::add_column(x, allRc, .after = "all")
    gcAll <- purrr::map(seq_len(nrow(x)), function(i) {
        toCalculate <- ifelse(!is.na(x[i, ]$majority), x[i, ]$all, x[i, ]$allRc)
        toCalculate <- unlist(toCalculate)
        purrr::map_dbl(toCalculate, ~ round(.gcContent(.x), 2))
    })
    tmAll <- purrr::map(seq_len(nrow(x)), function(i) {
        toCalculate <- ifelse(!is.na(x[i, ]$majority), x[i, ]$all, x[i, ]$allRc)
        toCalculate <- unlist(toCalculate)
        purrr::map_dbl(toCalculate, ~ round(.tm(.x), 2))
    })
    x <- tibble::add_column(x, gcAll, tmAll, .before = "alignmentStart")
    x
}

.getPrimers <- function(x,
                        lengthPrimer = 18:22,
                        maxGapFrequencyPrimer = 0.1,
                        maxDegeneracyPrimer = 4,
                        gcClampPrimer = TRUE,
                        avoid3EndRunsPrimer = TRUE,
                        minEndIdentityPrimer = 0.98,
                        gcRangePrimer = c(0.45, 0.55),
                        tmRangePrimer = c(55, 65),
                        concPrimer = 500,
                        concNa = 0.05) {
    allPrimers <- purrr::map_dfr(lengthPrimer, function(i) {
        primers <- .generateOligos(
            x,
            oligoLength = i, maxGapFrequency = maxGapFrequencyPrimer,
            maxDegeneracy = maxDegeneracyPrimer
        )
        primers <- .exclude(primers)
        primers <- .addGcContent(primers, gcRange = gcRangePrimer)
        primers <- .addTm(
            primers,
            concOligo = concPrimer,
            concNa = concNa, tmRange = tmRangePrimer
        )
        primers <- .addReverseComplement(primers)
        primers <- .filterOligos(
            primers,
            gcClamp = gcClampPrimer, avoid5EndG = FALSE,
            avoid3EndRuns = avoid3EndRunsPrimer,
            minEndIdentity = minEndIdentityPrimer
        )
        primers <- tibble::add_column(
            primers,
            type = rep("primer", nrow(primers)), .before = "start"
        )
        primers
    })
    allPrimers <- .expandOligos(
        allPrimers,
        concOligo = concPrimer, concNa = concNa
    )
    allPrimers
}

.getProbes <- function(x,
                       lengthProbe = 18:22,
                       maxGapFrequencyProbe = 0.1,
                       maxDegeneracyProbe = 4,
                       avoid5EndGProbe = TRUE,
                       gcRangeProbe = c(0.45, 0.55),
                       tmRangeProbe = c(55, 70),
                       concProbe = 250,
                       concNa = 0.05) {
    allProbes <- purrr::map_dfr(lengthProbe, function(i) {
        probes <- .generateOligos(
            x,
            oligoLength = i, maxGapFrequency = maxGapFrequencyProbe,
            maxDegeneracy = maxDegeneracyProbe
        )
        probes <- .exclude(probes)
        probes <- .addGcContent(probes, gcRange = gcRangeProbe)
        probes <- .addTm(
            probes,
            concOligo = concProbe,
            concNa = concNa, tmRange = tmRangeProbe
        )
        probes <- .addReverseComplement(probes)
        probes <- .filterOligos(
            probes,
            gcClamp = FALSE, avoid5EndG = avoid5EndGProbe,
            avoid3EndRuns = FALSE,
            minEndIdentity = NULL
        )
        probes <- tibble::add_column(
            probes,
            type = rep("probe", nrow(probes)), .before = "start"
        )
        probes
    })
    allProbes <- .expandOligos(
        allProbes,
        concOligo = concProbe, concNa = concNa
    )
    allProbes
}
