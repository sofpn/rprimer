#' Design primers and probes
#'
#' \code{oligos()} designs oligos (primers and probes)
#' from an \code{RprimerProfile} object.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency at the primer and probe binding sites in
#' the target alignment.
#' A number [0, 1], defaults to \code{0.01}.
#'
#' @param lengthPrimer
#' Primer length range. A numeric vector [15, 30],
#' defaults to \code{c(18, 22)}.
#'
#' @param maxDegeneracyPrimer
#' Maximum number of variants of each primer. A number [1, 64],
#' defaults to \code{4}.
#'
#' @param gcClampPrimer
#' If primers must have a GC-clamp.
#' A GC-clamp
#' is identified as two to three G or
#' C:s within the last five bases (3' end) of the oligo.
#' \code{TRUE} or \code{FALSE}, defaults to \code{TRUE}.
#'
#' @param avoidThreeEndRunsPrimer
#' If primers with more than two runs
#' of the same nucleotide at the terminal 3' end should be avoided.
#' \code{TRUE} or \code{FALSE}, defaults to \code{TRUE}.
#'
#' @param gcPrimer
#' GC-content range for primers.
#' A numeric vector [0, 1], defaults to \code{c(0.40, 0.65)}.
#'
#' @param tmPrimer
#' Tm range for primers (in Celcius degrees).
#' A numeric vector [30, 90], defaults to \code{c(55, 65)}.
#'
#' @param concPrimer
#' Primer concentration in nM, for Tm calculation. A number
#' [20, 2000], defaults to \code{500}.
#'
#' @param designStrategyPrimer
#' \code{"ambiguous"} or \code{"mixed"}. Defaults to \code{"ambiguous"}
#' (see details below).
#'
#' @param probe
#' If probes should be designed. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param lengthProbe
#' Probe length range. A numeric vector [15, 40],
#' defaults to \code{c(18, 22)}.
#'
#' @param maxDegeneracyProbe
#' Maximum number of variants of each probe. A number [1, 64],
#' defaults to \code{4}.
#'
#' @param avoidFiveEndGProbe
#' If probes with G
#' at the 5' end should be avoided. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param gcProbe
#' GC-content range for probes. A numeric vector [0, 1],
#' defaults to \code{c(0.40, 0.65)}.
#'
#' @param tmProbe
#' Tm range for probes (in Celcius degrees).
#' A numeric vector [30, 90], defaults to \code{c(55, 70)}.
#'
#' @param concProbe
#' Primer concentration in nM, for tm calculation. A numeric vector
#' [20, 2000], defaults to \code{250}.
#'
#' @param concNa
#' The sodium ion concentration in the PCR reaction (in M). For calculation of
#' tm and delta G.
#' A numeric vector [0.01, 1], defaults to \code{0.05} (50 mM).
#'
#' @return
#' An \code{RprimerOligo} object, which contains the following information:
#'
#' \describe{
#'   \item{type}{Whether the oligo is a primer or probe.}
#'   \item{fwd}{\code{TRUE} if the oligo is valid in forward direction,
#'     \code{FALSE} otherwise.}
#'   \item{rev}{\code{TRUE} if the oligo is valid in reverse direction,
#'     \code{FALSE} otherwise.}
#'   \item{start}{Start position of the oligo.}
#'   \item{end}{End positon of the oligo.}
#'   \item{length}{Oligo length.}
#'   \item{iupacSequence}{Oligo sequence in IUPAC format
#'   (i.e. with ambiguous bases).}
#'   \item{iupaSequenceRc}{The reverse complement of the iupacSequence.}
#'   \item{identity}{For ambiguous oligos: Average identity of the oligo.
#'     For mixed oligos: Average identity of the 5' (consensus) part of the
#'     oligo. The value can range from 0 to 1.}
#'   \item{coverage}{For ambiguous oligos: Average coverage of the oligo.
#'     For mixed oligos: Average coverage of the 3' (degenerate) part of the
#'     oligo. The value can range from 0 to 1.}
#'   \item{degeneracy}{Number of sequence variants of the oligo.}
#'   \item{gcContentMean}{Mean GC-content of all sequence variants of the oligo.
#'   }
#'   \item{gcContent}{Range in GC-content of all sequence variants of
#'     the oligo.}
#'   \item{tmMean}{Mean tm of all sequence variants of the oligo
#'   (in Celcius degrees).}
#'   \item{tm}{Range in tm of all sequence variants of the oligo
#'   (in Celcius degrees).}
#'   \item{deltaGMean}{Mean delta G of all sequence variants of the oligo
#'   (in kcal/mol).}
#'   \item{deltaG}{Range in delta G of all sequence variants of the oligo
#'   (in kcal/mol).}
#'   \item{sequence}{All sequence variants of the oligo.}
#'   \item{sequenceRc}{Reverse complements of all sequence variants.}
#'   \item{gcContent}{GC-content of all sequence variants.}
#'   \item{tm}{Tm of all sequence variants (in Celcius degrees).}
#'   \item{deltaG}{Delta G of all sequence variants (in kcal/mol).}
#'   \item{method}{Design method used to generate the oligo: "ambiguous",
#'   "mixedFwd" or "mixedRev".}
#'   \item{score}{Oligo score, the lower the better.}
#'   \item{roiStart}{First position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd}{Last position of the input \code{RprimerProfile} object.}
#' }
#'
#' An error message will return if no oligos are found. If so, a good idea
#' could be to re-run the design process with relaxed constraints.
#'
#' @details
#'
#' \strong{Valid oligos}
#'
#' For an oligo to be considered as valid, all sequence variants must fulfill
#' all the specified design constraints.
#'
#' Furthermore, oligos with at least one sequence variant containing
#' more than four consecutive runs
#' of the same
#' nucleotide (e.g. "AAAAA") and/or more than three consecutive runs
#' of the same di-nucleotide (e.g. "TATATATA") will be excluded
#' from consideration.
#'
#' \strong{Calculation of tm and delta G}
#'
#' Melting temperatures
#' are calculated for perfectly matching
#' DNA duplexes using the
#' nearest-neighbor
#' method (SantaLucia and Hicks, 2004), by using the following equation:
#'
#' \deqn{Tm = (\Delta H ^o \cdot 1000) / (\Delta S ^o + R \cdot \log [\mathrm{oligo}]) - 273.15}
#'
#' where \eqn{\Delta H ^o} is
#' the change in enthalpy (in cal/mol) and \eqn{\Delta S ^o} is the
#' change in entropy (in cal/K/mol) when an
#' oligo and a perfectly matching target sequence goes from random coil to
#' duplex formation.
#' \eqn{K} is the gas constant (1.9872 cal/mol \emph{K}).
#'
#' Delta G is calculated at 37 Celcius degrees, for when an oligo and a
#' perfectly matching target
#' sequence goes from random coil to duplex state, by using  the following
#' equation:
#'
#' \deqn{\Delta G ^o _T = (\Delta H ^o \cdot 1000 - T \cdot \Delta S ^o) / 1000}
#
#' For both tm and delta G, the following salt correction method is used
#' for \eqn{\Delta S^o}, as
#' described in SantaLucia and Hicks (2004):
#'
#' \deqn{\Delta S^o [\mathrm{Na^+}] = \Delta S^o [\mathrm{1 M NaCl}] + 0.368 \cdot N / 2 \cdot \log [\mathrm{Na^+}]}
#'
#' where \eqn{N} is the total number of phosphates in the duplex, and [Na+] is
#' the total
#' concentration of monovalent cations.
#'
#' Nearest neighbor table values for \eqn{\Delta S^o} and \eqn{\Delta H^o} are
#' from SantaLucia and Hicks, 2004, and can be
#' retrieved calling \code{rprimer:::lookup$nn}.
#'
#' \strong{Primer design strategies}
#'
#' Primers can be generated by using one of the two following strategies:
#'
#' \itemize{
#' \item{The \strong{ambiguous strategy} (default) generates primers from the
#' IUPAC consensus sequence, which means that ambiguous bases can
#' occur at any position in the primer.}
#'
#' \item{The \strong{mixed strategy} generates primers from both the majority
#' and the IUPAC consensus sequence. These primers consist of a shorter
#' degenerate part at the 3' end (approx. 1/3 of the primer, targeting a
#' conserved
#' region) and a longer consensus part at the 5' end (approx.
#' 2/3 of the primer),
#' which instead of having ambiguous bases contains the most probable
#' nucleotide at each position.
#' This strategy is based on the Consensus-Degenerate Hybrid
#' Oligonucleotide Primer (CODEHOP) principle (Rose et al., 1998), and aims to
#' to allow amplification of highly variable targets with minimal degeneracy.
#' The idea is that the degenerate 3' end part will bind specifically to
#' the target sequence in the initial PCR cycles, and promote amplification
#' in spite of eventual mismatches at the 5' consensus part
#' (since 5' end mismatches are generally less detrimental than
#' 3' end mismatches). In this way, the generated products
#' will match the 5' ends of all primers perfectly, which allows them
#' to be efficiently amplified in later PCR cycles.
#' To provide a sufficiently high tm in spite of mismatches, it is
#' recommended to design relatively long primers (at least 25 bases) when using
#' this strategy}
#' }
#'
#' Probes are always designed using the ambiguous strategy.
#'
#' \strong{Scoring system for oligos}
#'
#' All valid oligos are scored based on their identity, coverage,
#' degeneracy,
#' average GC content and tm range. The scoring system is presented below.
#'
#' \strong{Identity and coverage}
#'
#' \tabular{lr}{
#' Value range \tab Score \cr
#' \eqn{(0.99, 1]} \tab 0 \cr
#' \eqn{(0.95, 0.99]} \tab 1 \cr
#' \eqn{(0.90, 0.95]} \tab 2 \cr
#' \eqn{\leq 0.90} \tab 3
#' }
#'
#' \strong{Degeneracy}
#'
#' \tabular{lr}{
#' Value \tab Score \cr
#' 1 \tab 0 \cr
#' 2 or 3 \tab 1 \cr
#' 3 or 4 \tab 2 \cr
#' More than 4 \tab 3
#' }
#'
#' \strong{Average GC-content}
#'
#' This score is based on how much
#' the average GC-content deviates from 0.5 (in absolute value).
#'
#' \tabular{lr}{
#' Value range \tab Score \cr
#' \eqn{[0, 0.05)}
#' \tab 0 \cr
#' \eqn{[0.05, 0.1)}
#' \tab 1 \cr
#' \eqn{[0.1,  0.2)}
#' \tab 2 \cr
#' \eqn{\geq 0.2} \tab 3
#' }
#'
#' \strong{Tm range}
#'
#' \tabular{lr}{
#' Value range \tab Score \cr
#' \eqn{[0, 1)}
#' \tab 0 \cr
#' \eqn{[1, 2)}
#' \tab 1 \cr
#' \eqn{[2, 3)}
#' \tab 2 \cr
#' \eqn{\geq 3} \tab 3
#' }
#'
#' These scores are summarized to a total score, and
#' the weight of each individual score is 1. Thus, the lowest and best
#' possible score for an oligo is 0, and the worst possible score is 12.
#'
#' @references
#' Rose, TM., Schultz ER., Henikoff JG., Pietrokovski S.,
#' McCallum CM., and Henikoff S. 1998. Consensus-Degenerate
#' Hybrid Oligonucleotide Primers for
#' Amplification of Distantly Related Sequences. Nucleic Acids Research 26 (7):
#' 1628-35.
#'
#' SantaLucia Jr, J., & Hicks, D. (2004).
#' The thermodynamics of DNA structural motifs.
#' Annu. Rev. Biophys. Biomol. Struct., 33, 415-440.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- exampleRprimerProfile
#'
#' ## Design primers and probes with default values
#' oligos(x)
#'
#' ## Design primers and probes only within a specific region of interest
#' roi <- x[x$position >= 5000 & x$position <= 6000, ]
#' oligos(roi)
#'
#' ## Design primers only
#' oligos(roi, probe = FALSE)
#'
#' ## Allow higher degeneracy
#' oligos(roi, maxDegeneracyPrimer = 32, probe = FALSE)
#'
#' ## Use the mixed strategy for primers
#' oligos(roi, designStrategyPrimer = "mixed", probe = FALSE)
oligos <- function(x,
                   maxGapFrequency = 0.01,
                   lengthPrimer = c(18, 22),
                   maxDegeneracyPrimer = 4,
                   gcClampPrimer = TRUE,
                   avoidThreeEndRunsPrimer = TRUE,
                   gcPrimer = c(0.40, 0.65),
                   tmPrimer = c(50, 65),
                   concPrimer = 500,
                   designStrategyPrimer = "ambiguous",
                   probe = TRUE,
                   lengthProbe = c(18, 22),
                   maxDegeneracyProbe = 4,
                   avoidFiveEndGProbe = TRUE,
                   gcProbe = c(0.40, 0.65),
                   tmProbe = c(50, 70),
                   concProbe = 250,
                   concNa = 0.05) {
    if (!methods::is(x, "RprimerProfile")) {
        stop("'x' must be an RprimerProfile object.", call. = FALSE)
    }
    if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
        stop("'maxGapFrequency' must be from 0 to 1.", call. = FALSE)
    }
    if (!(min(lengthPrimer) >= 15 && max(lengthPrimer) <= 40)) {
        stop("'lengthPrimer' must be from 15 to 40.", call. = FALSE)
    }
    if (!(maxDegeneracyPrimer >= 1 && maxDegeneracyPrimer <= 64)) {
        stop("'maxDegeneracyPrimer' must be from 1 to 64.", call. = FALSE)
    }
    if (!is.logical(gcClampPrimer)) {
        stop("'gcClampPrimer' must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.logical(avoidThreeEndRunsPrimer)) {
        stop("'avoidThreeEndRunsPrimer' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(min(gcPrimer) >= 0 && max(gcPrimer) <= 1)) {
        stop(
            "'gcPrimer' must be from 0 to 1, e.g. c(0.45, 0.65).",
            call. = FALSE
        )
    }
    if (!(min(tmPrimer) >= 20 && max(tmPrimer) <= 90)) {
        stop(
            "'tmPrimer' must be from 20 to 90, e.g. c(55, 60).",
            call. = FALSE
        )
    }
    if (!(concPrimer >= 20 && concPrimer <= 2000)) {
        stop("'concPrimer' must be from 20 to 2000.", call. = FALSE)
    }
    if (!(
        designStrategyPrimer == "ambiguous" || designStrategyPrimer == "mixed")
    ) {
        stop(
            "'designStrategyPrimer' must be either 'ambiguous' or 'mixed'.",
            call. = FALSE
        )
    }
    if (!is.logical(probe)) {
        stop("'probe' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(min(lengthProbe) >= 15 && max(lengthProbe) <= 40)) {
        stop("'lengthProbe' must be from 15 to 40.", call. = FALSE)
    }
    if (!(maxDegeneracyProbe >= 1 && maxDegeneracyProbe <= 64)) {
        stop("'maxDegeneracyProbe' must be from 1 to 64.", call. = FALSE)
    }
    if (!is.logical(avoidFiveEndGProbe)) {
        stop("'avoidFiveEndGProbe' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(min(gcProbe) >= 0 && max(gcProbe) <= 1)) {
        stop(
            "'gcProbe' must be from 0 to 1, e.g. c(0.45, 0.65).",
            call. = FALSE
        )
    }
    if (!(min(tmProbe) >= 20 && max(tmProbe) <= 90)) {
        stop(
            "'tmProbe' must be from 20 to 90, e.g. c(55, 60).",
            call. = FALSE
        )
    }
    if (!(concProbe >= 20 && concProbe <= 2000)) {
        stop("'concProbe' must be from 20 to 2000.", call. = FALSE)
    }
    if (!(concNa >= 0.01 && concNa <= 1)) {
        stop("'concNa' must be from 0.01 to 1.", call. = FALSE)
    }
    lengthPrimer <- seq(min(lengthPrimer), max(lengthPrimer))
    lengthProbe <- seq(min(lengthProbe), max(lengthProbe))
    if (designStrategyPrimer == "mixed") {
        primers <- .designMixedPrimers(
            x,
            maxGapFrequency = maxGapFrequency,
            lengthPrimer = lengthPrimer,
            maxDegeneracyPrimer = maxDegeneracyPrimer,
            gcClampPrimer = gcClampPrimer,
            avoidThreeEndRunsPrimer = avoidThreeEndRunsPrimer,
            gcPrimer = gcPrimer,
            tmPrimer = tmPrimer,
            concPrimer = concPrimer,
            concNa = concNa
        )
        if (nrow(primers) == 0L) {
            stop("No primers were found.", call. = FALSE)
        }
        oligos <- primers
        if (probe) {
            probes <- .designAmbiguousOligos(
                x,
                primer = FALSE,
                lengthProbe = lengthProbe,
                maxDegeneracyProbe = maxDegeneracyProbe,
                avoidFiveEndGProbe = avoidFiveEndGProbe,
                gcProbe = gcProbe,
                tmProbe = tmProbe,
                concProbe = concProbe,
                concNa = concNa
            )
            if (nrow(probes) == 0L) {
                stop("No probes were found.", call. = FALSE)
            }
            oligos <- rbind(oligos, probes)
        }
    } else {
        oligos <- .designAmbiguousOligos(
            x,
            maxGapFrequency = maxGapFrequency,
            primer = TRUE,
            lengthPrimer = lengthPrimer,
            maxDegeneracyPrimer = maxDegeneracyPrimer,
            gcClampPrimer = gcClampPrimer,
            avoidThreeEndRunsPrimer = avoidThreeEndRunsPrimer,
            gcPrimer = gcPrimer,
            tmPrimer = tmPrimer,
            concPrimer = concPrimer,
            probe = probe,
            lengthProbe = lengthProbe,
            maxDegeneracyProbe = maxDegeneracyProbe,
            avoidFiveEndGProbe = avoidFiveEndGProbe,
            gcProbe = gcProbe,
            tmProbe = tmProbe,
            concProbe = concProbe,
            concNa = concNa
        )
        if (nrow(oligos[oligos$type == "primer", ]) == 0L) {
            stop("No primers were found.", call. = FALSE)
        }
        if (probe) {
            if (nrow(oligos[oligos$type == "probe", ]) == 0L) {
                stop("No probes were found.", call. = FALSE)
            }
        }
    }
    oligos <- .scoreOligos(oligos)
    oligos <- .beautifyOligos(oligos)
    RprimerOligo(oligos)
}

# Helpers ======================================================================

#' @noRd
#'
#' @examples
#' .nmers(c("A", "G", "T", "T", "C", "G"), n = 4)
.nmers <- function(x, n) {
    start <- seq_len(length(x) - n + 1)
    end <- start + n - 1
    nmers <- lapply(start, \(i) x[start[[i]]:end[[i]]])
    do.call("rbind", nmers)
}

#' @noRd
#'
#' @examples
#' .countDegeneracy(c("A", "R", "T", "T", "N", "G"))
.countDegeneracy <- function(x) prod(lookup$degeneracy[x])

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateAmbiguousOligos(exampleRprimerProfile, lengthOligo = 18)
.generateAmbiguousOligos <- function(x, lengthOligo = 20) {
    oligos <- list()
    oligos$iupacSequence <- .nmers(x$iupac, lengthOligo)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    oligos$coverage <- .nmers(x$coverage, lengthOligo) |> rowMeans()
    oligos$identity <- .nmers(x$identity, lengthOligo) |> rowMeans()
    oligos$method <- rep("ambiguous", nrow(oligos$iupacSequence))
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}

#' @noRd
#'
#' @examples
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))))
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))), combine = FALSE)
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))), rev = TRUE)
.splitAndPaste <- function(first, second, rev = FALSE, combine = TRUE) {
    n <- ncol(first)
    small <- seq_len(as.integer(n / 3))
    large <- seq(small[length(small)] + 1, n)
    if (rev) {
        first <- first[, small, drop = FALSE]
        second <- second[, large, drop = FALSE]
    } else {
        first <- first[, seq_along(large), drop = FALSE]
        second <- second[, small + n - length(small), drop = FALSE]
    }
    if (combine) {
        cbind(first, second)
    } else {
        list(first, second)
    }
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .mixFwd(exampleRprimerProfile, 20)
.mixFwd <- function(x, lengthOligo = 20) {
    oligos <- list()
    majority <- .nmers(x$majority, lengthOligo)
    iupac <- .nmers(x$iupac, lengthOligo)
    oligos$iupacSequence <- .splitAndPaste(majority, iupac)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    identityCoverage <- .splitAndPaste(
        .nmers(x$identity, lengthOligo), .nmers(x$coverage, lengthOligo),
        combine = FALSE
    )
    oligos$identity <- identityCoverage[[1]] |> rowMeans()
    oligos$coverage <- identityCoverage[[2]] |> rowMeans()
    oligos$method <- rep("mixedFwd", nrow(oligos$iupacSequence))
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .mixRev(exampleRprimerProfile, 20)
.mixRev <- function(x, lengthOligo = 20) {
    oligos <- list()
    majority <- .nmers(x$majority, lengthOligo)
    iupac <- .nmers(x$iupac, lengthOligo)
    oligos$iupacSequence <- .splitAndPaste(iupac, majority, rev = TRUE)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    coverageIdentity <- .splitAndPaste(
        .nmers(x$coverage, lengthOligo), .nmers(x$identity, lengthOligo),
        combine = FALSE, rev = TRUE
    )
    oligos$identity <- coverageIdentity[[2]] |> rowMeans()
    oligos$coverage <- coverageIdentity[[1]] |> rowMeans()
    oligos$method <- rep("mixedRev", nrow(oligos$iupacSequence))
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}


#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' fwd <- .mixFwd(exampleRprimerProfile)
#' rev <- .mixRev(exampleRprimerProfile)
#' .mergeLists(fwd, rev)
.mergeLists <- function(first, second) {
    x <- lapply(names(first), \(i) {
        if (is.matrix(first[[i]])) {
            rbind(first[[i]], second[[i]])
        } else {
            c(first[[i]], second[[i]])
        }
    })
    names(x) <- names(first)
    x
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateMixedOligos(exampleRprimerProfile, lengthOligo = 20)
.generateMixedOligos <- function(x, lengthOligo = 20) {
    fwd <- .mixFwd(x, lengthOligo)
    rev <- .mixRev(x, lengthOligo)
    .mergeLists(fwd, rev)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateMixedOligos(exampleRprimerProfile)
#' .filterOligos(x)
.filterOligos <- function(x, maxGapFrequency = 0.1, maxDegeneracy = 4) {
    invalidCharacters <- apply(x$iupacSequence, 1, \(x) {
        any(x == "-") | any(is.na(x))
    })
    invalid <- unique(c(
        which(x$degeneracy > maxDegeneracy),
        which(x$gapFrequency > maxGapFrequency),
        which(invalidCharacters)
    ))
    if (length(invalid) > 0L) {
        x <- lapply(x, \(x) {
            if (is.matrix(x)) x[-invalid, , drop = FALSE] else x[-invalid]
        })
    }
    x
}

#' @noRd
#'
#' @examples
#' .expandDegenerates(c("A", "R", "T", "T", "N", "G"))
.expandDegenerates <- function(x) {
    bases <- lapply(x, \(i) {
        degen <- unname(lookup$degenerates[[i]])
        unlist(strsplit(degen, split = ","))
    })
    all <- expand.grid(bases[seq_along(bases)], stringsAsFactors = FALSE)
    all <- as.matrix(all)
    colnames(all) <- NULL
    all
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterOligos(.generateMixedOligos(exampleRprimerProfile))
#' x$sequence <- apply(x$iupacSequence, 1, .expandDegenerates)
#' .oligoMatrix(x$sequence)
.oligoMatrix <- function(x) {
    degeneracy <- vapply(x, nrow, integer(1L))
    id <- lapply(seq_along(degeneracy), \(x) rep(x, degeneracy[[x]]))
    id <- unlist(id)
    x <- do.call("rbind", x)
    rownames(x) <- id
    x
}

#' @noRd
#'
#' @examples
#' .reverseComplement(matrix(c("A", "R", "T", "T", "N", "G")))
.reverseComplement <- function(x) {
    rc <- x[, rev(seq_len(ncol(x))), drop = FALSE]
    rc[] <- lookup$complement[rc]
    rc
}

#' @noRd
#'
#' @examples
#' seq <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0))
#' .gcClamp(gc, rev = FALSE)
#' .gcClamp(gc, rev = TRUE)
.gcClamp <- function(x, rev = FALSE) {
    if (rev) {
        end <- x[, seq_len(5), drop = FALSE]
        5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3 ## Because of complement..
    } else {
        end <- x[, seq((ncol(x) - 4), ncol(x)), drop = FALSE]
        rowSums(end) >= 2 & rowSums(end) <= 3
    }
}

#' @noRd
#'
#' @examples
#' seq <- matrix(c("A", "C", "G", "G", "T", "T", "A", "A"))
#' .endRuns(seq, rev = FALSE)
.endRuns <- function(x, rev = FALSE) {
    if (rev) {
        end <- x[, seq_len(3), drop = FALSE]
    } else {
        end <- x[, seq((ncol(x) - 2), ncol(x)), drop = FALSE]
    }
    apply(end, 1, \(x) {
        all(x == "A") | all(x == "C") | all(x == "T") | all(x == "G")
    })
}

#' @noRd
#'
#' @examples
#' .repeats(c("ACTTTTCT", "ACTTTTTCT", "ATCTCTCTCA"))
.repeats <- function(x) {
    di <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(CT){4,}|(TC){4,}|)"
    mono <- "([A-Z])\\1\\1\\1\\1"
    vapply(x, \(y) {
        grepl(di, y) | grepl(mono, y)
    }, logical(1L))
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterOligos(.generateMixedOligos(exampleRprimerProfile))
#' .allVariants(x)
.allVariants <- function(x,
                         concPrimer = 500,
                         concProbe = 250,
                         concNa = 0.05) {
    all <- list()
    all$sequence <- apply(x$iupacSequence, 1, .expandDegenerates, simplify = FALSE)
    ## If there is only one variant of each oligo,
    ## (and apply returns a matrix instead of a list):
    if (is.matrix(all$sequence)) {
        all$sequence <- t(all$sequence)
        all$sequence <- lapply(seq_len(nrow(all$sequence)), \(i) {
            all$sequence[i, , drop = FALSE]
        })
    }
    all$sequence <- .oligoMatrix(all$sequence)
    all$sequenceRc <- .reverseComplement(all$sequence)
    gc <- all$sequence == "C" | all$sequence == "G"
    n <- rowSums(
        all$sequence == "A" | all$sequence == "C" |
            all$sequence == "G" | all$sequence == "T"
    )
    all$gcContent <- rowSums(gc) / n
    all$gcClampFwd <- .gcClamp(gc)
    all$gcClampRev <- .gcClamp(gc, rev = TRUE)
    all$threeEndRunsFwd <- .endRuns(all$sequence)
    all$threeEndRunsRev <- .endRuns(all$sequence, rev = TRUE)
    all$fiveEndGPlus <- all$sequence[, 1] == "G"
    all$fiveEndGMinus <- all$sequence[, ncol(all$sequence)] == "C"
    tmParam <- .tmParameters(all$sequence, concNa)
    all$tmPrimer <- .tm(tmParam, concPrimer)
    all$tmProbe <- .tm(tmParam, concProbe)
    all$deltaG <- .deltaG(tmParam)
    all$sequence <- apply(all$sequence, 1, paste, collapse = "")
    all$sequenceRc <- apply(all$sequenceRc, 1, paste, collapse = "")
    all$repeats <- .repeats(all$sequence)
    lapply(all, \(x) unname(split(unname(x), f = as.integer(names(x)))))
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .allVariants(.filterOligos(.generateMixedOligos(exampleRprimerProfile)))
#' .meanRange(x)
.meanRange <- function(x) {
    x <- x[c("gcContent", "tmPrimer", "tmProbe", "deltaG")]
    means <- lapply(x, \(y) {
        vapply(y, \(z) {
            sum(z) / length(z)
        }, double(1L))
    })
    means <- do.call("cbind.data.frame", means)
    names(means) <- paste0(names(means), "Mean")
    ranges <- lapply(x, \(y) {
        vapply(y, \(z) {
            max(z) - min(z)
        }, double(1L))
    })
    ranges <- do.call("cbind.data.frame", ranges)
    names(ranges) <- paste0(names(ranges), "Range")
    cbind(means, ranges)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterOligos(.generateMixedOligos(exampleRprimerProfile))
#' .makeOligoDf(x)
.makeOligoDf <- function(x) {
    x <- within(x, rm("gapFrequency"))
    x$iupacSequenceRc <- .reverseComplement(x$iupacSequence)
    x$iupacSequence <- apply(x$iupacSequence, 1, paste, collapse = "")
    x$iupacSequenceRc <- apply(x$iupacSequenceRc, 1, paste, collapse = "")
    do.call("cbind.data.frame", x)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .allVariants(.filterOligos(.generateMixedOligos(exampleRprimerProfile)))
#' .isWithinRange(x$gcContent, c(0.4, 0.6))
.isWithinRange <- function(x, range) {
    lapply(x, \(y) y >= min(range) & y <= max(range))
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .allVariants(.filterOligos(.generateMixedOligos(exampleRprimerProfile)))
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' x <- data.frame(cbind(gcInRange))
#' .convertToMatrices(x["gcInRange"])
.convertToMatrices <- function(x) {
    lapply(seq_len(nrow(x)), \(i) {
        y <- lapply(x[i, , drop = FALSE], unlist)
        do.call("cbind", y)
    })
}

#' Column means represent the proportion
#' of sequence variants that fulfill a specific criteria (e.g. GC-clamp),
#' and row means represent the proportion of the desired design criteria that
#' are fulfilled by specific sequence variants.
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .allVariants(.filterOligos(.generateMixedOligos(exampleRprimerProfile)))
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' x <- data.frame(cbind(gcInRange))
#' check <- .convertToMatrices(x["gcInRange"])
#' .isValid(check, rowThreshold = 0.5, colThreshold = 0.5)
.isValid <- function(x, rowThreshold = 1, colThreshold = 1) {
    valid <- vapply(x, \(y) {
        toInvert <- c(
            "repeats", "threeEndRunsFwd", "threeEndRunsRev",
            "fiveEndGPlus", "fiveEndGMinus"
        )
        select <- colnames(y) %in% toInvert
        y[, select] <- !y[, select]
        col <- colMeans(y)
        row <- rowMeans(y)
        all(col >= colThreshold) & all(row >= rowThreshold)
    }, logical(1L))
    valid
}

#' @noRd
.checkPrimers <- function(x,
                          gcClampPrimer = TRUE,
                          avoidThreeEndRunsPrimer = TRUE) {
    selectFwd <- c("repeats", "tmInRange", "gcInRange")
    selectRev <- c("repeats", "tmInRange", "gcInRange")
    if (gcClampPrimer) {
        selectFwd <- c(selectFwd, "gcClampFwd")
        selectRev <- c(selectRev, "gcClampRev")
    }
    if (avoidThreeEndRunsPrimer) {
        selectFwd <- c(selectFwd, "threeEndRunsFwd")
        selectRev <- c(selectRev, "threeEndRunsRev")
    }
    fwdPrimers <- x[selectFwd]
    fwdPrimers <- .convertToMatrices(fwdPrimers)
    fwd <- .isValid(fwdPrimers)
    revPrimers <- x[selectRev]
    revPrimers <- .convertToMatrices(revPrimers)
    rev <- .isValid(revPrimers)
    x <- cbind(x, fwd, rev)
    x[x$fwd | x$rev, , drop = FALSE]
}

#' @noRd
.filterPrimers <- function(x,
                           lengthPrimer = 18:22,
                           maxDegeneracyPrimer = 4,
                           gcClampPrimer = TRUE,
                           avoidThreeEndRunsPrimer = TRUE,
                           gcPrimer = c(0.45, 0.55),
                           tmPrimer = c(55, 65)) {
    x <- x[x$length %in% lengthPrimer, , drop = FALSE]
    x <- x[x$degeneracy <= maxDegeneracyPrimer, , drop = FALSE]
    gcInRange <- .isWithinRange(x$gcContent, gcPrimer)
    tmInRange <- .isWithinRange(x$tmPrimer, tmPrimer)
    x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
    x <- .checkPrimers(
        x,
        gcClampPrimer,
        avoidThreeEndRunsPrimer
    )
    x$rev[x$method == "mixedFwd" & x$rev] <- FALSE
    x$fwd[x$method == "mixedRev" & x$fwd] <- FALSE
    x <- x[x$fwd | x$rev, , drop = FALSE]
    remove <- c(
        "gcInRange", "tmInRange",
        "okFwd", "okRev", "tmProbeMean", "tmProbeRange", "tmProbe"
    )
    x <- x[!names(x) %in% remove]
    oldnames <- c("tmPrimerMean", "tmPrimerRange", "tmPrimer")
    newnames <- c("tmMean", "tmRange", "tm")
    names(x)[names(x) %in% oldnames] <- newnames
    type <- rep("primer", nrow(x))
    cbind(type, x)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .designMixedPrimers(exampleRprimerProfile)
.designMixedPrimers <- function(x,
                                maxGapFrequency = 0.01,
                                lengthPrimer = 18:22,
                                maxDegeneracyPrimer = 4,
                                gcClampPrimer = TRUE,
                                avoidThreeEndRunsPrimer = TRUE,
                                gcPrimer = c(0.45, 0.55),
                                tmPrimer = c(55, 65),
                                concPrimer = 500,
                                concNa = 0.05) {
    lengthPrimer <- lengthPrimer[order(lengthPrimer)]
    all <- lapply(lengthPrimer, \(i) {
        mixed <- .generateMixedOligos(x, lengthOligo = i)
        mixed <- .filterOligos(
            mixed, maxGapFrequency, maxDegeneracyPrimer
        )
        if (length(mixed[[1]] > 0L)) {
            allVariants <- .allVariants(
                mixed, concPrimer,
                concProbe = NA, concNa
            )
            meansAndRanges <- .meanRange(allVariants)
            allVariants <- data.frame(do.call("cbind", allVariants))
            mixed <- .makeOligoDf(mixed)
            cbind(mixed, meansAndRanges, allVariants)
        } else {
            NULL
        }
    })
    all <- do.call("rbind", all)
    all <- .filterPrimers(
        all,
        lengthPrimer,
        maxDegeneracyPrimer,
        gcClampPrimer,
        avoidThreeEndRunsPrimer,
        gcPrimer,
        tmPrimer
    )
    all
}

#' @noRd
.checkProbes <- function(x,
                         avoidFiveEndGProbe) {
    selectFwd <- c("repeats", "tmInRange", "gcInRange")
    selectRev <- c("repeats", "tmInRange", "gcInRange")
    if (avoidFiveEndGProbe) {
        selectFwd <- c(selectFwd, "fiveEndGPlus")
        selectRev <- c(selectRev, "fiveEndGMinus")
    }
    fwdProbes <- x[selectFwd]
    fwdProbes <- .convertToMatrices(fwdProbes)
    fwd <- .isValid(fwdProbes)
    revProbes <- x[selectRev]
    revProbes <- .convertToMatrices(revProbes)
    rev <- .isValid(revProbes)
    x <- cbind(x, fwd, rev)
    x[x$fwd | x$rev, , drop = FALSE]
}

#' @noRd
.filterProbes <- function(x,
                          lengthProbe = 18:22,
                          maxDegeneracyProbe = 4,
                          avoidFiveEndGProbe = TRUE,
                          gcProbe = c(0.45, 0.55),
                          tmProbe = c(55, 65)) {
    x <- x[x$length %in% lengthProbe, , drop = FALSE]
    x <- x[x$degeneracy <= maxDegeneracyProbe, , drop = FALSE]
    gcInRange <- .isWithinRange(x$gcContent, gcProbe)
    tmInRange <- .isWithinRange(x$tmProbe, tmProbe)
    x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
    x <- .checkProbes(
        x,
        avoidFiveEndGProbe
    )
    remove <- c(
        "gcInRange", "tmInRange", "tmPrimerMean",
        "tmPrimerRange", "tmPrimer"
    )
    x <- x[!names(x) %in% remove]
    oldnames <- c("tmProbeMean", "tmProbeRange", "tmProbe")
    newnames <- c("tmMean", "tmRange", "tm")
    names(x)[names(x) %in% oldnames] <- newnames
    type <- rep("probe", nrow(x))
    cbind(type, x)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .designAmbiguousOligos(exampleRprimerProfile)
.designAmbiguousOligos <- function(x,
                                   maxGapFrequency = 0.01,
                                   primer = TRUE,
                                   lengthPrimer = c(18, 22),
                                   maxDegeneracyPrimer = 4,
                                   gcClampPrimer = TRUE,
                                   avoidThreeEndRunsPrimer = TRUE,
                                   gcPrimer = c(0.40, 0.65),
                                   tmPrimer = c(50, 65),
                                   concPrimer = 500,
                                   probe = TRUE,
                                   lengthProbe = c(18, 22),
                                   maxDegeneracyProbe = 4,
                                   avoidFiveEndGProbe = TRUE,
                                   gcProbe = c(0.40, 0.65),
                                   tmProbe = c(50, 70),
                                   concProbe = 250,
                                   concNa = 0.05) {
    if (probe && primer) {
        lengthOligo <- unique(c(lengthPrimer, lengthProbe))
        maxDegeneracy <- max(c(maxDegeneracyPrimer, maxDegeneracyProbe))
    } else if (!probe && primer) {
        lengthOligo <- lengthPrimer
        maxDegeneracy <- maxDegeneracyPrimer
    } else if (probe && !primer) {
        lengthOligo <- lengthProbe
        maxDegeneracy <- maxDegeneracyProbe
    }
    lengthOligo <- lengthOligo[order(lengthOligo)]
    ambiguous <- lapply(lengthOligo, \(i) {
        amb <- .generateAmbiguousOligos(x, lengthOligo = i)
        amb <- .filterOligos(
            amb,
            maxGapFrequency = maxGapFrequency,
            maxDegeneracy = maxDegeneracy
        )
        if (length(amb[[1]] > 0L)) {
            allVariants <- .allVariants(
                amb,
                concPrimer = concPrimer,
                concProbe = concProbe,
                concNa = concNa
            )
            meanAndRange <- .meanRange(allVariants)
            allVariants <- data.frame(do.call("cbind", allVariants))
            amb <- .makeOligoDf(amb)
            cbind(amb, meanAndRange, allVariants)
        } else {
            NULL
        }
    })
    ambiguous <- do.call("rbind", ambiguous)
    if (primer) {
        primers <- .filterPrimers(
            ambiguous,
            lengthPrimer,
            maxDegeneracyPrimer,
            gcClampPrimer,
            avoidThreeEndRunsPrimer,
            gcPrimer,
            tmPrimer
        )
    } else {
        primers <- NULL
    }
    if (probe) {
        probes <- .filterProbes(
            ambiguous,
            lengthProbe,
            maxDegeneracyProbe,
            avoidFiveEndGProbe,
            gcProbe,
            tmProbe
        )
    } else {
        probes <- NULL
    }
    rbind(primers, probes)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- head(exampleRprimerOligo$identity)
#' .scoreIdentityCoverage(x)
.scoreIdentityCoverage <- function(x) {
    score <- vector(mode = "double", length = length(x))
    score[x <= 1 & x > 0.99] <- 0
    score[x <= 0.99 & x > 0.95] <- 1
    score[x <= 0.95 & x > 0.90] <- 2
    score[x <= 0.90] <- 3
    score
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- head(exampleRprimerOligo$degeneracy)
#' .scoreDegeneracy(x)
.scoreDegeneracy <- function(x) {
    score <- vector(mode = "double", length = length(x))
    score[x == 1] <- 0
    score[x == 2 | x == 3] <- 1
    score[x == 3 | x == 4] <- 2
    score[x > 4] <- 3
    score
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- head(exampleRprimerOligo$gcContentMean)
#' .scoreGcContent(x)
.scoreGcContent <- function(x) {
    deviation <- abs(x - 0.5)
    score <- vector(mode = "double", length = length(deviation))
    score[deviation >= 0 & deviation < 0.05] <- 0
    score[deviation >= 0.05 & deviation < 0.1] <- 1
    score[deviation >= 0.1 & deviation < 0.2] <- 2
    score[deviation >= 0.2] <- 3
    score
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- head(exampleRprimerOligo$tmRange)
#' .scoreTmRange(x)
.scoreTmRange <- function(x) {
    score <- vector(mode = "double", length = length(x))
    score[x >= 0 & x < 1] <- 0
    score[x >= 1 & x < 2] <- 1
    score[x >= 2 & x < 3] <- 2
    score[x >= 3] <- 3
    score
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- head(exampleRprimerOligo)
#' .scoreOligos(x)
.scoreOligos <- function(x) {
    score <- list()
    score$identity <- .scoreIdentityCoverage(x$identity)
    score$coverage <- .scoreIdentityCoverage(x$coverage)
    score$degeneracy <- .scoreDegeneracy(x$degeneracy)
    score$gcContent <- .scoreGcContent(x$gcContentMean)
    score$tm <- .scoreTmRange(x$tmRange)
    score <- do.call("cbind", score)
    score <- rowSums(score)
    cbind(x, score)
}

#' @noRd
.beautifyOligos <- function(x) {
    keep <- c(
        "type", "fwd", "rev", "start", "end", "length",
        "iupacSequence", "iupacSequenceRc", "identity",
        "coverage", "degeneracy", "gcContentMean", "gcContentRange",
        "tmMean", "tmRange", "deltaGMean", "deltaGRange", "sequence",
        "sequenceRc", "gcContent", "tm", "deltaG", "method", "score",
        "roiStart", "roiEnd"
    )
    x <- x[keep]
    x <- x[order(x$start), ]
    rownames(x) <- NULL
    x
}
