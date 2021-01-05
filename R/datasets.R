#' An alignment of hepatitis E virus sequences
#'
#' An alignment with 100 hepatitis E virus sequences.
#'
#' @format A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @usage data("exampleRprimerAlignment")
#'
#' @source
#' The sequences were collected from NCBI GenBank and the alignment
#' was imported to the package using \code{Biostrings::readDNAStringSet()}.
#'
#' @references
#' H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
"exampleRprimerAlignment"

#' Consensus profile of an alignment of hepatitis E virus sequences
#'
#' @format An \code{RprimerProfile} object.
#'
#' @description
#'
#' \describe{
#'   \item{position}{Position in the alignment.}
#'   \item{a}{Proportion of A.}
#'   \item{c}{Proportion of C.}
#'   \item{g}{Proportion of G.}
#'   \item{t}{Proportion of T.}
#'   \item{other}{Proportion of bases other than A, C, G, T.}
#'   \item{gaps}{proportion of gaps (recognized as "-" in the alignment).}
#'   \item{majority}{Majority consensus sequence.
#'   The most frequently occurring nucleotide.}
#'   \item{identity}{Proportion of the most common nucleotide.
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from this
#'   calculation.}
#'   \item{iupac}{
#'   The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
#'   Note that the IUPAC consensus sequence only
#'   takes 'A', 'C', 'G', 'T' and '-' as input. Degenerate bases
#'   present in the alignment will be skipped. If a position only contains
#'   degenerate/invalid bases, the IUPAC consensus will be \code{NA} at that
#'   position.}
#'   \item{entropy}{Shannon entropy.
#'   Shannon entropy is a measurement of
#'   variability.
#'   First, for each nucleotide that occurs at a specific position,
#'   \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#'   that nucleotide. Then, these values are summarized,
#'   and multiplied by \code{-1}.
#'   A value of \code{0} indicate no variability and a high value
#'   indicate high variability.
#'   Gaps (-), as well as bases other than
#'   A, C, G and T are excluded from the calculation.}
#' }
#'
#' @usage data("exampleRprimerProfile")
#'
#' @source The sequences were collected from NCBI GenBank.
"exampleRprimerProfile"

#' Oligonucleotides for hepatitis E virus
#'
#' @format An \code{RprimerOligo} object.
#'
#' @description
#' The object contains the following information:
#'
#' \describe{
#'   \item{type: whether the oligo is a primer or probe.}
#'   \item{fwd: \code{TRUE} if the oligo is valid in forward direction,
#'     \code{FALSE} otherwise.}
#'   \item{rev: \code{TRUE} if the oligo is valid in reverse direction,
#'     \code{FALSE} otherwise.}
#'   \item{start: start position of the oligo.}
#'   \item{end: end positon of the oligo.}
#'   \item{length: oligo length.}
#'   \item{iupacSequence: oligo sequence, with wobble bases (if any).}
#'   \item{iupaSequenceRc: the reverse complement of the iupacSequence.}
#'   \item{identity: average identity score of the oligo, can range from 0 to 1.
#'     The identity is the proportion of the most common base at each position
#'     in the input alignment.}
#'   \item{degeneracy: number of sequence variants of the oligo.}
#'   \item{gcContentMean: mean GC-content of all sequence variants of the oligo.
#'   }
#'   \item{gcContentRange: range in GC-content of all sequence variants of
#'     the oligo.}
#'   \item{tmMean: mean tm of all sequence variants of the oligo.}
#'   \item{tmRange: range in tm of all sequence variants of the oligo.}
#'   \item{sequence: all sequence variants of the oligo.}
#'   \item{sequenceRc: reverse complements of all sequence variants.}
#'   \item{gcContent: GC-content of all sequence variants.}
#'   \item{tm: tm of all sequence variants.}
#'   \item{roiStart: first position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd: last position of the input \code{RprimerProfile} object}
#' }
#'
#' @usage data("exampleRprimerOligo")
"exampleRprimerOligo"


#' PCR assays for hepatitis E virus
#'
#' @format An \code{RprimerAssay} object.
#'
#' @description The object contains the following information:
#'
#' \describe{
#'   \item{start}{Position where the assay starts.}
#'   \item{end}{Position where the assay ends.}
#'   \item{ampliconLength}{Length of the amplicon.}
#'   \item{tmDifferencePrimer}{Difference in Tm between
#'   the forward and reverse primer, absolute value.}
#'   \item{meanIdentity}{Average identity score of the primers
#'   (and probe if selected)}.
#'   \item{totalDegeneracy}{Total number of oligos in the assay.}
#'   \item{startFwd}{Position where the forward primer starts.}
#'   \item{endFwd}{Position where the reverse primer ends.}
#'   \item{lengthFwd}{Length of the forward primer.}
#'   \item{majorityFwd}{Majority sequence of the forward primer.}
#'   \item{gcMajorityFwd}{GC-content of the forward primer
#'   (majority sequence), proportion.}
#'   \item{identityFwd}{Average identity of the forward primer.}
#'   \item{tmMajorityFwd}{Tm of the forward primer
#'   (majority sequence).}
#'   \item{iupacFwd}{IUPAC sequence (i.e. with degenerate bases)
#'   of the forward primer.}
#'   \item{degeneracyFwd}{Number of variants of the forward primer}.
#'   \item{startRev}{Position where the reverse primer starts.}
#'   \item{endRev}{Position where the reverse primer ends.}
#'   \item{lengthRev}{Length of the reverse primer.}
#'   \item{majorityRev}{Majority sequence of the reverse primer.}
#'   \item{gcMajorityRev}{GC-content of the reverse primer
#'   (majority sequence), proportion.}
#'   \item{tmMajorityRev}{Tm of the reverse primer
#'   (majority sequence).}
#'   \item{identityRev}{Average identity of the reverse primer.}
#'   \item{iupacRev}{IUPAC sequence (i.e. with degenerate bases)
#'   of the reverse primer.}
#'   \item{degeneracyRev}{Number of variants of the reverse primer.}
#'   \item{allFwd}{Lists with all sequence variants of the forward primer.}
#'   \item{gcAllFwd}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{tmAllFwd}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#'   \item{allRev}{Lists with all sequence variants of the reverse primer.}
#'   \item{gcAllRev}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{tmAllRev}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#'   \item{tmDifferencePrimerProbe}{Difference in Tm between the average
#'   Tm of the primer pair and the probe, majority sequences.}
#'   \item{startPr}{Position where the probe starts.}
#'   \item{endPr}{Position where the probe ends.}
#'   \item{lengthPr}{Length of the probe.}
#'   \item{majorityPr}{Majority sequence of the probe.}
#'   \item{gcMajorityPr}{GC-content of the probe
#'   (majority sequence), proportion.}
#'   \item{tmMajorityPr}{Tm of the probe
#'   (majority sequence).}
#'   \item{identityPr}{Average identity of the probe.}
#'   \item{iupacPr}{IUPAC sequence (i.e. with degenerate bases)
#'   of the probe.}
#'   \item{degeneracyPr}{Number of variants of the probe.}
#'   \item{sensePr}{Sense of the probe (pos or neg). If both probes are valid,
#'   the probe with the least G:s is selected.}
#'   \item{allPr}{Lists with all sequence variants of the probe.}
#'   \item{gcAllPr}{Lists with the GC content of all
#'   sequence variants of the probe.}
#'   \item{tmAllPr}{Lists with the Tm of all sequence variants of
#'   the probe.}
#'   \item{alignmentStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{alingnmentEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"
