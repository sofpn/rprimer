#' An alignment of hepatitis E virus sequences
#'
#' An alignment with 200 hepatitis E virus sequences (see the file
#' "documentation_example_alignment.txt" within the inst/script folder
#' for more details).
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
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
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
#'   \item{position}{Position in the alignment. Note that masked columns
#'   from the original alignment are removed, and
#'   hence not taken into account when position is determined.}
#'   \item{a}{Proportion of A.}
#'   \item{c}{Proportion of C.}
#'   \item{g}{Proportion of G.}
#'   \item{t}{Proportion of T.}
#'   \item{other}{Proportion of bases other than A, C, G, T.}
#'   \item{gaps}{Proportion of gaps (recognized as "-" in the alignment).}
#'   \item{majority}{Majority consensus sequence.
#'   The most frequently occurring nucleotide.
#'   If two or more bases occur with the same frequency,
#'   the consensus nucleotide will be randomly selected among these bases.}
#'   \item{identity}{Proportion of the most common nucleotide.
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from the
#'   calculation.}
#'   \item{iupac}{
#'   The consensus sequence expressed in IUPAC format.
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
#'   \item{coverage}{The proportion of bases which are included the
#'   consensus/ambiguous (IUPAC) base.
#'   Will be one if there are no "remaining" bases (and if
#'   \code{ambiguityThreshold = 0}).
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from the
#'   calculation.}
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
#'   \item{type}{Whether the oligo is a primer or probe.}
#'   \item{fwd}{\code{TRUE} if the oligo is valid in forward direction,
#'     \code{FALSE} otherwise.}
#'   \item{rev}{\code{TRUE} if the oligo is valid in reverse direction,
#'     \code{FALSE} otherwise.}
#'   \item{start}{Start position of the oligo.}
#'   \item{end}{End positon of the oligo.}
#'   \item{length}{Oligo length.}
#'   \item{iupacSequence}{Oligo sequence, with amiguous bases (if any).}
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
#'   \item{gcContentRange}{Range in GC-content of all sequence variants of
#'     the oligo.}
#'   \item{tmMean}{Mean tm of all sequence variants of the oligo.}
#'   \item{tmRange}{Range in tm of all sequence variants of the oligo.}
#'   \item{deltaGMean}{Mean delta G (in kcal/mole)
#'   of all sequence variants of the oligo.}
#'   \item{deltaGRange}{Range in delta G (in kcal/mole)
#'   of all sequence variants of the oligo.}
#'   \item{sequence}{All sequence variants of the oligo.}
#'   \item{sequenceRc}{Reverse complements of all sequence variants.}
#'   \item{gcContent}{GC-content of all sequence variants.}
#'   \item{tm}{Tm of all sequence variants.}
#'   \item{deltaG}{deltaG (in kcal/mole) of all sequence variants.}
#'   \item{method}{Design method used to generate the oligo: "ambiguous",
#'   "mixedFwd" or "mixedRev".}
#'   \item{roiStart}{First position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd}{Last position of the input \code{RprimerProfile} object.}
#' }
#'
#' Tm was calculated for a sodium ion concentration of 0.05 M, a primer
#' concentration of 500 nM and a probe concentration of 250 nM. Delta G was
#' calculated for a temperature of 60C, and a sodium ion concentration
#' of 0.05 M.
#'
#' @usage data("exampleRprimerOligo")
"exampleRprimerOligo"

#' PCR assays for hepatitis E virus
#'
#' @format An \code{RprimerAssay} object.
#'
#' @description
#' The object contains the following information:
#'
#' \describe{
#'   \item{start}{Position where the assay starts.}
#'   \item{end}{Position where the assay ends.}
#'   \item{ampliconLength}{Length of the amplicon.}
#'   \item{tmDifferencePrimer}{Difference in tm between
#'   the forward and reverse primer, absolute value.}
#'   \item{totalDegeneracy}{Total number of oligos in the assay.}
#'   \item{startFwd}{Position where the forward primer starts.}
#'   \item{endFwd}{Position where the forward primer ends.}
#'   \item{lengthFwd}{Length of the forward primer.}
#'   \item{iupacSequenceFwd}{IUPAC sequence of the forward primer.}
#'   \item{coverageFwd}{Average coverage of the forward primer.}
#'   \item{degeneracyFwd}{Number of variants of the forward primer.}
#'   \item{gcContentMeanFwd}{Mean GC-content of the forward primer.}
#'   \item{gcContentRangeFwd}{Range in GC-content of the forward primer.}
#'   \item{tmMeanFwd}{Mean tm of the forward primer.}
#'   \item{tmRangeFwd}{Range in tm of the forward primer.}
#'   \item{sequenceFwd}{Sequence of the forward primer, all variants.}
#'   \item{gcContentFwd}{GC-content of the forward primer, all variants.}
#'   \item{tmFwd}{Tm of the forward primer, all variants.}
#'   \item{startRev}{Position where the reverse primer starts.}
#'   \item{endRev}{Position where the reverse primer ends.}
#'   \item{lengthRev}{Length of the reverse primer.}
#'   \item{iupacSequenceRev}{IUPAC sequence of the reverse primer.}
#'   \item{coverageRev}{Average coverage of the reverse primer.}
#'   \item{degeneracyRev}{Number of variants of the reverse primer.}
#'   \item{gcContentMeanRev}{Mean GC-content of the reverse primer.}
#'   \item{gcContentRangeRev}{Range in GC-content of the reverse primer.}
#'   \item{tmMeanRev}{Mean tm of the reverse primer.}
#'   \item{tmRangeRev}{Range in tm of the reverse primer.}
#'   \item{sequenceRev}{Sequence of the reverse primer, all variants.}
#'   \item{gcContentRev}{GC-content of the reverse primer, all variants.}
#'   \item{tmRev}{Tm of the reverse primer, all variants.}
#'   \item{plusPr}{If the probe is valid in positive sense.}
#'   \item{minusPr}{If the probe is valid in negative sense.}
#'   \item{startPr}{Position where the probe starts.}
#'   \item{endPr}{Position where the probe ends.}
#'   \item{lengthPr}{Length of the probe.}
#'   \item{iupacSequencePr}{IUPAC sequence of the probe.}
#'   \item{coveragePr}{Average coverage of the probe.}
#'   \item{degeneracyPr}{Number of variants of the probe.}
#'   \item{gcContentMeanPr}{Mean GC-content of the probe.}
#'   \item{gcContentRangePr}{Range in GC-content of the probe.}
#'   \item{tmMeanPr}{Mean tm of the probe.}
#'   \item{tmRangePr}{Range in tm of the probe.}
#'   \item{sequencePr}{Sequence of the probe, all variants.}
#'   \item{gcContentPr}{GC-content of the probe, all variants.}
#'   \item{tmPr}{Tm of the probe, all variants.}
#'   \item{roiStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{roiEnd}{End position of the input consensus profile used
#'   for oligo design.}
#'  }
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"
