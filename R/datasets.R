#' An multiple DNA alignment of hepatitis E virus sequences
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
#' was imported to the package using \code{Biostrings::readDNAStringSet()}
#' (Pages et al., 2020).
#'
#' @references
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
"exampleRprimerAlignment"

#' A consensus profile from an alignment of hepatitis E virus sequences
#'
#' @format An \code{RprimerProfile} object.
#'
#' @description
#' The object contains the following information:
#'
#' \describe{
#'   \item{position}{Position in the alignment.}
#'   \item{a}{Proportion of A.}
#'   \item{c}{Proportion of C.}
#'   \item{g}{Proportion of G.}
#'   \item{t}{Proportion of T.}
#'   \item{other}{Proportion of bases other than A, C, G, T.}
#'   \item{gaps}{Proportion of gaps.}
#'   \item{majority}{Majority consensus sequence.}
#'   \item{identity}{Proportion of the most frequently occurring nucleotide.}
#'   \item{iupac}{IUPAC consensus character.}
#'   \item{entropy}{Shannon entropy.}
#'   \item{coverage}{Proportion of sequences covered by the
#'   IUPAC consensus character.}
#' }
#'
#' @usage data("exampleRprimerProfile")
"exampleRprimerProfile"

#' Primers and probes targeted to hepatitis E virus
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
#'   \item{iupacSequence}{Oligo sequence in IUPAC format
#'   (i.e. with ambiguous bases).}
#'   \item{iupaSequenceRc}{The reverse complement of the iupacSequence.}
#'   \item{identity}{Average identity of the oligo.}
#'   \item{coverage}{Average coverage of the oligo.}
#'   \item{degeneracy}{Number of sequence variants of the oligo.}
#'   \item{gcContentMean}{Mean GC-content of all sequence variants of the oligo.
#'   }
#'   \item{gcContentRange}{Range in GC-content of all sequence variants of
#'     the oligo.}
#'   \item{tmMean}{Mean tm of all sequence variants of the oligo.}
#'   \item{tmRange}{Range in tm of all sequence variants of the oligo.}
#'   \item{sequence}{All sequence variants of the oligo.}
#'   \item{sequenceRc}{Reverse complements of all sequence variants.}
#'   \item{gcContent}{GC-content of all sequence variants.}
#'   \item{tm}{Tm of all sequence variants (in Celcius degrees).}
#'   \item{dH}{delta H of all sequence variants (in cal/mol).}
#'   \item{dS}{delta S of all sequence
#'   variants (in cal/K/mol).}
#'   \item{method}{Design method used to generate the oligo: "ambiguous",
#'   "mixedFwd" or "mixedRev".}
#'   \item{score}{Oligo score, the lower the better.
#'   See "Score" for more details.}
#'   \item{roiStart}{First position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd}{Last position of the input \code{RprimerProfile} object.}
#' }
#'
#' Tm was calculated for perfectly matching oligo-target duplexes at a
#' sodium ion
#' concentration of 0.05 M, a primer
#' concentration of 500 nM and a probe concentration of 250 nM.
#'
#' Delta H was calculated for perfectly matching oligo-target duplexes at a
#' sodium ion
#' concentration of 0.05 M.
#'
#' @usage data("exampleRprimerOligo")
"exampleRprimerOligo"

#' PCR assays targeted to hepatitis E virus
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
#'   \item{score}{Summarized oligo score.}
#'   \item{plusPr}{If the probe is valid in positive sense.}
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
#'   \item{sequenceFwd}{Sequence of the forward primer, all sequence variants.}
#'   \item{gcContentFwd}{GC-content of the forward primer,
#'   all sequence variants.}
#'   \item{tmFwd}{Tm of the forward primer, all sequence variants.}
#'   \item{dHFwd}{delta H of the forward
#'   primer (in cal/mol), all sequence variants.}
#'   \item{dSFwd}{delta S of
#'   the forward primer (in cal/K/mol), all sequence variants.}
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
#'   \item{sequenceRev}{Sequence of the reverse primer, all sequence variants.}
#'   \item{gcContentRev}{GC-content of the reverse primer, all sequence variants.}
#'   \item{tmRev}{Tm of the reverse primer, all sequence variants.}
#'   \item{dHRev}{delta H of the reverse
#'   primer (in cal/mol), of all sequence variants.}
#'   \item{dSRev}{delta S of the reverse primer (in cal/K/mol),
#'   all sequence variants.}
#'   \item{plusPr}{If the probe is valid in negative sense.}
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
#'   \item{sequencePr}{Sequence of the probe, all sequence variants.}
#'   \item{gcContentPr}{GC-content of the probe, all sequence variants.}
#'   \item{tmPr}{Tm of the probe, all sequence variants.}
#'   \item{dHPr}{delta H of the probe (in cal/mol), all sequence variants.}
#'   \item{dSPr}{delta S of the probe (in cal/K/mol), all sequence variants.}
#'   \item{roiStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{roiEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' Tm was calculated for perfectly matching oligo-target duplexes at a
#' sodium ion
#' concentration of 0.05 M, a primer
#' concentration of 500 nM and a probe concentration of 250 nM.
#'
#' Delta H was calculated for perfectly matching oligo-target duplexes at a
#' sodium ion
#' concentration of 0.05 M.
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"

#' Proportion of matching sequences to a target
#'
#' This dataset describes the proportion of oligo sequences that
#' matches to a set of HEV target sequences.
#'
#' @format An \code{RprimerMatchOligo} object.
#'
#' @description
#' The object contains the following information:
#'
#' \describe{
#'   \item{iupacSequence}{The oligo sequence in IUPAC format.}
#'   \item{perfectMatch}{Proportion of target sequences that matches perfectly
#'   to the oligo within the intended oligo binding region.}
#'   \item{oneMismatch}{Proportion of target sequences with one mismatch to
#'   the oligo
#'   within the intended oligo binding region.}
#'   \item{twoMismatches}{Proportion of target sequences with two mismatches
#'   to the oligo
#'   within the intended oligo binding region.}
#'   \item{threeMismatches}{Proportion of target sequences with
#'   three mismatches to the oligo
#'   within the intended oligo binding region.}
#'   \item{fourOrMoreMismatches}{Proportion of target sequences with
#'   four or more mismatches to the oligo
#'   within the intended oligo binding region.}
#'   \item{offTargetMatch}{Proportion of target sequences that matches
#'   to the oligo with
#'   no more than four mismatches to all other regions of the target alignment.}
#'  }
#'
#' @usage data("exampleRprimerMatchOligo")
"exampleRprimerMatchOligo"
