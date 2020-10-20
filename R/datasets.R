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
#' }
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"
