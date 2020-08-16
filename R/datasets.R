#' An alignment of hepatitis E virus sequences
#'
#' An alignment with 621 hepatitis E virus sequences.
#'
#' @format A list with class attribute 'rprimer_alignment'.
#'
#' @source The sequences were collected from NCBI GenBank,
#' by querying “[a-z]”[porgn:__txid1678143].
#' Only sequences with more than 6000 bases were included.
#' \url{https://www.ncbi.nlm.nih.gov/}
"example_rprimer_alignment"

#' A sequence profile from hepatitis E virus sequences
#'
#' @format A numeric matrix with class attribute 'rprimer__profile'.
#'
#' @source The sequences were collected from NCBI GenBank,
#' by querying “[a-z]”[porgn:__txid1678143].
#' Only sequences with more than 6000 bases were included.
#' \url{https://www.ncbi.nlm.nih.gov/}
"example_rprimer_profile"

#' Sequence properties of hepatitis E virus sequences
#'
#' @format A data frame (tibble) with class attribute ´
#' rprimer_properties'.
#' \describe{
#'   \item{position}{position in the alignment}
#'   \item{majority}{majority consensus sequence}
#'   \item{iupac}{iupac consensus sequence}
#'   \item{gaps}{proportion of gaps}
#'   \item{identity}{proportion of the most common nucleotide}
#'   \item{entropy}{Shannon entropy}
#' }
#'
#' @source The sequences were collected from NCBI GenBank,
#' by querying “[a-z]”[porgn:__txid1678143].
#' Only sequences with more than 6000 bases were included.
#' \url{https://www.ncbi.nlm.nih.gov/}
"example_rprimer_properties"


#' Oligos generated from hepatitis E virus sequences
#'
#' @format A data frame (tibble) with class attribute ´rprimer_oligo'.
#' \describe{
#'   \item{begin}{position where the oligo begins}
#'   \item{end}{position where the oligo ends}
#'   \item{length}{length of the oligo}
#'   \item{majority}{majority sequence}
#'   \item{iupac}{iupac sequence (i.e. with degenerate bases)}
#'   \item{majority_rc}{majority sequence, reverse complement}
#'   \item{iupac_rc}{iupac sequence, reverse complement}
#'   \item{degenerates}{number of degenerate bases}
#'   \item{degeneracy}{number of variants}
#'   \item{gc_majority}{gc-content (majority sequence), proportion}
#'   \item{tm_majority}{melting temperature (majority sequence),
#'   degrees Celcius}
#'   \item{pm_majority}{proportion of perfectly matching sequences,
#'   (to the target alignment) majority sequence}
#'   \item{pm_iupac}{proportion of perfectly matching sequences,
#'   iupac sequence}
#'   \item{match_matrix}{a logical matrix describing which sequences
#'   in the target alignment the oligo matches perfectly to}
#' }
#'
#' @source The sequences were collected from NCBI GenBank,
#' by querying “[a-z]”[porgn:__txid1678143].
#' Only sequences with more than 6000 bases were included.
#' \url{https://www.ncbi.nlm.nih.gov/}
"example_rprimer_oligo"

#' Assays generated from hepatitis E virus sequences
#'
#' @format A data frame (tibble) with class attribute ´rprimer_assay'.
#'  \describe{
#'   \item{begin}{position where the assay begins}
#'   \item{end}{position where the assay ends}
#'   \item{amplicon_length}{length of the amplicon}
#'   \item{tm_difference_primer}{difference in melting temperature between
#'   the forward and reverse primer, absolute value, degrees Celcius}
#'   \item{total_degeneracy}{total number of oligos in the assay}
#'   \item{begin_fwd}{position where the forward primer begins}
#'   \item{end_fwd}{position where the reverse primer ends}
#'   \item{length_fwd}{length of the forward primer}
#'   \item{majority_fwd}{majority sequence of the forward primer}
#'   \item{iupac_fwd}{iupac sequence (i.e. with degenerate bases)
#'   of the forward primer}
#'   \item{degenerates_fwd}{number of degenerate bases of the forward primer}
#'   \item{degeneracy_fwd}{number of variants of the forward primer}
#'   \item{gc_majority_fwd}{gc-content of the forward primer
#'   (majority sequence), proportion}
#'   \item{tm_majority_fwd}{melting temperature of the forward primer
#'   (majority sequence), degrees Celcius}
#'   \item{pm_majority_fwd}{proportion of sequences in the target alignment
#'   that matches perfectly with the forward primer, majority sequence}
#'   \item{pm_iupac_fwd}{proportion of sequences in the target alignment
#'   that matches perfectly with the forward primer, iupac sequence}
#'   \item{begin_rev}{position where the reverse primer begins}
#'   \item{end_rev}{position where the reverse primer ends}
#'   \item{length_rev}{length of the reverse primer}
#'   \item{majority_rev}{majority sequence of the reverse primer}
#'   \item{iupac_rev}{iupac sequence (i.e. with degenerate bases)
#'   of the reverse primer}
#'   \item{degenerates_rev}{number of degenerate bases of the reverse primer}
#'   \item{degeneracy_rev}{number of variants of the reverse primer}
#'   \item{gc_majority_rev}{gc-content of the reverse primer
#'   (majority sequence), proportion}
#'   \item{tm_majority_rev}{melting temperature of the reverse primer
#'   (majority sequence), degrees Celcius}
#'   \item{pm_majority_rev}{proportion of sequences in the target alignment
#'   that matches perfectly with the reverse primer, majority sequence}
#'   \item{pm_iupac_rev}{proportion of sequences in the target alignment
#'   that matches perfectly with the reverse primer, iupac sequence}
#'   \item{pm_majority_all}{proportion of sequences in the target alignment
#'   that matches perfectly with both the forward and reverse primer,
#'   majority sequence}
#'   \item{pm_iupac_all}{proportion of sequences in the target alignment
#'   that matches perfectly with both the forward and reverse primer,
#'   iupac sequence}
#'   \item{match_matrix}{a logical matrix describing which sequences
#'   in the target alignment the assay matches perfectly to}
#' }
#'
#' @source The sequences were collected from NCBI GenBank,
#' by querying “[a-z]”[porgn:__txid1678143].
#' Only sequences with more than 6000 bases were included.
#' \url{https://www.ncbi.nlm.nih.gov/}
"example_rprimer_assay"
