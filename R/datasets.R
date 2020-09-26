#' An alignment of hepatitis E virus sequences
#'
#' An alignment with 621 hepatitis E virus sequences.
#'
#' @format A Biostrings::DNAMultipleAlignment object.
#'
#' @usage data(exampleRprimerAlignment)
#'
#' @source
#' The sequences were collected from NCBI GenBank and the alignment
#' was imported to the package using Biostrings::readDNAStringSet().
#'
#' @references
#' H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
"exampleRprimerAlignment"

#' A sequence profile from hepatitis E virus sequences
#'
#' @format An RprimerProfile object.
#'
#' @usage data(exampleRprimerProfile)
#'
#' @source
#' The sequences were collected from NCBI GenBank.
"exampleRprimerProfile"

#' Sequence properties of hepatitis E virus sequences
#'
#' @format An RprimerProperties object.
#'
#' \describe{
#'   \item{Position}{position in the alignment}
#'   \item{Majority}{majority consensus sequence}
#'   \item{IUPAC}{iupac consensus sequence}
#'   \item{Gaps}{proportion of gaps}
#'   \item{Identity}{proportion of the most common nucleotide}
#'   \item{Entropy}{Shannon entropy}
#' }
#'
#' @usage data(exampleRprimerProperties)
#'
#' @source The sequences were collected from NCBI GenBank.
"exampleRprimerProperties"
