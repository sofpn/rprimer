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

#' A sequence profile from hepatitis E virus sequences
#'
#' @format An \code{RprimerProfile} object.
#'
#' @usage data("exampleRprimerProfile")
#'
#' @source
#' The sequences were collected from NCBI GenBank.
"exampleRprimerProfile"

#' Sequence properties of hepatitis E virus sequences
#'
#' @format An \code{RprimerProperties} object.
#'
#' \describe{
#'   \item{Position}{Position in the alignment.}
#'   \item{Majority}{Majority consensus sequence.}
#'   \item{IUPAC}{IUPAC consensus sequence.}
#'   \item{Gaps}{Proportion of gaps.}
#'   \item{Identity}{Proportion of the most common nucleotide.}
#'   \item{Entropy}{Shannon entropy}
#' }
#'
#' @usage data("exampleRprimerProperties")
#'
#' @source The sequences were collected from NCBI GenBank.
"exampleRprimerProperties"
