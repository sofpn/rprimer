#' An multiple DNA alignment of hepatitis E virus sequences
#'
#' An alignment with 200 hepatitis E virus sequences.
#'
#' @format A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @usage data("exampleRprimerAlignment")
#'
#' @source
#' The sequences were collected from NCBI GenBank and the alignment
#' was imported to the package using \code{Biostrings::readDNAStringSet()}
#' (Pages et al., 2020). See the file
#' "documentation_example_alignment.txt" within the inst/script folder
#' for more details on how it was generated.
#'
#' @references
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
"exampleRprimerAlignment"

#' A consensus profile from an alignment of hepatitis E virus sequences
#'
#' @format An example of an \code{RprimerProfile} object.
#'
#' @usage data("exampleRprimerProfile")
"exampleRprimerProfile"

#' Primers and probes targeted to hepatitis E virus
#'
#' @format An example of an \code{RprimerOligo} object.
#'
#' @describeIn exampleData
#'
#' @usage data("exampleRprimerOligo")
"exampleRprimerOligo"

#' PCR assays targeted to hepatitis E virus
#'
#' @format An example of an \code{RprimerAssay} object.
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"

#' Proportion of sequences that matches to oligos
#'
#' This dataset describes the proportion of HEV target sequences that
#' matches to a set of oligos.
#'
#' @format An example of an \code{RprimerMatchOligo} object.
#'
#' @usage data("exampleRprimerMatchOligo")
"exampleRprimerMatchOligo"


#' Proportion of sequences that matches to assays
#'
#' This dataset describes the proportion of HEV target sequences that
#' matches to a set of assays.
#'
#' @format An example of an \code{RprimerMatchAssay} object.
#'
#' @usage data("exampleRprimerMatchAssay")
"exampleRprimerMatchAssay"
