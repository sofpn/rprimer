#' Example datasets
#'
#' @name example-datasets
#'
#' @description
#' The purpose of these datasets is to illustrate the functionality of
#' rprimer. The following datasets are provided:
#'
#' \itemize{
#' \item \code{exampleRprimerAlignment} - a \code{Biostrings::DNAMultipleAlignment}
#' object (Pages et al., 2020)
#' containing an alignment of 198 hepatitis E virus
#' sequences collected from NCBI GenBank. See
#' "documentation_example_alignment.txt" within the inst/script folder of this
#' package for more details.
#' \item \code{exampleRprimerProfile} - an \code{RprimerProfile} object, generated
#' from the alignment above.
#' \item \code{exampleRprimerOligo} - an \code{RprimerOligo} object, generated from
#' the consensus profile above.
#' \item \code{exampleRprimerAssay} - an \code{RprimerAssay} object, generated from
#' the oligos above.
#' \item \code{exampleRprimerMatchOligo} - an \code{RprimerMatchOligo} object,
#' describing how well some oligos match with the sequences in
#' exampleRprimerAlignment.
#' \item \code{exampleRprimerMatchAssay} - an \code{RprimerMatchAssay} object,
#' describing how well some assays match with the seuqences in
#' exampleRprimerAlignment.
#' }
#'
#' @references
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
NULL

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerAlignment")
"exampleRprimerAlignment"

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerProfile")
"exampleRprimerProfile"

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerOligo")
"exampleRprimerOligo"

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerAssay")
"exampleRprimerAssay"

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerMatchOligo")
"exampleRprimerMatchOligo"

#' @rdname example-datasets
#'
#' @usage data("exampleRprimerMatchAssay")
"exampleRprimerMatchAssay"
