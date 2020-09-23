
#' Check if an object is as RprimerProperties-object
#'
#' @param x An 'RprimerProperties'-like object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
#'
#' @noRd
is.RprimerProperties <- function(x) inherits(x, "RprimerProperties")


#mask
#as matrix

# matchProbePair function in the Biostrings package.
#the context of a computer-simulated PCR experiment, one wants
#to find the amplicons mapped to a given primer pair. The
#'matchProbePair' function can be used for this: given a forward
#and a reverse probe (i.e. the chromosome-specific sequences of the
#                     forward and reverse primers used for the experiment) and a target
#sequence (generally a chromosome sequence), the 'matchProbePair'
#function will return all the "theoretical amplicons" mapped to
#this probe pair.

#Note that matchProbePair will only return the locations of the amplicons so,
#in order to complete the validation process, some extra work is required to
#map those locations to the features that you are targetting.
