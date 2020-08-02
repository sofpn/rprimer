#' Read an alignment in fasta format
#'
#' \code{read_fasta_alignment} reads a .txt-file with an alignment of
#' DNA sequences in fasta format.
#'
#' @param x
#' The name of the file of which the alignment is to be read from
#' (a character vector of length one).
#'
#' @details
#' The file must be a .txt file and contain one or more aligned
#' DNA sequences in fasta format. Sequences in fasta format always
#' begins with a '>' symbol, followed by a single-line name, with the
#' sequence present on the next line.
#' All sequences (including gaps) must be of the same length.
#' Valid nucleotides are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
#' n', 'h', 'd', 'v', 'b' and '-'. The alignment can be in either
#' upper- or lowercase format.
#'
#' @note
#' This function is partly inspired by \code{seqinr::read.fasta}
#' (see references).
#'
#' @return
#' The alignment from the input
#' file (a named list with class attribute 'rprimer_alignment').
#' Each object (DNA sequence) has the same name as in the input file,
#' except for the '>' symbol.
#' Sequences are presented as character vectors of length one,
#' in lowercase format.
#'
#' @examples
#' read_fasta_alignment(
#' system.file('extdata', 'example_alignment.txt', package = 'rprimer')
#' )
#'
#' @references
#' Charif, D. and Lobry, J.R. (2007)
#' SeqinR 1.0-2: a contributed package to the R project for statistical
#' computing devoted to biological sequences retrieval and analysis.
#' Structural approaches to sequence evolution: Molecules, networks,
#' populations, pp 207-232.
#'
#' @export
read_fasta_alignment <- function(x) {
  if (!is.character(x)) {
    stop("'filename' must be a character vector", call. = FALSE)
  }
  # Check if we have permission to read the file named x, and stop if not
  access <- file.access(x, mode = 4)
  if (access != 0) {
    stop(paste("File", x, "was not found/is not readable"), call. = FALSE)
  }
  # Import the fasta file
  infile <- readLines(x, warn = FALSE)
  # Identify where all the sequence names are
  index <- grep(">", infile)
  # Stop if the file does not appear to be in fasta format
  if (length(index) == 0L) {
    stop(
      "The file does not appear to be in fasta format. \n
      No line starts with '>'.", call. = FALSE
    )
  }
  # Identify where the sequences start
  begin <- index + 1
  # Identify where the sequences end
  end <- index - 1
  # Remove the first 'end', as there is no sequence above the first sequence
  end <- c(end[-1], length(infile))
  # Make a list with all sequences
  sequences <- purrr::map(seq_along(index), function(i) {
    sequence <- paste(infile[begin[[i]]:end[[i]]], collapse = "")
    sequence <- tolower(sequence)
    sequence
  })
  # Get the name of each sequence
  name <- infile[index]
  name <- gsub(">", "", name)
  names(sequences) <- name
  # Assign and validate a class attribute (rprimer_alignment)
  sequences <- new_rprimer_alignment(sequences)
  sequences
}
