# Import alignment ============================================================

#' Read an alignment in fasta format
#'
#' \code{read_fasta_alignment} reads a fasta file with aligned DNA
#' sequences.
#'
#' @param x The name of the file of which the
#' alignment is to be read from (a character vector of length one).
#'
#' @details The file must contain one or more aligned
#' DNA sequences in fasta format. Sequences in fasta format
#' begins with a '>' symbol, followed by a single-line name.
#' The sequence is present on the next line.
#' All sequences (including gaps) must be of the same length,
#' and the minimum length of the alignment is 200.
#' Valid nucleotides are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
#' n', 'h', 'd', 'v', 'b' and '-'. The alignment can be in either
#' upper- or lowercase format.
#'
#' @return The alignment from the input
#' file (a named list with class attribute 'rprimer_alignment').
#' Each object (DNA sequence) will have the same name as its
#' sequence name in the input file, except for the '>' symbol.
#' Sequencees are will be presented as character vectors of length one,
#' in lowercase format.
#'
#' @examples
#' read_fasta_alignment(
#' system.file('extdata', 'example_alignment.txt', package = 'rprimer')
#' )
#'
#' @export
read_fasta_alignment <- function(x) {
    if (!is.character(x) || length(x) != 1) {
        stop(
          "A filename (character vector of length one) is expected",
          call. = FALSE
        )
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
          "The file does not appear to be in fasta format
          (no line starts with '>').", call. = FALSE
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
        return(sequence)
    })
    # Get the name of each sequence
    name <- infile[index]
    name <- gsub(">", "", name)
    names(sequences) <- name
    # Assign and validate a class attribute (rprimer_alignment)
    sequences <- new_rprimer_alignment(sequences)
    return(sequences)
}

#' Remove positions with high gap frequency
#'
#' \code{remove_gaps} removes all positions with a gap frequency
#' higher than a stated threshold in an alignment of DNA sequences.
#'
#' @param x An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @param threshold A number between 0.5 and <1 (the default is 0.5).
#'
#' @details Gaps are recognised as "-". The alignment (with gaps removed)
#' must contain at least 200 bases (an error message will return if not).
#' Note that the positions will
#' not be kept from the input alignment. The positions in the new
#' alignment will always start at 1, and increase by 1 for every new
#' position.
#'
#' @return An alignment (an list of class 'rprimer_alignment'), where
#' all positions with a gap proportion higher than the stated threshold
#' are removed.
#'
#' @examples
#' remove_gaps(example_rprimer_alignment, threshold = 0.5)
#'
#' @export
#'
#' @seealso \code{read_fasta_alignment}
remove_gaps <- function(x, threshold = 0.5) {
    if (!(threshold >= 0.5 && threshold < 1)) {
        stop(paste0(
          "The threshold was set to ", threshold, ".
          treshold must be between 0.5 and <1"
        ), call. = FALSE)
    }
    if (!inherits(x, "rprimer_alignment")) {
        stop("An rprimer_alignment object is expected.", call. = FALSE)
    }
    splitted <- purrr::map(x, split_sequence)
    # Make a matrix
    matr <- do.call("rbind", splitted)
    # gaps will be represented as 1
    matr <- gsub("-", 1, matr)
    # Other characters will be represented as 0
    matr <- gsub("[a-z]", 0, matr)
    # Convert to integer
    matr <- apply(matr, 2, as.integer)
    # Keep positions with gap frequency below or at threshold
    splitted <- purrr::map(splitted, ~.x[which(colMeans(matr) <= threshold)])
    # And paste them together
    x <- purrr::map(splitted, paste, collapse = "")
    # Set class attribute
    x <- new_rprimer_alignment(x)
    return(x)
}

#' Select region of interest
#'
#' \code{select_roi} selects a specified region of interest within
#' an alignment of DNA sequences.
#'
#' @param  x An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @param from Where the roi begins (an integer). The default is 1.
#'
#' @param to Where the roi ends (an integer). The default is \code{NULL}.
#' In that case, \code{to} will be set as the last position in the alignment.
#'
#' @details The roi must be at least 200 bases (an error will return if not).
#' Note that the positions will not be kept from the input alignment.
#' The new alignment will always start at position 1.
#'
#' @return The roi of the alignment
#' (an object of class 'rprimer_alignment').
#'
#' @examples
#' # Select the first 1000 bases
#' select_roi(example_rprimer_alignment, from = 1, to = 1000)
#'
#' @export
select_roi <- function(x, from = 1, to = NULL) {
    if (!inherits(x, "rprimer_alignment")) {
        stop("An rprimer_alignment object is expected.", call. = FALSE)
    }
    splitted <- purrr::map(x, split_sequence)
    # Get the length of the alignment
    # All sequences are of equal length in an rprimer_alignment object
    aln_length <- length(splitted[[1]])
    if (to > aln_length) to <- NULL
    if (is.null(to)) to <- aln_length
    if (!(is.numeric(from) && is.numeric(to))) {
      stop("from and to must be positive integers", call. = FALSE)
    }
    if (to < from) {
      stop("to must be greater than from", call. = FALSE)
    }
    if (!(from >= 1 && to >= 1)) {
      stop("from and to must be positive integers", call. = FALSE)
    }
    if (!(from <= aln_length - 199)) {
      stop("The roi must be at least 200 bases.", call. = FALSE)
    }
    if (!(to - from >= 199)) {
      stop("The roi must be at least 200 bases.", call. = FALSE)
    }
    splitted <- purrr::map(splitted, ~.x[from:to])
    x <- purrr::map(splitted, paste, collapse = "")
    x <- new_rprimer_alignment(x)
    return(x)
}

# Get sequence profile ========================================================

#' Get the sequence profile of an alignment
#'
#' \code{sequence_profile} returns a matrix with the
#' proportion of each nucleotide at each position within an alignment
#' of DNA sequences.
#'
#' @param x An alignment of DNA sequences (an object of class
#' 'rprimer_alignment').
#'
#' @return The sequence profile (an object of class
#' 'rprimer_sequence_profile'). A numeric m x n matrix,
#' where m is the number of unique bases in the alignment, and n is the
#' number of positions in the alignment.
#'
#' @examples
#' sequence_profile(example_rprimer_alignment)
#'
#' @export
sequence_profile <- function(x) {
    if (!inherits(x, "rprimer_alignment")) {
        stop("An rprimer_alignment object is expected.", call. = FALSE)
    }
    splitted <- purrr::map(x, split_sequence)
    # Make a matrix
    matr <- do.call("rbind", splitted)
    # Get all unique bases in the dataset and sort them in alphabetical order
    bases <- unique(sort(unlist(apply(matr, 1, unique), use.names = FALSE)))
    # Count the occurence of each base at each position
    count_base <- function(x, base) length(x[which(x == base)])
    counts <- apply(matr, 2, function(x) {
        purrr::map_int(bases, ~count_base(x, base = .x))
    })
    # Present the data as proportions instead of counts
    proportions <- counts / nrow(matr)
    # Base as rowname
    rownames(proportions) <- bases
    # Nucleotide positon as colname
    colnames(proportions) <- seq_len(dim(matr)[[2]])
    proportions <- new_rprimer_sequence_profile(proportions)
    return(proportions)
}

# Get sequence properties from a nucleotide profile ===========================

#' Get sequence properties
#'
#' \code{sequence_properties} returns sequence information from an alignment
#' of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @param iupac_threshold A number between 0 and 0.2 (the default is 0).
#' At each position, all nucleotides with a proportion
#' higher than or equal to the stated threshold will be included in
#' the iupac consensus sequence.
#'
#' @section Majority consensus sequence:
#' The most frequently occuring base at each position.
#' If two or more bases occur with the same frequency,
#' the consensus base will be randomly selected among these bases.
#'
#' @section IUPAC consensus sequence:
#' The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
#' Note that the iupac consensus sequence only
#' takes 'a', 'c', 'g', 't' and '-' as input. Degenerate bases
#' present in the alignment will be skipped. If a position only contains
#' degenerate/invalid bases, the iupac consensus will be \code{NA} at that
#' position.
#'
#' @section Gaps:
#' Gaps are recognised as "-".
#'
#' @section Identity:
#' The nucleotide identity is the proportion of
#' the most common base at each position in the alignment.  Gaps (-),
#' as well as nucleotides other than a, c, g and t, are excluded from the
#' calculation.
#'
#' @section Shannon entropy:
#' Shannon entropy is a measurement of
#' variability. First, for each nucleotide that occurs at the position in
#' matter, \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#' that nucleotide. Then, the shannon entropy is calculated by summarising
#' these values for each nucleotide at the position in matter,
#' followed by multiplication by \code{-1}.
#' A value of \code{0} indicate no variability and a high value
#' indicate high variability.
#' Gaps (-), as well as bases other than
#' a, c, g and t, are excluded from the calculation.
#'
#' @return A tibble (data frame) of class 'rprimer_sequence_properties',
#' with information about majority and iupac consensus sequence, gap frequency,
#' nucleotide identity and shannon entropy.
#'
#' @examples
#' sequence_properties(example_rprimer_sequence_profile)
#'
#' @export
sequence_properties <- function(x, iupac_threshold = 0) {
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    position <- seq_len(ncol(x))
    majority <- majority_consensus(x)
    iupac <- iupac_consensus(x, threshold = iupac_threshold)
    gaps <- gap_frequency(x)
    identity <- nucleotide_identity(x)
    entropy <- shannon_entropy(x)
    sequence_properties <- tibble::tibble(
      position, majority, iupac, gaps, identity, entropy
    )
    sequence_properties <- tibble::new_tibble(
      sequence_properties, nrow = nrow(sequence_properties),
      class = "rprimer_sequence_properties"
    )
    return(sequence_properties)
}

# Get oligos from sequence data ===============================================

#' Get oligos from sequence data
#'
#' \code{get_oligos} identifies oligos (primers and probes) from sequence data
#' based on user-specfied constraints.
#'
#' @param x an object of class 'rprimer_sequence_properties'.
#'
#' @param max_gap_frequency Maximum allowed gap frequency.
#' A number between 0 and 1 (default is 0.1, which means that only
#' positions with a gap frequency equal to or less than 0.1 will be
#' considered as an oligo region).
#'
#' @param length Oligo length. One or more numbers. The minimum allowed
#' value is 6 and the maximum allowed value is 30.
#' The default is \code{18:22}.
#'
#' @param max_degenerates Maximum number of degenerate positions in each
#' oligo. An integer. The minimum allowed value is 0 and the maximum
#' allowed value is 6 (default is 2).
#'
#' @param max_degeneracy Maximum number of degenerate variants of each
#' oligo. A number. The minimum allowed value is 1 and the maximum
#' allowed value is 64 (default is 4).
#'
#' @param avoid_3end_ta \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with NA. The default is \code{TRUE}.
#'
#' @param avoid_5end_g \code{TRUE} or \code{FALSE}.
#' If oligos with g at the 5' end
#' should be replaced with NA. The default is \code{FALSE}.
#'
#' @param avoid_3end_runs \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with NA.
#' The default is \code{TRUE}.
#'
#'  @param gc_range The GC-content range of each oligo. Can range between
#' 0 and 1. The default is \code{c(0.45, 0.55)}.
#'
#' @param tm_range The Tm-range of each oligo. Can range between 20 and 90.
#' The default is \code{c(48, 70)}.
#'
#' @param conc_oligo The concentration of oligonucleotide in M,
#' ranging from 0.2e-07 M (20 nM) to 2e-06 M (2000 nM).
#' The default value is 5e-07 M (500 nM) (for Tm calculation)
#'
#' @param conc_na The sodium ion concentration in M, ranging
#' from 0.01 M to 1 M. The default value is 0.05 M (50 mM)
#' (for Tm calculation).
#'
#' @section Excluded oligos:
#' The function excludes oligos with
#' more than than three consecutive runs of the same dinucleotide
#' (e.g. 'tatatata'), oligos with more than four consecutive runs of the
#' same nucleotide, and oligos that are duplicated.
#'
#' @section Tm:
#' The melting temperature is calculated using the nearest-neigbour method.
#' The oligo concentration is set to 500 nM and the sodium ion concentration
#' is set to 50 mM.
#'
#' Assumptions for Tm calculation:
#'
#' Oligos are not expected to be self-complementary, so no symmetry
#' correction is done.
#'
#' We assume that the oligo concentration is much higher
#' than the target concentration.
#'
#' See references for table values and equations.
#'
#' #Warning:
#' GC-content and Tm are calculated based on the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return A tibble (a data frame) of class 'rprimer_oligo', with all oligo
#' candidates. An error message will return if no oligos are found.
#'
#' @examples
#' get_oligos <- function(
#' example_rprimer_sequence_properties,
#' length = 18:22,
#' max_gap_frequency = 0.1,
#' max_degenerates = 2,
#' max_degeneracy = 4,
#' avoid_3end_ta = TRUE,
#' avoid_5end_g = FALSE,
#' avoid_3end_runs = TRUE,
#' gc_range = c(0.45, 0.55),
#' tm_range = c(48, 70)
#' )
#'
#' @references
#' Tm-calculation:
#'
#' SantaLucia, J, et al. (1996)
#' Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability.
#' Biochemistry, 35: 3555-3562 (Formula and salt correction are from here)
#'
#' Allawi, H. & SantaLucia, J. (1997)
#' Thermodynamics and NMR of Internal G·T Mismatches in DNA.
#' Biochemistry, 34: 10581?\200?10594
#' (Duplex initiation parameters are from here)
#'
#' #' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbour thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @export
get_oligos <- function(x, length = 18:22, max_gap_frequency = 0.1, max_degenerates = 2, max_degeneracy = 4, avoid_3end_ta = TRUE,
    avoid_5end_g = FALSE, avoid_3end_runs = TRUE, gc_range = c(0.45, 0.55), tm_range = c(48, 70), conc_oligo = 5e-07, conc_na = 0.05) {
    if (!inherits(x, "rprimer_sequence_properties")) {
        stop("An rprimer_sequence_properties object is expected for x.", call. = FALSE)
    }
    if (!(min(length) >= 6 && max(length) <= 30 && is.numeric(length))) {
        stop("length must be between 4 and 30", call. = FALSE)
    }
    if (!(max_gap_frequency >= 0 && max_gap_frequency <= 1 && is.numeric(max_gap_frequency) && length(max_gap_frequency) ==
        1)) {
        stop("max_gap_frequency must be between 0 and 1", call. = FALSE)
    }
    if (!(max_degenerates <= 30 && max_degenerates >= 0 && is.numeric(max_degenerates))) {
        stop("max_degenerates must be between 0 and 6", call. = FALSE)
    }
    if (!(max_degeneracy >= 1 && max_degeneracy <= 64 && length(max_degeneracy == 1) && is.numeric(max_degeneracy))) {
        stop("max_degeneracy must be between 1 and 64", call. = FALSE)
    }
    if (!(min(gc_range) >= 0 && max(gc_range) <= 1 && is.numeric(gc_range))) {
        stop("gc_range must be between 0 and 1, e.g. c(0.45, 0.65)", call. = FALSE)
    }
    if (!(min(tm_range) >= 20 && max(tm_range) <= 90 && is.numeric(tm_range))) {
        stop("tm_range must be between 20 and 90, e.g. c(55, 60)", call. = FALSE)
    }
    if (!is.double(conc_oligo) || conc_oligo < 2e-08 || conc_oligo > 2e-06) {
        stop("The oligo concentration must be between
         0.2e-07 M (20 nM) and 2e-08 M (2000 nM)", call. = FALSE)
    }
    if (!is.double(conc_na) || conc_na < 0.01 || conc_na > 1) {
        stop("The Na+ concentration must be between 0.01 and 1 M", call. = FALSE)
    }
    all_oligos <- purrr::map_dfr(length, function(y) {
        # Find all possible oligos of length y
        majority <- get_nmers(x$majority, n = y)
        iupac <- get_nmers(x$iupac, n = y)
        majority_rc <- purrr::map_chr(majority, ~reverse_complement(.x))
        iupac_rc <- purrr::map_chr(iupac, ~reverse_complement(.x))
        degenerates <- purrr::map_int(iupac, ~count_degenerates(.x))
        degeneracy <- purrr::map_dbl(iupac, ~count_degeneracy(.x))
        begin <- seq(1, length(majority))
        end <- seq(y, length(majority) + y - 1)
        length <- y
        # Combine the information in a data frame
        oligos <- tibble::tibble(begin, end, length, majority, iupac, majority_rc, iupac_rc, degenerates, degeneracy)
        # Exclude oligos with positions with a gap frequency higher than the stated threshold
        gap_bin <- ifelse(x$gaps > max_gap_frequency, 1L, 0L)
        gap_penalty <- running_sum(gap_bin, n = y)
        oligos <- oligos[which(gap_penalty == 0), ]
        # Exclude oligos with more degenerate bases than the stated threshold
        oligos <- oligos[which(oligos$degenerates <= max_degenerates), ]
        # Exclude oligos with higher degeneracy than the stated threshold
        oligos <- oligos[which(oligos$degeneracy <= max_degeneracy), ]
        # Exclude oligos with '-' (gaps)
        if (any(grepl("-", oligos$majority))) {
            oligos <- oligos[-grep("-", oligos$majority), ]
        }
        # Identify and exclude oligos that are duplicated
        unique_oligos <- match(oligos$majority, unique(oligos$majority))
        oligos <- oligos[unique_oligos, ]
        # Calculate GC content of all majority oligos
        gc_majority <- purrr::map_dbl(oligos$majority, ~gc_content(.x))
        oligos <- tibble::add_column(oligos, gc_majority)
        # Exclude oligos with GC content outside the stated thresholds
        oligos <- oligos[which(gc_majority >= min(gc_range) & gc_majority <= max(gc_range)), ]
        # Calculate Tm of all majority oligos
        tm_majority <- tm(oligos$majority)
        oligos <- tibble::add_column(oligos, tm_majority)
        # Exclude oligos with Tm outside the stated thresholds
        oligos <- oligos[which(tm_majority >= min(tm_range) & tm_majority <= max(tm_range)), ]
        # Indentify non-complex oligos and (if specified) oligos with 3' t-a or 5' g, and set to NA
        oligos$majority <- exclude_unwanted_oligos(oligos$majority, avoid_3end_ta = avoid_3end_ta, avoid_5end_g = avoid_5end_g,
            avoid_3end_runs = avoid_3end_runs)
        oligos$majority_rc <- exclude_unwanted_oligos(oligos$majority_rc, avoid_3end_ta = avoid_3end_ta, avoid_5end_g = avoid_5end_g,
            avoid_3end_runs = avoid_3end_runs)
        oligos$iupac[which(is.na(oligos$majority))] <- NA
        oligos$iupac_rc[which(is.na(oligos$majority_rc))] <- NA
        # Identify oligos where both the sense and antisense sequence is NA
        invalid_oligos <- dplyr::intersect(which(is.na(oligos$majority)), which(is.na(oligos$majority_rc)))
        if (length(invalid_oligos) > 0L)
            oligos <- oligos[-invalid_oligos, ]
        return(oligos)
    })
    if (nrow(all_oligos) == 0L)
        stop("No oligos were found.", call. = FALSE)
    all_oligos <- dplyr::arrange(all_oligos, begin)
    all_oligos <- tibble::new_tibble(all_oligos, nrow = nrow(all_oligos), class = "rprimer_oligo")
    return(all_oligos)
}

# Get PCR assays from oligos ==================================================

#' Get (RT)-PCR assays from oligos
#'
#' \code{get_assays} combines forward and reverse primers to (RT)-PCR assays.
#'
#' @param x an object of class 'rprimer_oligo'.
#'
#' @param length Desired amplicon length,
#' ranging from 40 to 5000 base pairs. The
#' default is \code{65:120}, which means that assays with amplicon lengths from
#' 65 to 120 base-pairs will be considered as acceptable.
#'
#' @param max_tm_difference Maximum Tm difference (in C) between the two primers
#' (absolute value). A number between 0 and 30. The default is 1.
#'
#' @details
#' The Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return A tibble (a data frame) of class 'rprimer_assay' with all candidate
#' assays. An error message will return if no assays are found.
#'
#' @examples
#' get_assays(example_rprimer_oligo, length = 60:150, max_tm_difference = 1)
#'
get_assays <- function(x, length = 65:120, max_tm_difference = 1) {
    if (!inherits(x, "rprimer_oligo")) {
        stop("An rprimer_oligo object is expected.", call. = FALSE)
    }
    if (!(max_tm_difference > 0 && max_tm_difference < 30 && is.numeric(max_tm_difference))) {
        stop("max_tm_difference must be between 0 and 30", call. = FALSE)
    }
    if (!(min(length) >= 40 && max(length) <= 5000 && is.numeric(length))) {
        stop("length must be between 40 and 5000", call. = FALSE)
    }
    # Remove columns with match information (if any)
    if (any(grepl("^pm_", names(x)))) {
        drop <- grep("^pm_|match_matrix", names(x))
        x <- x[, -drop] # # # # # Logical subsetting ??? include match in get ols?
    }
    # Get all fwd and rev primers
    fwd <- x[which(!is.na(x$majority)), ]
    rev <- x[which(!is.na(x$majority_rc)), ]
    # Get all possible combinations of fwd and rev primers
    combinations <- expand.grid(fwd$majority, rev$majority_rc, stringsAsFactors = FALSE)
    names(combinations) <- c("fwd_majority", "rev_majority")
    # Get the indexes of the primer sequences
    index_fwd <- match(combinations$fwd_majority, x$majority)
    index_rev <- match(combinations$rev_majority, x$majority_rc)
    # Expand and combine
    fwd <- x[index_fwd, ]
    rev <- x[index_rev, ]
    # Add a tag on the colnames
    colnames(fwd) <- paste0(colnames(fwd), "_fwd")
    colnames(rev) <- paste0(colnames(rev), "_rev")
    # Combine in a tibble
    assays <- list(fwd, rev)
    assays <- do.call("cbind", assays)
    assays <- tibble::as_tibble(assays)
    # Add amplicon length, tm difference and total degeneracy to the data
    amplicon_length <- assays$end_rev - assays$begin_fwd + 1
    amplicon_length <- as.integer(amplicon_length)
    tm_difference_primer <- abs(assays$tm_majority_fwd - assays$tm_majority_rev)
    begin <- assays$begin_fwd
    end <- assays$end_rev
    total_degeneracy <- assays$degeneracy_fwd + assays$degeneracy_rev
    assays <- tibble::add_column(assays, begin, end, amplicon_length, tm_difference_primer, total_degeneracy, .before = "begin_fwd")
    # Drop columns that we do not longer need
    drop <- c("majority_rc_fwd", "iupac_rc_fwd", "majority_rev", "iupac_rev")
    assays <- assays[, !(names(assays) %in% drop)]
    # Rename columns
    names(assays)[grep("_rc", names(assays))] <- c("majority_rev", "iupac_rev")
    # Collect assays with correct amplicon length
    assays <- assays[which(assays$amplicon_length >= min(length) & assays$amplicon_length <= max(length)), ]
    assays <- assays[which(assays$tm_difference_primer <= max_tm_difference), ]
    if (nrow(assays) == 0L)
        stop("No assays were found.", call. = FALSE)
    assays <- dplyr::arrange(assays, begin)
    assays <- tibble::new_tibble(assays, class = "rprimer_assay")
    return(assays)
}

# Add probes to PCR assays ====================================================

#' Add probes to (RT)-PCR assays
#'
#' \code{add_probes} adds probes to (RT)-PCR assays.
#'
#' @param x Assays to add probes to (an object of class 'rprimer_assay').
#'
#' @param y Candidate probes (an object of class 'rprimer_oligo').
#'
#' @param tm_difference The Tm-difference range between primers and probe.
#' The default values are \code{c(0, 20)}.The minimum allowed
#' value is -20 and the maxiumum allowed value  is 20.
#'
#' @details The Tm-difference is calculated by subtracting the
#' Tm of the probe with the average Tm of the
#' primer pair. Thus, a negative Tm-difference
#' means that the Tm of the probe is lower than the average Tm of the
#' primer pair.
#'
#' #Warning
#' The Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return Assays with probes. A tibble (a data frame) of
#' class 'rprimer_assay'.
#'
#' @examples
#' add_probes(
#' example_rprimer_assay, example_rprimer_oligo, tm_difference = c(-2, 10)
#' )
#'
#' @export
add_probes <- function(x, y, tm_difference = c(0, 20)) {
    if (!inherits(x, "rprimer_assay")) {
        stop("An rprimer_assay object is expected for x.", call. = FALSE)
    }
    if (!inherits(y, "rprimer_oligo")) {
        stop("An rprimer_oligo object is expected for y.", call. = FALSE)
    }
    if (!(min(tm_difference) >= -20 && max(tm_difference) <= 20 && is.numeric(tm_difference))) {
        stop("tm_difference must be between -30 and 30, e.g. c(-1, 5)", call. = FALSE)
    }
    assays <- x
    probes <- y
    # For each assay, take all probes that bind whithin the assay region
    probe_candidates <- purrr::map(seq_len(nrow(assays)), function(i) {
        # The probe has to begin after the fwd primer ends and we want at least one base inbetween (hence the + 1)
        from <- assays$end_fwd[[i]] + 1
        # The probe has to end before the rev primer begins and we want at least one base inbetween
        to <- assays$begin_rev[[i]] - 1
        # Take all probes that begin and end 'within' the assay region
        probe <- probes[which(probes$begin >= from & probes$end <= to), ]
        probe
    })
    # For each assay, we check how many probe candidates we have
    number_of_probes <- purrr::map_int(probe_candidates, nrow)
    # Then, we need to pick all the assays that can harbour a probe.  For assays that has n number of probes, we need to
    # repeat that row n times.
    rows_to_select <- purrr::map(seq_along(number_of_probes), ~rep(.x, number_of_probes[[.x]]))
    rows_to_select <- unlist(rows_to_select, use.names = FALSE)
    assays <- assays[rows_to_select, ]
    if (nrow(assays) == 0L)
        stop("No assays could be generated.", call. = FALSE)
    # Now, we can make a data frame of the probe candidates
    probe_candidates <- do.call("rbind", probe_candidates)
    # Check the sense of the probe, if the plus sense probe is ok, (i.e. not NA), we take that one.
    sense <- ifelse(!is.na(probe_candidates$majority), "pos", "neg")
    probe_candidates$majority <- ifelse(sense == "pos", probe_candidates$majority, probe_candidates$majority_rc)
    probe_candidates$iupac <- ifelse(sense == "pos", probe_candidates$iupac, probe_candidates$iupac_rc)
    probe_candidates <- tibble::add_column(probe_candidates, sense)
    drop <- c("majority_rc", "iupac_rc")
    probe_candidates <- probe_candidates[, -which(names(probe_candidates) %in% drop)]
    names(probe_candidates) <- paste0(names(probe_candidates), "_pr")
    assays <- dplyr::bind_cols(assays, probe_candidates)
    assays$total_degeneracy <- assays$total_degeneracy + probe_candidates$degeneracy_pr
    tm_difference_primer_probe <- purrr::map_dbl(seq_len(nrow(assays)), ~assays$tm_majority_pr[[.x]] - mean(assays$tm_majority_fwd[[.x]],
        assays$tm_majority_rev[[.x]]))
    assays <- tibble::add_column(assays, tm_difference_primer_probe, .after = "tm_difference_primer")
    assays <- assays[which(assays$tm_difference_primer_probe >= min(tm_difference) & assays$tm_difference_primer_probe <=
        max(tm_difference)), ]
    if (nrow(assays) == 0L)
        stop("No assays could be generated.", call. = FALSE)
    assays <- tibble::new_tibble(assays, class = "rprimer_assay")
    return(assays)
}

# Check if oligos and assays matches their targets ============================

#' Check if oligos or assays matches their targets (generic)
#'
#' \code{check_match} checks if oligos or assays matches with their
#' intended target sequences.
#'
#' @param x An object of class
#' 'rprimer_oligo' or 'rprimer_assay'.
#'
#' @param y The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @return A tibble (
#' a data frame) of class 'rprimer_oligo' or 'rprimer_assay',
#' with columns describing the proportion of perfectly matching
#' sequences for each oligo or assay, and a column named
#' 'match_report', which contains a matrix with information about
#' which sequences the oligo/assay matches perfectly to.
#'
#' @examples
#' check_match(example_rprimer_oligo, example_rprimer_alignment)
#' check_match(example_rprimer_assay, example_rprimer_alignment)
#'
#' @export
check_match <- function(x, y) {
    if (!inherits(y, "rprimer_alignment")) {
        stop("An rprimer_alignment object is expected for y.", call. = FALSE)
    }
    UseMethod("check_match")
}

#' @describeIn check_match Check match of an object of class
#' 'rprimer_oligo'.
#' @export
check_match.rprimer_oligo <- function(x, y) {
    # If the positive-sense sequence is NA, then take the negative sense and make a reverse complement Make regular
    # expressions of majority and iupac sequences
    majority <- ifelse(!is.na(x$majority), x$majority, purrr::map_chr(x$majority_rc, reverse_complement))
    majority <- purrr::map_chr(majority, make_regex)
    iupac <- ifelse(!is.na(x$iupac), x$iupac, purrr::map_chr(x$iupac_rc, reverse_complement))
    iupac <- purrr::map_chr(iupac, make_regex)
    # Shorten the sequence names to only accession numbers
    names(y) <- purrr::map_chr(names(y), truncate_name)
    # Make a matrix that describes exactly which sequences the oligo matches perfectly to
    match_matrix <- purrr::map(seq_len(nrow(x)), function(i) {
        match_majority <- grepl(majority[[i]], y)
        match_iupac <- grepl(iupac[[i]], y)
        match <- cbind(match_majority, match_iupac)
        colnames(match) <- c("match_majority", "match_iupac")
        rownames(match) <- names(y)
        return(match)
    })
    names(match_matrix) <- x$iupac
    # Calculate the match percentage for each oligo
    match_percentage <- purrr::map(match_matrix, colMeans)
    match_percentage <- do.call("rbind", match_percentage)
    colnames(match_percentage) <- c("pm_majority", "pm_iupac")
    match_percentage <- tibble::as_tibble(match_percentage)
    # Add this information to the rprimer_oligo object
    x <- dplyr::bind_cols(x, match_percentage)
    match_matrix <- tibble::tibble(match_matrix)
    x <- dplyr::bind_cols(x, match_matrix)
    x <- tibble::new_tibble(x, nrow = nrow(x), class = "rprimer_oligo")
    return(x)
}

#' @describeIn check_match Check match of an object of class
#' 'rprimer_assay'.
#' @export
check_match.rprimer_assay <- function(x, y) {
    # Shorten the names of y to only accession numbers
    names(y) <- purrr::map_chr(names(y), truncate_name)
    # Get match matrices for each assay
    primer_match_matrix <- check_primer_match(x, y)
    if (any(grepl("_pr$", names(x)))) {
        probe_match_matrix <- check_probe_match(x, y)
        match_matrix <- purrr::map(seq_len(nrow(x)), function(x) {
            match <- cbind(primer_match_matrix[[x]], probe_match_matrix[[x]])
            return(match)
        })
    } else match_matrix <- primer_match_matrix
    match_matrix <- purrr::map(match_matrix, function(x) {
        majority_cols <- seq(1, ncol(x), 2)
        iupac_cols <- seq(2, ncol(x), 2)
        majority_all <- rowSums(x[, majority_cols])
        iupac_all <- rowSums(x[, iupac_cols])
        majority_all <- ifelse(majority_all == length(majority_cols), TRUE, FALSE)
        iupac_all <- ifelse(iupac_all == length(iupac_cols), TRUE, FALSE)
        x <- cbind(x, majority_all, iupac_all)
        return(x)
    })
    names(match_matrix) <- paste0(x$majority_fwd, "_", x$majority_rev, "_", x$majority_pr)
    match_percentage <- purrr::map(match_matrix, colMeans)
    match_percentage <- do.call("rbind", match_percentage)
    match_percentage <- tibble::as_tibble(match_percentage)
    names(match_percentage) <- paste0("pm_", names(match_percentage))
    if (any(grepl("^pm_", names(x)))) {
        stop("matches have already been checked in x", call. = FALSE)
    }
    x <- dplyr::bind_cols(x, match_percentage)
    match_matrix <- tibble::tibble(match_matrix)
    x <- dplyr::bind_cols(x, match_matrix)
    x <- tibble::new_tibble(x, nrow = nrow(x), class = "rprimer_assay")
    return(x)
}

# Plot results ================================================================

#' Plot an rprimer-object (generic)
#'
#' @param x An rprimer-object (see methods for details).
#'
#' @return A plot.
#'
#' @examples
#' rp_plot(example_rprimer_alignment)
#' rp_plot(example_rprimer_sequence_profile, from = 1, to = 20, rc = TRUE)
#' rp_plot(example_rprimer_sequence_properties)
#'
#' @export
rp_plot <- function(x, ...) {
    object_name <- as.character(substitute(x))
    if (!exists(object_name)) {
        stop(paste("object", object_name, "does not exist."), call. = FALSE)
    }
    UseMethod("rp_plot")
}

#' @describeIn rp_plot Plot an object of class 'rprimer_alignment'.
#' @export
rp_plot.rprimer_alignment <- function(x, ...) {
    op <- par(mar = c(5, 3, 3, 3))
    on.exit(par(op))
    # Make a matrix of the alignment
    index <- seq_along(x)
    sequences_as_integers <- purrr::map(index, function(i) {
        sequence <- split_sequence(x[[i]])
        sequence <- gsub("-", NA, sequence)
        sequence <- as.integer(gsub("[a-z]", i, sequence))
        return(sequence)
    })
    to_plot <- do.call("rbind", sequences_as_integers)
    plot(seq_len(ncol(to_plot)), rep(1, ncol(to_plot)), ylab = "", xlab = "position in consensus sequence", yaxt = "n",
        pch = NA, ylim = c(0, max(nrow(to_plot)) + 1))
    invisible(apply(to_plot, 1, function(x) lines(seq_along(x), x, col = "grey")))
}

#' @describeIn rp_plot Plot an object of class 'rprimer_sequence_profile'.
#'
#' @param from At which position the plot begins (an integer).
#'
#' @param to At which position the plot ends (an integer).
#'
#' @param rc \code{TRUE/FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. The default is \code{FALSE}.
#'
#' @export
rp_plot.rprimer_sequence_profile <- function(x, ..., from = NULL, to = NULL, rc = FALSE) {
    if (is.null(from))
        from <- 1
    if (is.null(to))
        to <- ncol(x)
    if (!(is.numeric(from) && is.numeric(to) && length(from) == 1 && length(to) == 1)) {
        stop("from and to must be numbers", call. = FALSE)
    }
    if (!(from <= ncol(x) && to <= ncol(x))) {
        stop("from and to cannot be greater than the number of columns in x", call. = FALSE)
    }
    if (!(from >= 1 && to >= 1)) {
        stop("from and to must be 1 or greater", call. = FALSE)
    }
    if (!(is.logical(rc)))
        stop("rc must be set to TRUE or FALSE", call. = FALSE)
    op <- par(mar = c(0.75, 4.57, 0.75, 0.75))
    on.exit(par(op))
    # Get the data
    selection <- x[, from:to]
    selection <- selection[which(rownames(selection) != "-"), ]
    if (rc == TRUE) {
        selection <- selection[, ncol(selection):1]  # rev seq len instead
        rownames(selection) <- unname(complement_lookup[rownames(selection)])
    }
    sequence_barplot(selection)
}

#' @describeIn rp_plot Plot an object of class 'rprimer_sequence_properties'.
#' @export
rp_plot.rprimer_sequence_properties <- function(x, ...) {
    op <- par(mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1), xpd = FALSE, oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1))
    on.exit(par(op))
    identity_plot <- plot(x$position, x$identity, type = "h", ylim = c(0, 1), ylab = "identity", xlab = "", xaxt = "n",
        col = ifelse(x$identity < 1, "gray80", "gray60"))
    identity_line <- lines(running_average(x$identity))
    entropy_plot <- plot(x$position, x$entropy, type = "h", ylim = c(0, max(x$entropy, na.rm = TRUE) * 1.1), ylab = "shannon entropy",
        xlab = "", xaxt = "n", col = ifelse(x$entropy > 0, "gray80", "gray60"))
    entropy_line <- lines(running_average(x$entropy))
    gc_plot <- plot(x$position, pch = NA, ylab = "gc content", ylim = c(0, 1), xlab = "", xaxt = "n")
    clip(0, nrow(x), -1, 2)
    abline(h = 0.5, col = "gray80")
    gc_line <- lines(gc_running_average(x$majority))
    gap_plot <- plot(x$position, x$gaps, type = "h", ylim = c(0, 1), ylab = "gaps", xlab = "")
    identity_plot
    identity_line
    entropy_plot
    entropy_line
    gc_plot
    gc_line
    gap_plot
    mtext(side = 1, outer = TRUE, "position in consensus sequence", line = 3, cex = 0.7)
}

# Save results ================================================================

#' Save an rprimer object to a file (generic)
#'
#' @param x An rprimer object (see methods for details).
#'
#' @param filename The name of the file. A character vector of length one.
#'
#' @details The file format will depend on object.
#' Objects of class rprimer_alignment will be saved as .txt-files,
#' whereas objects of class rprimer_sequence_properties, rprimer_oligo
#' and rprimer_assay will be saved as .csv-files.
#'
#' @return A saved file.
#'
#' @examples
#' \dontrun{
#' rp_save(example_rprimer_alignment, 'my_alignment')
#' rp_save(example_rprimer_sequence_properties, 'my_sequence_properties')
#' rp_save(example_rprimer_oligo, 'my_oligos')
#' rp_save(example_rprimer_assay, 'my_assays)
#' }
#'
#' @export
rp_save <- function(x, filename) {
    if (!is.character(filename) || length(filename) != 1) {
        stop("filename must be a character vector of length one", call. = FALSE)
    }
    object_name <- as.character(substitute(x))
    if (!exists(object_name)) {
        stop(paste("object", object_name, "does not exist."), call. = FALSE)
    }
    UseMethod("rp_save")
}

#' @describeIn rp_save Save an object of class 'rprimer_alignment'.
#' @export
rp_save.rprimer_alignment <- function(x, filename) {
    seq_names <- purrr::map_chr(names(x), ~paste0(">", .x))
    my_file <- file(paste0(filename, ".txt"), open = "w")
    purrr::walk2(seq_names, x, ~writeLines(text = c(.x, .y), con = my_file))
    close(my_file)
}

#' @describeIn rp_save Save an object of class 'rprimer_sequence_properties'.
#' @export
rp_save.rprimer_sequence_properties <- function(x, filename) {
    write.csv(x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE)
}

#' @describeIn rp_save Save an object of class 'rprimer_oligo'.
#' @export
rp_save.rprimer_oligo <- function(x, filename) {
    if (any(grepl("match_matrix", names(x)))) {
        x <- x[, which(names(x) != "match_matrix")]  # Can't save a list
    }
    write.csv(x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE)
}

#' @describeIn rp_save Save an object of class 'rprimer_assay'.
#' @export
rp_save.rprimer_assay <- function(x, filename) {
    if (any(grepl("match_matrix", names(x)))) {
        x <- x[, which(names(x) != "match_matrix")]  # Can't save a list
    }
    write.csv(x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE)
}

# Generate report ============================================================

#' Write assay report
#'
#' @export
write_assay_report <- function(filename = "my_assay_report", assay_selection, alignment, sequence_profile, sequence_properties) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("rmkardown is needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("kableExtra", quietly = TRUE)) {
        stop("kableExtra is needed for this function to work. Please install it.", call. = FALSE)
    }
    filename <- paste0(filename, ".html")
    # if if if inputs
    rmarkdown::render(input = "assay_report.Rmd", output_file = filename, params = list(assay_selection = assay_selection,
        alignment = alignment, sequence_profile = sequence_profile, sequence_properties = sequence_properties))
}
