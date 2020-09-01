# Custom functions needed to generate the assay report

# x - a DNA sequence (e.g. "cggttrt")
find_g_repeats <- function(x) {
  g_repeats <- grepl("(g)\\1\\1\\1", x)
  if (g_repeats) "Yes"
  else "No"
}

# x - a DNA sequence (e.g. "cggttrt")
expand_degenerate <- function(x) {
  x <- split_sequence(x)
  # Go through each base of the DNA sequence
  expanded <- purrr::map(x, function(i) {
    # Check which bases the IUPAC base at position i correspond to
    all_bases <- unname(degenerate_lookup[[i]])
    all_bases <- unlist(strsplit(all_bases, split = ","))
    return(all_bases)
  })
  # Get all possible combinations of DNA sequences
  expanded <- expand.grid(
    expanded[seq_along(expanded)], stringsAsFactors = FALSE
  )
  expanded <- purrr::map(
    seq_len(nrow(expanded)), ~paste(expanded[.x, ], collapse = "")
  )
  expanded <- unlist(expanded, use.names = FALSE)
  return(expanded)
}

# x - an object of class rprimer_assay (one row)
print_assay_report <- function(x) {
  x <- x[which(!grepl("match_matrix", names(x)))]
  all <- x[which(!grepl("_fwd$|_rev$|_pr$", names(x)))]
  names(all) <- gsub("_all$", "", names(all))
  fwd <- x[which(grepl("_fwd$", names(x)))]
  names(fwd) <- gsub("_fwd$", "", names(fwd))
  rev <- x[which(grepl("_rev$", names(x)))]
  names(rev) <- gsub("_rev$", "", names(rev))
  if (any(grepl("_pr$", names(x)))) {
    pr <- x[which(grepl("_pr$", names(x)))]
    names(pr) <- gsub("_pr$", "", names(pr))
    assay <- list(general = all, forward = fwd, reverse = rev, probe = pr)
  } else {
    assay <- list(general = all, forward = fwd, reverse = rev)
  }
  iupac_oligos <- x[grep("^iupac", names(x))]
  names(iupac_oligos) <- gsub("iupac_", "", names(iupac_oligos))
  degenerates <- purrr::map(iupac_oligos, function(x) {
    variant <- expand_degenerate(x)
    gc_content <- purrr::map_dbl(variant, ~gc_content(.x))
    tm <- tm(variant)
    table <- tibble::tibble(variant, gc_content, tm)
    mean_values <- round(c(gc_content = mean(gc_content), tm = mean(tm)), 2)
    return(list(all_variants = table, mean_values = mean_values))
  })
  return(list(assay, degenerates))
}

# x - an object of class rprimer_sequence_properties
# y - an object of class rprimer_assay (one row)
plot_assay_overview <- function(x, y) {
  op <- graphics::par(
    mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1),
    xpd = FALSE, oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1)
  )
  on.exit(graphics::par(op))
  sequence_detail_plot(x)
  rectangle(y$begin_fwd, y$end_fwd)
  rectangle(y$begin_rev, y$end_rev)
  if (any(grepl("_pr$", names(y)))) {
    rectangle(y$begin_pr, y$end_pr)
  }
}

# x - an object of class rprimer_sequence_profile
# y - an object of class rprimer_assay (one row)
plot_assay_details <- function(x, y) {
  x <- x[which(rownames(x) != "-"), ]
  fwd <- x[, y$begin_fwd:y$end_fwd]
  rev <- x[, y$begin_rev:y$end_rev]
  rev <- rev[, ncol(rev):1]
  rownames(rev) <- unname(complement_lookup[rownames(rev)])
  if (any(grepl("_pr$", names(y)))) {
    pr <- x[, y$begin_pr:y$end_pr]
    if (y$sense_pr == "neg") {
      pr <- pr[, ncol(pr):1]
      rownames(pr) <- unname(complement_lookup[rownames(pr)])
    }
  }
  if (any(grepl("_pr$", names(y)))) {
    op <- graphics::par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow = c(1, 3))
    on.exit(graphics::par(op))
    sequence_barplot(fwd, main = "forward")
    sequence_barplot(rev, main = "reverse")
    sequence_barplot(pr, main = "probe")
  } else {
    op <- graphics::par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow = c(1, 2))
    on.exit(graphics::par(op))
    sequence_barplot(fwd, main = "forward")
    sequence_barplot(rev, main = "reverse")
  }
}

# x - a matrix
split_matrix <- function(x) {
  n <- nrow(x)
  out <- split.data.frame(x, rep(seq_len(ceiling(n / 100)), each = 50)[1:n])
  return(out)
}

# x - an object of class rprimer_assay (one row)
plot_match_matrix <- function(assay_selection) {
  x <- assay_selection$match_matrix[[1]]
  # Convert from logical to integer
  x[] <- as.integer(x)
  # Remove columns with "all"
  x <- x[ , !grepl("_all", colnames(x))]
  colnames(x) <- gsub("match_", "", colnames(x))
  # Split the matrix to several matrices to make neater plots
  matrices <- split_matrix(x)
  purrr::walk(matrices, function(x) {
    # Make subsets, one with majority and one with IUPAC
    majority <- x[, grepl("majority", colnames(x))]
    colnames(majority) <- gsub("majority_", "", colnames(majority))
    iupac <- x[, grepl("iupac", colnames(x))]
    colnames(iupac) <- gsub("iupac_", "", colnames(iupac))
    # Set plot colors
    colors <- c(mismatch = "#E07A5F", match = "#81B29A")
    # Set plot layout
    graphics::layout(
      matrix(data = c(1, 2, 3), nrow = 1, ncol = 3),
      widths = c(0.6, 0.6, 0.2), heights = c(1, 1)
    )
    # Plot majority
    op <- graphics::par(mar = c(3, 5, 3, 3))
    on.exit(graphics::par(op))
    graphics::image(
      seq_along(colnames(majority)), seq_along(rownames(majority)), t(majority),
      col = colors, xlab = "", ylab = "", axes = FALSE, main = "Majority"
    )
    graphics::axis(
      side = 1, at = seq_along(colnames(majority)),
      labels = colnames(majority), cex.axis = 1, tick = FALSE
    )
    graphics::axis(
      side = 2, at = seq_along(rownames(x)), labels = rownames(x), las = 1,
      cex.axis = 0.7, tick = FALSE
    )
    # Plot IUPAC
    graphics::image(
      seq_along(colnames(iupac)), seq_along(rownames(iupac)), t(iupac),
      col = colors, xlab = "", ylab = "", axes = FALSE, main = "IUPAC"
    )
    graphics::axis(
      side = 1, at = seq_along(colnames(iupac)),
      labels = colnames(iupac), cex.axis = 1, tick = FALSE
    )
    op <- graphics::par(mar = c(3, 1, 3, 3))
    on.exit(graphics::par(op))
    # Add a legend
    graphics::image(matrix(-4:5, ncol = 10, nrow = 1),
                    col = c(rep("white", 4), colors, rep("white", 4)),
                    axes = FALSE, xlab = "", ylab = ""
    )
    graphics::text(0, 0.45, "mismatch", cex = 1)
    graphics::text(0, 0.55, "match", cex = 1)
  })
}
