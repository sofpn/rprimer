#' @noRd
print_assay_report <- function(x) { # x rprimer assay obj - one row
  if (any(grepl("match_matrix", names(x)))) {
    return(
      list(
        print_assay(x),
        print_degenerates(x),
        print_match_matrix(x$match_matrix[[1]])
      )
    )
  } else {
    return(list(
      print_assay(x),
      print_degenerates(x)
    )
    )
  }
}

assay_detail_plot <- function(x, y) { # x = rprimer sequence profile y = assay object (one row)
  x <- x[which(rownames(x) != "-"), ]
  fwd <- x[, y$begin_fwd:y$end_fwd]
  rev <- x[, y$begin_rev:y$end_rev]
  rev <- rev[, ncol(rev):1]
  rownames(rev) <- unname(
    complement_lookup[rownames(rev)]
  )
  if (any(grepl("_pr", names(y)))) {
    pr <- x[, y$begin_pr:y$end_pr]
    if (y$sense_pr == "neg") {
      pr <- pr[, ncol(pr):1]
      rownames(pr) <- unname(
        complement_lookup[rownames(pr)]
      )
    }
  }
  if (any(grepl("_pr", names(y)))) {
    op <- par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow =  c(1, 3))
    on.exit(par(op))
    sequence_barplot(fwd, main = "forward")
    sequence_barplot(rev, main = "reverse")
    sequence_barplot(pr, main = "probe")
  } else {
    op <- par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow =  c(1, 2))
    on.exit(par(op))
    sequence_barplot(fwd, main = "forward")
    sequence_barplot(rev, main = "reverse")
  }
}

#' @noRd
assay_overview_plot <- function(x, y) { # x rprimer_sequence properties # y assay object (one row)
  op <- par(
    mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1), xpd = FALSE,
    oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1)
  )
  on.exit(par(op))
  identity_plot <- plot(
    x$position, x$identity, type = "h",
    ylim = c(0, 1), ylab = "identity",
    xlab = "", xaxt = "n",
    col = ifelse(x$identity < 1, "gray80", "gray60")
  )
  identity_line <- lines(running_average(x$identity))
  entropy_plot <- plot(
    x$position, x$entropy, type = "h",
    ylim = c(0, max(x$entropy, na.rm = TRUE)*1.1),
    ylab = "shannon entropy", xlab = "", xaxt = "n",
    col = ifelse(x$entropy > 0, "gray80", "gray60")
  )
  entropy_line <- lines(running_average(x$entropy))
  gc_plot <- plot(
    x$position, pch = NA, ylab = "gc content",
    ylim = c(0, 1), xlab = "", xaxt = "n"
  )
  clip(0, nrow(x), -1, 2)
  abline(h = 0.5, col = "gray80")
  gc_line <- lines(gc_running_average(x$majority))
  gap_plot <- plot(
    x$position, x$gaps, type = "h", ylim = c(0, 1),
    ylab = "gaps", xlab = ""
  )
  identity_plot
  identity_line
  entropy_plot
  entropy_line
  gap_plot
  rectangle(y$begin_fwd, y$end_fwd)
  rectangle(y$begin_rev, y$end_rev)
  if (any(grepl("_pr$", names(x)))) rectangle(y$begin_pr, y$end_pr)
  mtext(
    side = 1, outer = TRUE, "position in consensus sequence",
    line = 3, cex = 0.7
  )
}

#' @noRd
print_assay <- function(x) { # x assay object (one row)
  if (any(grepl("match_matrix", names(x)))) {
    x <- x[which(!grepl("match_matrix", names(x)))]
  }
  all <- x[which(!grepl("_fwd$|_rev$|_pr$", names(x)))]
  names(all) <- gsub("_all$", "", names(all))
  fwd <- x[which(grepl("_fwd$", names(x)))]
  names(fwd) <- gsub("_fwd$", "", names(fwd))
  rev <- x[which(grepl("_rev$", names(x)))]
  names(rev) <- gsub("_rev$", "", names(rev))
  if (any(grepl("_pr", names(x)))) {
    pr <- x[which(grepl("_pr$", names(x)))]
    names(pr) <- gsub("_pr$", "", names(pr))
    my_list <- list(
      "general" = all,
      "forward" = fwd,
      "reverse" = rev,
      "probe" = pr
    )
  } else {
    my_list <- list(
      "general" = all,
      "forward" = fwd,
      "reverse" = rev
    )
  }
  return(my_list)
}

#' @noRd
print_degenerates <- function(x) { # x assay object (one row)
  iupac_oligos <- x[grep("^iupac", names(x))]
  names(iupac_oligos) <- gsub("iupac_", "", names(iupac_oligos))
  to_print <- purrr::map(iupac_oligos, function(x) {
    variant <- expand_degenerate(x)
    gc_content <- purrr::map_dbl(variant, ~ gc_content(.x))
    tm <- tm(variant)
    table <- tibble::tibble(variant, gc_content, tm)
    mean_values <- round(
      c("gc_content" = mean(gc_content), "tm" = mean(tm)), 2
    )
    return(list("all_variants" = table, "mean_values"= mean_values))
  })
  return(to_print)
}

#' @noRd
print_match_matrix <- function(x) { # assay obj (one row)
  iupac <- x[, grep("iupac", colnames(x))]
  majority <- x[, grep("majority", colnames(x))]
  pm_majority <-  majority[ , which(colnames(majority) == "majority_all")]
  pm_majority <- names(pm_majority)[which(pm_majority == TRUE)]
  if (length(pm_majority) != 0) {
    to_remove_majority <- match(pm_majority, rownames(majority))
    mm_majority <- majority[-to_remove_majority, -ncol(majority)]
  } else {
    mm_majority <- majority
  }
  colnames(mm_majority) <- gsub("iupac", "", colnames(mm_majority))
  pm_iupac <-  iupac[ , which(colnames(iupac) == "iupac_all")]
  pm_iupac <- names(pm_iupac)[which(pm_iupac == TRUE)]
  if (length(pm_majority) != 0) {
    to_remove_iupac <- match(pm_iupac, rownames(iupac))
    mm_iupac <- iupac[-to_remove_iupac, -ncol(iupac)]
  } else {
    mm_iupac <- iupac
  }
  colnames(mm_iupac) <- gsub("majority", "", colnames(mm_iupac))
  result <- list(
    "matches_perfectly_majority" = pm_majority,
    "mismatches_majority" = mm_majority,
    "matches_perfectly_iupac" = pm_iupac,
    "mismatches_iupac" = mm_iupac)
  return(result)
}
