#' Plot an rprimer-object (generic)
#'
#' @param x An rprimer-object (see methods for details).
#'
#' @param ...
#' Additional arguments that should be passed to the plot, when
#' plotting an object of class 'rprimer_profile'.
#'
#' @return A plot.
#'
#' @examples
#' rp_plot(example_rprimer_alignment)
#' rp_plot(example_rprimer_profile, from = 1, to = 20, rc = TRUE)
#' rp_plot(example_rprimer_properties)
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
  op <- graphics::par(mar = c(5, 3, 3, 3))
  on.exit(graphics::par(op))
  # Make a matrix of the alignment
  index <- seq_along(x)
  sequences_as_integers <- purrr::map(index, function(i) {
    sequence <- split_sequence(x[[i]])
    sequence <- gsub("-", NA, sequence)
    sequence <- as.integer(gsub("[a-z]", i, sequence))
    return(sequence)
  })
  to_plot <- do.call("rbind", sequences_as_integers)
  graphics::plot(
    seq_len(ncol(to_plot)), rep(1, ncol(to_plot)),
    ylab = "", xlab = "position in consensus sequence", yaxt = "n",
    pch = NA, ylim = c(0, max(nrow(to_plot)) + 1)
  )
  invisible(
    apply(to_plot, 1, function(x) {
      graphics::lines(seq_along(x), x, col = "grey")
    })
  )
}

#' @describeIn rp_plot Plot an object of class 'rprimer_profile'.
#'
#' @param from At which position the plot begins (an integer).
#'
#' @param to At which position the plot ends (an integer).
#'
#' @param rc \code{TRUE/FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. The default is \code{FALSE}.
#'
#' @export
rp_plot.rprimer_profile <- function(
  x, ..., from = NULL, to = NULL, rc = FALSE
) {
  if (is.null(from))
    from <- 1
  if (is.null(to))
    to <- ncol(x)
  if (!(is.numeric(from) && is.numeric(to))) {
    stop("from and to must be numbers", call. = FALSE)
  }
  if (!(from <= ncol(x) && to <= ncol(x))) {
    stop(
      "from and to cannot be greater than the number of columns in x",
      call. = FALSE
    )
  }
  if (!(from >= 1 && to >= 1)) {
    stop("from and to must be 1 or greater", call. = FALSE)
  }
  if (!(is.logical(rc)))
    stop("rc must be set to TRUE or FALSE", call. = FALSE)
  op <- graphics::par(mar = c(0.75, 4.57, 0.75, 0.75))
  on.exit(graphics::par(op))
  # Get the data
  selection <- x[, from:to]
  selection <- selection[which(rownames(selection) != "-"), ]
  if (rc == TRUE) {
    selection <- selection[, rev(seq_len(ncol(selection)))]
    rownames(selection) <- unname(complement_lookup[rownames(selection)])
  }
  sequence_barplot(selection)
}

#' @describeIn rp_plot Plot an object of class 'rprimer_properties'.
#' @export
rp_plot.rprimer_properties <- function(x, ...) {
  op <- graphics::par(
    mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1),
    xpd = FALSE, oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1)
  )
  on.exit(graphics::par(op))
  sequence_detail_plot(x)
}
