#' Plot an rprimer object (generic)
#'
#' @param x An rprimer object (see methods for details).
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
#'
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
    ylab = "", xlab = "position", yaxt = "n",
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
#'
#' @export
rp_plot.rprimer_properties <- function(x, ...) {
  op <- graphics::par(
    mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1),
    oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1)
  )
  on.exit(graphics::par(op))
  sequence_detail_plot(x)
}

#' Calculate running average
#'
#' \code{running_average} calculates the centered running average of a
#' numeric vector.
#'
#' @param x A numeric vector.
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c(1, 3, 2, 1, 3, 4, 5), size = 2)
#' running_average(c(10, 22.4, 14.2, 44, 32))
#'
#' @keywords internal
#'
#' @noRd
running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.numeric(x) && size <= length(x) && size >= 1)
  sums <- c(0, cumsum(x))
  average <- (sums[((size + 1):length(sums))] - sums[(1:(length(sums) - size))]) / size
  # Set the position to get a trailed moving average
  position <- seq(size, length(x))
  # To get a centered running average,
  # I align each moving average to the midpoint of the range of observations
  # that each average includes
  if (size > 1) {
    midpoint <- size / 2
    position <- position - midpoint
  }
  df <- tibble::tibble(position, average)
  df
}

#' Calculate running average of GC content
#'
#' \code{gc_running_average} calculates the centered running average of
#' the GC-content of a DNA sequence
#'
#' @param x a DNA sequence (a character vector of length one).
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c("gtttgtttctt"), size = 2)
#' running_average(c("cttggggtttctttggtt-ttagg"))
#'
#' @seealso
#' gc_content
#'
#' @keywords internal
#'
#' @noRd
gc_running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.character(x) && size <= length(x) && size >= 1)
  begins <- seq_len(length(x) - size + 1)
  ends <- begins + size - 1
  # Calculate gc content of chunks of x
  average <- purrr::map2_dbl(begins, ends, function(y, z) {
    string <- paste(x[y:z], collapse = "")
    gc_content(string)
  })
  # Set the position to get a trailed moving average
  position <- seq(size, length(x))
  # To get a centered running average,
  # I align each moving average to the midpoint of the range of observations
  # that each average includes
  if (size > 1) {
    midpoint <- size / 2
    position <- position - midpoint
  }
  df <- tibble::tibble(position, average)
  df
}

#' Draw a rectangle
#'
#' @param from Where (at the x-axis) the rectangle starts.
#'
#' @param to Where (at the x-axis) the rectangle ends.
#'
#' @return A rectangle.
#'
#' @keywords internal
#'
#' @noRd
rectangle <- function(from, to) {
  graphics::rect(
    from, 0, to, 5,
    border = NA,
    col = grDevices::rgb(
      123, 149, 169, alpha = 100, maxColorValue = 200
    ), xpd = NA
  )
}

#' Sequence detail plot
#'
#' @param x An object of class 'rprimer_properties'.
#'
#' @return A visual representation of \code{x}.
#'
#' @keywords internal
#'
#' @noRd
sequence_detail_plot <- function(x) {
  graphics::plot(
    x$position, x$identity,
    type = "h", ylim = c(0, 1),
    ylab = "identity", xlab = "", xaxt = "n",
    col = ifelse(x$identity < 1, "gray80", "gray60")
  )
  graphics::lines(running_average(x$identity))
  graphics::plot(
    x$position, x$entropy,
    type = "h",
    ylim = c(0, max(x$entropy, na.rm = TRUE) * 1.1),
    ylab = "shannon entropy",
    xlab = "", xaxt = "n", col = ifelse(x$entropy > 0, "gray80", "gray60")
  )
  graphics::lines(running_average(x$entropy))
  graphics::plot(
    x$position, rep(0.5, length(x$position)),
    pch = NA, ylab = "gc content", ylim = c(0, 1),
    xlab = "", xaxt = "n"
  )
  graphics::clip(0, max(x$position, na.rm = TRUE), 0, 3)
  graphics::abline(h = 0.5, col = "gray80")
  gc_line <- graphics::lines(gc_running_average(x$majority))
  graphics::plot(
    x$position, x$gaps,
    type = "h", ylim = c(0, 1), ylab = "gaps", xlab = ""
  )
  graphics::mtext(
    side = 1, outer = TRUE, "position",
    line = 3, cex = 0.7
  )
}

#' Sequence barplot
#'
#' @param x A selection of a sequence profile
#'
#' @param ... Additional arguments that should be
#' passed to the barplot, e.g. a title.
#'
#' @return A barplot.
#'
#' @keywords internal
#'
#' @noRd
sequence_barplot <- function(x, ...) {
  base_cols <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038", "white")
  names(base_cols) <- c("a", "c", "g", "t", "-")
  gray_cols <- grDevices::gray.colors(nrow(x) - 5, start = 0.6) ######## ALL POSSIBLE WBBL BASES HERE
  names(gray_cols) <- setdiff(rownames(x), c("a", "c", "g", "t", "-"))
  cols <- c(base_cols, gray_cols)
  # Make the plot
  graphics::barplot(
    x,
    space = 0, xaxt = "n", font.main = 1, border = "grey80",
    col = cols[rownames(x)], legend = TRUE, ylab = "Proportion",
    ylim = c(0, 1), ...,
    args.legend = list(
      x = "right", box.col = NA, border = "grey80",
      bg = grDevices::rgb(
        60, 60, 60, alpha = 50, maxColorValue = 200
      ), inset = c(0, -0.3)
    )
  )
}
