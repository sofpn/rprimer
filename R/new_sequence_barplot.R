newbarp <- function(x, rc = FALSE) {
  if (!(is.logical(rc))) {
    stop("'rc' must be set to TRUE or FALSE", call. = FALSE)
  }
  x <- assay(x)
  x <- x[rownames(x) != "-", ]
  if (rc == TRUE) {
    x <- x[, rev(seq_len(ncol(x)))]
    rownames(x) <- unname(complementLookup[rownames(x)])
  }
  bases <- c("A", "C", "G", "T")
  Other <- colSums(x[!rownames(x) %in% bases, ])
  x <- x[rownames(x) %in% bases, ]
  x <- x[sort(rownames(x)), ]
  x <- rbind(x, Other)
  x <- reshape2::melt(x)
  x[, 1] <- factor(x[, 1])
  x[, 2] <- factor(x[, 2])
  names(x) <- c("Base", "Position", "Frequency")
  cols <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038", "gray")
  p <- ggplot2::ggplot(x, aes(x = Position, y = Frequency, fill = Base))
  barchart <- ggplot2::geom_bar(stat = "identity")
  coloring <- ggplot2::scale_fill_manual(values = cols)
  appearance <- ggplot2::theme_light()
  options <- ggplot2::theme(legend.title = element_blank())
  #ordering <- ggplot2::scale_x_discrete(limits= x$Position)
  p + barchart + coloring + appearance + options #+ ordering
}

