.histPlot <- function(data, y, title = "", color = "#7B95A9") {
  ggplot2::ggplot(data, ggplot2::aes(y)) +
    ggplot2::geom_histogram(
      fill = color, color = color, alpha = 0.4, binwidth = 1
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("Count") +
    ggplot2::labs(title = title) +
    .themeRprimer(showXAxis = TRUE)
}
