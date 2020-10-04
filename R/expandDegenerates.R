# to to: write tests
# x - a DNA sequence (e.g. "cggttrt")
expandDegenerates <- function(x) {
  x <- splitSequence(x)
  expanded <- purrr::map(x, function(i) {
    allBases <- unname(degenerateLookup[[i]])
    allBases <- unlist(strsplit(allBases, split = ","))
    allBases
  })
  expanded <- expand.grid(
    expanded[seq_along(expanded)], stringsAsFactors = FALSE
  )
  expanded <- purrr::map(
    seq_len(nrow(expanded)), ~paste(expanded[.x, ], collapse = "")
  )
  expanded <- unlist(expanded, use.names = FALSE)
  expanded
}
