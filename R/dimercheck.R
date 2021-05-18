.slide <- function(upper, lower) {
    upper <- unlist(strsplit(upper, split = ""))
    lower <- unlist(strsplit(lower, split = ""))
    grid <- length(lower) - 1
    firstRow <- c(
        rep("", length(lower) - 1), upper, rep("", length(lower) - 1)
    )
    frame <- seq_len(length(upper) + length(lower) - 1) - 1
    toAdd <- lapply(frame, function(i) {
        j <- i + 1
        c(rep("", grid + length(upper) - j), lower, rep("", i))
    })
    toAdd <- do.call("rbind", toAdd)
    rbind(unname(firstRow), toAdd)
}

.extract <- function(x) {
    upper <- x[1, ]
    lower <- x[-1, ]
    binds <- t(apply(lower, 1, function(y) upper == y & y != "")) ### Match to table degen - if matches at least one....
    nucleotide <- t(apply(lower, 1, function(y) upper != ""))
    out <- ifelse(nucleotide, "*", "")
    out <- ifelse(binds, x, out)
    out
}

.collapse <- function(x) {
    x <- apply(x, 1, paste, collapse = "")
    extract <- grepl("[A-Z]{3,}", x)
    x[extract]
}

.longest <- function(x) {
    x <- strsplit(x, split = "")
    match <- lapply(x, function(y) {
        out <- y != "*"
        names(out) <- seq_along(y)
        as.integer(names(out)[out])
    })
}

x <- .slide("ACTTTTAATTATATTT", "TTTAAGCTTC")
y <- .extract(x)
z <- .collapse(y)
x <- z

