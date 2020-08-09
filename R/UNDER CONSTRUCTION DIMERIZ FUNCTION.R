library(tidyverse)

# Functions
split_seq <- function(x) {
  x <- strsplit(x, split = "")
  x <- unlist(x, use.names = FALSE)
  return(x)
}

# Enter the sequences of the primer pair to check
x <- "tata"
x <- complement(x)
x <- split_seq(x) # Longest
y <- split_seq("tata") # Shortest

# Make a matrix where the shortest primer pair "slides" along the longest primer pair
# to get all possible base-pairing combinations
M <- matrix(nrow = length(x) + length(y) - 1, ncol = length(x) + 2*(length(y) - 1))
M <- map(seq_len(nrow(M)), function(i) {
  row_val <- c(rep(".", i - 1), y)
  row_val <- c(row_val, rep(".", ncol(M) - length(row_val)))
  M[i, ] <- row_val
})
M <- do.call("rbind", M)
# Add the complement sequence of the plus strand as first row
M <- rbind(c(rep("-", length(y) - 1), x, rep("-", length(y) - 1)), M)

# Split the primers into nearest neigbours
N <- matrix(nrow = nrow(M), ncol = ncol(M) - 1)
N[] <- paste(M[, 1:(ncol(M) - 1)], M[, 2:ncol(M)], sep = "")
N
# Match each row with lookup .....

matching <- apply(M, 2, function(x) ifelse(x[[1]] == x[-1], 1, 0))
M
matching
