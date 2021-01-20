## Import all genotyping results
results <- read.table("typing_results.txt", sep = "\t", header = TRUE)

## Identify all Ortohepevirus A sequences
valid <- results$BLAST.result == "Hepeviridae Orthohepevirus A"
hev <- results$name[valid]

## Randomly select 200 sequences to include as example data in the package
set.seed(1)
selection <- sample(hev, 200)

## Save accession numbers for the selected sequences so that they can be
## downloaded using batch entrez
write.table(
    selection, file = "accession_valid_sequences.txt", quote = FALSE,
    row.names = FALSE, col.names = FALSE, sep = "\t"
)
