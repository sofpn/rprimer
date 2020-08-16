# Code to generate example data

infile <- "example_alignment_100seq.txt"

example_rprimer_alignment <- read_fasta_alignment(infile) %>%
 remove_gaps(., threshold = 0.5)

example_rprimer_profile <- sequence_profile(example_rprimer_alignment)

example_rprimer_properties <- sequence_properties(
  example_rprimer_profile, iupac_threshold = 0.1
)

example_rprimer_oligo <- get_oligos(
  example_rprimer_properties,
  example_rprimer_alignment,
  max_gap_frequency = 0.1,
  length = 18:22,
  max_degenerates = 2,
  max_degeneracy = 4,
  gc_range = c(0.45, 0.56),
  tm_range = c(48, 65),
  avoid_3end_ta = TRUE
)

example_rprimer_assay <- get_assays(
  example_rprimer_oligo, length = 60:150, max_tm_difference = 1.5
 )

save(example_rprimer_alignment, file = "example_rprimer_alignment.RData")
save(example_rprimer_profile, file = "example_rprimer_profile.RData")
save(example_rprimer_properties, file = "example_rprimer_properties.RData")
save(example_rprimer_oligo, file = "example_rprimer_oligo.RData")
save(example_rprimer_assay, file = "example_rprimer_assay.RData")

tools::resaveRdaFiles(".")
