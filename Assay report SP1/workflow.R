# How to generate the SP1 assay (valdigt mkt efterhandskonstruktion=)

# Step 1: Import alignment
my_alignment <- "example_alignment.txt" %>%
  read_fasta_alignment %>%
  remove_gaps(., threshold = 0.5)  %>%
  select_roi(., from = 4000, to = 6000)


# Step 2: Get the sequence profile
my_sequence_profile <- sequence_profile(my_alignment)

# Step 3: Get sequence properties
my_sequence_properties <- sequence_properties(
  my_sequence_profile, iupac_threshold = 0.1
)

# Step 4: Get primers
my_primers <- get_oligos(
  my_sequence_properties,
  target = my_alignment,
  max_gap_frequency = 0.05,
  length = 18:24,
  max_degenerates = 1,
  max_degeneracy = 2,
  avoid_3end_ta = TRUE,
  avoid_3end_runs = TRUE,
  avoid_gc_rich_3end = FALSE,
  avoid_5end_g = FALSE,
  gc_range = c(0.40, 0.62),
  tm_range = c(50, 70),
  conc_oligo = 5e-07,
  conc_na = 0.05
)

# Step 5: Get assays
my_assays <- get_assays(my_primers, length = 70, max_tm_difference = 10)

# Step 6: Get probes
my_probes <- get_oligos(
  my_sequence_properties,
  target = my_alignment,
  max_gap_frequency = 0.05,
  length = 19,
  max_degenerates = 0,
  max_degeneracy = 1,
  avoid_3end_ta = FALSE,
  avoid_3end_runs = FALSE,
  avoid_gc_rich_3end = FALSE,
  avoid_5end_g = TRUE,
  gc_range = c(0.40, 0.65),
  tm_range = c(50, 70),
  conc_oligo = 2.5e-07,
  conc_na = 0.05
)

# Step 7: Add probes to assays
my_assays <- add_probes(my_assays, my_probes, tm_difference = c(-3, 1))

# Step 8: Save the data (if you want to)
# rp_save(my_alignment, filename = "my_alignment")
# rp_save(my_sequence_properties, filename = "my_sequence_properties")
# rp_save(my_assays, filename = "my_assays")

# Step 9: Generate a report

# Valj SP1=)
selected_assay <- my_assays[
  which(my_assays$majority_rev == "aggggttggttggatgaatatag")
, ]
selected_assay <- selected_assay[
  which(selected_assay$iupac_fwd == "rgtggtttctggggtgac")
, ]
selected_assay <- selected_assay[
  which(selected_assay$iupac_pr == "cgaagggctgagaatcaac")
, ]

# You can either select one assay and generate a report...
# selected_assay <- my_assays[1, ]

write_report(
 filename = "assay_report_SP1",
 selected_assay,
 my_sequence_profile,
 my_sequence_properties,
 comment = "my new hepatitis E virus assay :) (SP1)"
)

# ...or write several reports, e.g. one for each assay candidate:
# purrr::walk(seq_len(nrow(my_assays)), function(i) {
#  write_report(
#    filename = paste0("my_assay_report_", i),
#    my_assays[i, ],
#    my_sequence_profile,
#    my_sequence_properties,
#    comment = paste("my new hepatitis E virus assay, number", i)
#  )
#})
