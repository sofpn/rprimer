
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rprimer

<!-- badges: start -->

<!-- badges: end -->

rprimer is designed to automate primer, probe and PCR assay design as
much as possible.

## Installation

You can install the development version of rprimer from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sofpn/rprimer")
```

## Example

The primer, probe and assay design workflow is summarized below.

``` r
library(rprimer)

# Step 1: Import alignment 
my_alignment <- "example_alignment.txt" %>% 
  read_fasta_alignment %>%
   remove_gaps(., threshold = 0.5)  %>% 
    select_roi(., from = 4000, to = 6000) 

# Step 2: Get the sequence profile 
my_sequence_profile <- sequence_profile(my_alignment)

# Step 3: Get sequence properties 
my_sequence_properties <- sequence_properties(
  my_sequence_profile, iupac_threshold = 0.05 
) 

# Step 4: Get primers
my_primers <- get_oligos(
  my_sequence_properties, 
  target = my_alignment,
  max_gap_frequency = 0.05, 
  length = 18:22,
  max_degenerates = 1,
  max_degeneracy = 2, 
  avoid_3end_ta = TRUE, 
  avoid_3end_runs = TRUE,
  avoid_gc_rich_3end = TRUE,
  avoid_5end_g = FALSE,
  gc_range = c(0.45, 0.56),
  tm_range = c(55, 70), 
  conc_oligo = 5e-07, 
  conc_na = 0.05   
)

# Step 5: Get assays 
my_assays <- get_assays(my_primers, length = 65:75, max_tm_difference = 1.5) 
 
# Step 6: Get probes 
my_probes <- get_oligos(
  my_sequence_properties,
  target = my_alignment,
  max_gap_frequency = 0.05,
  length = 18:24,
  max_degenerates = 1,
  max_degeneracy = 2,
  avoid_3end_ta = FALSE,
  avoid_3end_runs = FALSE,
  avoid_gc_rich_3end = FALSE,
  avoid_5end_g = TRUE,
  gc_range = c(0.45, 0.56),
  tm_range = c(55, 70),
  conc_oligo = 2.5e-07, 
  conc_na = 0.05
)

# Step 7: Add probes to assays
my_assays <- add_probes(my_assays, my_probes, tm_difference = c(-1, 5))  

# Step 8: Save the data (if you want to)
# rp_save(my_alignment, filename = "my_alignment")
# rp_save(my_sequence_properties, filename = "my_sequence_properties")
# rp_save(my_assays, filename = "my_assays")

# Step 9: Generate a report

# You can either select one assay and generate a report...
selected_assay <- my_assays[1, ]

# write_report(
#  filename = "my_assay_report",
#  selected_assay,
#  my_sequence_profile,
#  my_sequence_properties,
#  comment = "my new hepatitis E virus assay :)"
# )
```

The package vignette contains a lot more information\!
