
# script1_calculate_sCNA_burden.R

# Calculate cell-fraction-weighted somatic copy number alteration (sCNA) burden from MoChA calls.

# Burden is calculated as the sum(event length × estimated cell fraction) / genomic length × 100

# Burdens are calculated:
#   1. genome-wide across autosomes 1-22
#   2. separately for each autosome
#   3. for all alterations combined
#   4. separately for gains, losses, and CN-LOH

# Samples with no retained alterations are assigned a burden of zero.

# Input:
#   mocha_calls.csv

# Output:
#   scna_burden_per_sample.csv

# Packages

library(tidyverse)

# Input

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

# Ensure the directory exists (creates it if it doesn't)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

mocha_calls <- read_csv(file.path(INPUT_DIR, "mocha_calls.csv"),  show_col_types = FALSE)


# Reference genome
#
# GRCh37 / hg19 chromosome lengths for autosomes 1-22.
# The genome-wide denominator is calculated from the same chromosome lengths
# used for chromosome-specific burden calculations.

chr_lengths <- tibble(
  seqnames = paste0("chr", 1:22),
  chr_length = c(
    249250621,
    243199373,
    198022430,
    191154276,
    180915260,
    171115067,
    159138663,
    146364022,
    141213431,
    135534747,
    135006516,
    133851895,
    115169878,
    107349540,
    102531392,
    90354753,
    81195210,
    78077248,
    59128983,
    63025520,
    48129895,
    51304566
  )
)

total_autosome_length <- sum(chr_lengths$chr_length)

# GRCh37 autosomes 1-22:
stopifnot(total_autosome_length == 2881033286)


# Sample manifest
#
# Create this before filtering alteration calls so that samples with no qualifying alterations remain in the final dataset with burden = 0.

all_samples <- mocha_calls %>%
  distinct(
    sample_id,
    Donor,
    Tissue,
    isCancer3
  ) %>%
  rename(tissue_class = isCancer3)


# Prepare retained sCNA calls
#
# Detailed MoChA QC and germline-event removal are described in the Methods.
# mocha.nongermline identifies calls retained as non-germline for burden
# calculation.
#
# Only autosomal Gain, Loss, and CN-LOH events are included.

scna_calls <- mocha_calls %>%
  filter(
    mocha.nongermline == TRUE,
    seqnames %in% chr_lengths$seqnames,
    type %in% c("Gain", "Loss", "CN-LOH")
  ) %>%
  mutate(
    # MoChA reports event length directly in base pairs.
    event_length_bp = parse_number(as.character(length)),
    cf = as.numeric(cf)
  )


# Basic integrity checks

if (anyNA(scna_calls$event_length_bp)) {
  stop("Missing event lengths detected among retained sCNA calls.")
}

if (any(scna_calls$event_length_bp <= 0)) {
  stop("Non-positive event lengths detected among retained sCNA calls.")
}

if (anyNA(scna_calls$cf)) {
  stop("Missing cell-fraction estimates detected among retained sCNA calls.")
}

if (any(scna_calls$cf < 0 | scna_calls$cf > 1)) {
  stop("Cell-fraction estimates outside the expected range [0, 1] detected.")
}


# Calculate cell-fraction-weighted event length

scna_calls <- scna_calls %>%
  mutate(
    cf_weighted_length_bp = event_length_bp * cf
  )


# Genome-wide burden

genome_burden <- scna_calls %>%
  group_by(
    sample_id,
    Donor,
    Tissue
  ) %>%
  summarise(
    n_events = n(),
    n_gain = sum(type == "Gain"),
    n_loss = sum(type == "Loss"),
    n_cnloh = sum(type == "CN-LOH"),
    
    cf_weighted_total_bp =
      sum(cf_weighted_length_bp),
    
    cf_weighted_gain_bp =
      sum(cf_weighted_length_bp[type == "Gain"]),
    
    cf_weighted_loss_bp =
      sum(cf_weighted_length_bp[type == "Loss"]),
    
    cf_weighted_cnloh_bp =
      sum(cf_weighted_length_bp[type == "CN-LOH"]),
    
    .groups = "drop"
  ) %>%
  mutate(
    burden_total_pct =
      100 * cf_weighted_total_bp / total_autosome_length,
    
    burden_gain_pct =
      100 * cf_weighted_gain_bp / total_autosome_length,
    
    burden_loss_pct =
      100 * cf_weighted_loss_bp / total_autosome_length,
    
    burden_cnloh_pct =
      100 * cf_weighted_cnloh_bp / total_autosome_length
  )


# Retain samples with no qualifying alterations

genome_burden <- all_samples %>%
  left_join(
    genome_burden,
    by = c("sample_id", "Donor", "Tissue")
  ) %>%
  mutate(
    across(
      c(
        n_events,
        n_gain,
        n_loss,
        n_cnloh
      ),
      ~ replace_na(.x, 0L)
    ),
    across(
      c(
        cf_weighted_total_bp,
        cf_weighted_gain_bp,
        cf_weighted_loss_bp,
        cf_weighted_cnloh_bp,
        burden_total_pct,
        burden_gain_pct,
        burden_loss_pct,
        burden_cnloh_pct
      ),
      ~ replace_na(.x, 0)
    )
  )


# Chromosome-specific burden

chr_burden <- scna_calls %>%
  left_join(
    chr_lengths,
    by = "seqnames"
  ) %>%
  group_by(
    sample_id,
    Donor,
    Tissue,
    seqnames
  ) %>%
  summarise(
    chr_length = first(chr_length),
    
    total_pct =
      100 * sum(cf_weighted_length_bp) /
      first(chr_length),
    
    gain_pct =
      100 * sum(cf_weighted_length_bp[type == "Gain"]) /
      first(chr_length),
    
    loss_pct =
      100 * sum(cf_weighted_length_bp[type == "Loss"]) /
      first(chr_length),
    
    cnloh_pct =
      100 * sum(cf_weighted_length_bp[type == "CN-LOH"]) /
      first(chr_length),
    
    .groups = "drop"
  )


# Generate every sample × chromosome combination.
# This guarantees explicit zeros for unaffected chromosomes.

sample_chr_grid <- all_samples %>%
  select(
    sample_id,
    Donor,
    Tissue
  ) %>%
  crossing(chr_lengths)


chr_burden <- sample_chr_grid %>%
  left_join(
    chr_burden,
    by = c(
      "sample_id",
      "Donor",
      "Tissue",
      "seqnames",
      "chr_length"
    )
  ) %>%
  mutate(
    across(
      c(
        total_pct,
        gain_pct,
        loss_pct,
        cnloh_pct
      ),
      ~ replace_na(.x, 0)
    )
  ) %>%
  select(
    -chr_length
  )


# Convert chromosome-specific results to one row per sample

chr_burden_wide <- chr_burden %>%
  pivot_wider(
    id_cols = c(
      sample_id,
      Donor,
      Tissue
    ),
    names_from = seqnames,
    values_from = c(
      total_pct,
      gain_pct,
      loss_pct,
      cnloh_pct
    ),
    names_glue = "{seqnames}_{.value}"
  )


# Combine genome-wide and chromosome-specific burdens

burden_all <- genome_burden %>%
  left_join(
    chr_burden_wide,
    by = c(
      "sample_id",
      "Donor",
      "Tissue"
    )
  )


# Internal consistency check
#
# Total burden should equal gain + loss + CN-LOH burden, apart from
# floating-point precision.

stopifnot(
  all(
    near(
      burden_all$burden_total_pct,
      burden_all$burden_gain_pct +
        burden_all$burden_loss_pct +
        burden_all$burden_cnloh_pct
    )
  )
)


# Summary for inspection

burden_summary <- burden_all %>%
  group_by(tissue_class) %>%
  summarise(
    n_samples = n(),
    median_total_pct = median(burden_total_pct),
    median_gain_pct = median(burden_gain_pct),
    median_loss_pct = median(burden_loss_pct),
    median_cnloh_pct = median(burden_cnloh_pct),
    .groups = "drop"
  )

print(burden_summary)


# Save the final file

write_csv(burden_all, file.path(OUT_DIR, "scna_burden_per_sample.csv"))