#!/usr/bin/env Rscript

# step2_filter_forcecalled_snvs.R

# Purpose:
#   Convert the force-called patient VCF produced by Step 1 into a high-confidence candidate somatic SNV dataset.

# The script:
# 1. extracts per-sample AD and DP with bcftools
# 2. validates sample metadata
# 3. calculates REF/ALT counts, allelic depth and VAF
# 4. applies a matched-blood filter
# 5. defines mutation presence consistently in tissue samples
# 6. retains variants present in >=1 tumor/adjacent sample
# 7. writes long-format data and QC tables

# Usage:

# Rscript step2_filter_forcecalled_snvs.R \
#   PATIENT_ID \
#   FORCECALLED_VCF \
#   SAMPLE_METADATA_TSV \
#   OUTPUT_DIR

# Example:

# Rscript step2_filter_forcecalled_snvs.R \
#   ABCDE \
#   ABCDE_mpileup.vcf.gz \
#   config/sample_metadata.tsv \
#   results/02_filtered_snvs

# Required sample metadata columns:
#
#   patient
#   vcf_sample_id
#   sample_id
#   role


suppressPackageStartupMessages({
  library(tidyverse)
})


# 1. Command-line arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste(
      "Usage:",
      "Rscript step2_filter_forcecalled_snvs.R",
      "PATIENT_ID FORCECALLED_VCF SAMPLE_METADATA_TSV OUTPUT_DIR"
    ),
    call. = FALSE
  )
}

PATIENT   <- args[[1]]
VCF_PATH  <- normalizePath(args[[2]], mustWork = TRUE)
META_PATH <- normalizePath(args[[3]], mustWork = TRUE)
OUT_DIR   <- args[[4]]

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Filtering thresholds

# Matched-blood filter

# A candidate somatic mutation must have:
#   - total depth (FORMAT/DP) >=10 in blood
#   - <=2 ALT-supporting reads in blood
#   - blood VAF <2%

BLOOD_VAF_MAX <- 0.02
BLOOD_ALT_MAX <- 2L
BLOOD_DP_MIN <- 10L


# Tissue presence

# A mutation is considered present in a tumor or adjacent sample
# when:
#   - total depth (FORMAT/DP) >=20
#   - >=4 ALT-supporting reads
#   - VAF >=5%

TISSUE_VAF_MIN <- 0.05
TISSUE_ALT_MIN <- 4L
TISSUE_DP_MIN <- 20L


# 3. Helper functions

parse_ad_component <- function(x, index) {
  
  if (is.na(x) || x == "" || x == ".") {
    return(NA_integer_)
  }
  
  values <- strsplit(
    x,
    ",",
    fixed = TRUE
  )[[1]]
  
  if (length(values) != 2L) {
    return(NA_integer_)
  }
  
  suppressWarnings(
    as.integer(values[[index]])
  )
}


safe_median <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  median(
    x,
    na.rm = TRUE
  )
}


# 4. Check bcftools availability

BCFTOOLS <- Sys.which("bcftools")

if (!nzchar(BCFTOOLS)) {
  stop(
    "bcftools was not found in PATH.",
    call. = FALSE
  )
}


# 5. Load and validate sample metadata

meta_all <- read_tsv(
  META_PATH,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

required_columns <- c(
  "patient",
  "vcf_sample_id",
  "sample_id",
  "role"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_all)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}


meta <- meta_all %>%
  mutate(
    patient       = str_trim(patient),
    vcf_sample_id = str_trim(vcf_sample_id),
    sample_id     = str_trim(sample_id),
    role          = str_to_lower(str_trim(role))
  ) %>%
  filter(patient == PATIENT)


if (nrow(meta) == 0) {
  stop(
    "No metadata found for patient ",
    PATIENT,
    call. = FALSE
  )
}


allowed_roles <- c(
  "blood",
  "tumor",
  "adjacent"
)

bad_roles <- setdiff(
  unique(meta$role),
  allowed_roles
)

if (length(bad_roles) > 0) {
  stop(
    "Unexpected role(s): ",
    paste(bad_roles, collapse = ", "),
    call. = FALSE
  )
}


if (anyDuplicated(meta$vcf_sample_id)) {
  stop(
    "Duplicated vcf_sample_id values for ",
    PATIENT,
    call. = FALSE
  )
}


if (anyDuplicated(meta$sample_id)) {
  stop(
    "Duplicated sample_id values for ",
    PATIENT,
    call. = FALSE
  )
}


blood_meta <- meta %>%
  filter(role == "blood")

if (nrow(blood_meta) != 1) {
  stop(
    "Expected exactly one matched blood sample for ",
    PATIENT,
    "; found ",
    nrow(blood_meta),
    call. = FALSE
  )
}


tissue_meta <- meta %>%
  filter(
    role %in% c(
      "tumor",
      "adjacent"
    )
  )

if (nrow(tissue_meta) == 0) {
  stop(
    "No tumor or adjacent samples found for ",
    PATIENT,
    call. = FALSE
  )
}


cat("Patient:      ", PATIENT, "\n")
cat("Blood sample: ", blood_meta$sample_id, "\n")
cat(
  "Tissues:      ",
  paste(tissue_meta$sample_id, collapse = ", "),
  "\n"
)


# 6. Verify VCF sample names

sample_err <- tempfile(fileext = ".log")

vcf_samples <- system2(
  BCFTOOLS,
  args = c(
    "query",
    "-l",
    shQuote(VCF_PATH)
  ),
  stdout = TRUE,
  stderr = sample_err
)

sample_status <- attr(
  vcf_samples,
  "status"
)

if (is.null(sample_status)) {
  sample_status <- 0L
}

if (sample_status != 0L) {
  
  err <- readLines(
    sample_err,
    warn = FALSE
  )
  
  stop(
    "Could not read sample names from VCF:\n",
    paste(err, collapse = "\n"),
    call. = FALSE
  )
}

unlink(sample_err)

vcf_samples <- str_trim(vcf_samples)
vcf_samples <- vcf_samples[nzchar(vcf_samples)]


missing_metadata <- setdiff(
  vcf_samples,
  meta$vcf_sample_id
)

missing_vcf <- setdiff(
  meta$vcf_sample_id,
  vcf_samples
)


if (length(missing_metadata) > 0) {
  stop(
    "VCF sample(s) absent from metadata: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}


if (length(missing_vcf) > 0) {
  stop(
    "Metadata sample(s) absent from VCF: ",
    paste(missing_vcf, collapse = ", "),
    call. = FALSE
  )
}


cat(
  "VCF samples matched to metadata:",
  length(vcf_samples),
  "\n"
)

# 7. Extract CHROM/POS/REF/ALT + AD + DP

query_tmp <- tempfile(fileext = ".tsv")
query_err <- tempfile(fileext = ".log")

query_format <-
  "[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%AD\\t%DP\\n]"


query_status <- system2(
  BCFTOOLS,
  args = c(
    "query",
    "-f",
    shQuote(query_format),
    shQuote(VCF_PATH)
  ),
  stdout = query_tmp,
  stderr = query_err
)


if (query_status != 0L) {
  
  err <- readLines(
    query_err,
    warn = FALSE
  )
  
  stop(
    "bcftools query failed:\n",
    paste(err, collapse = "\n"),
    call. = FALSE
  )
}


raw <- read_tsv(
  query_tmp,
  
  col_names = c(
    "chr",
    "pos",
    "ref",
    "alt",
    "vcf_sample_id",
    "AD",
    "total_depth"
  ),
  
  col_types = cols(
    chr           = col_character(),
    pos           = col_integer(),
    ref           = col_character(),
    alt           = col_character(),
    vcf_sample_id = col_character(),
    AD            = col_character(),
    total_depth   = col_integer()
  ),
  
  na = c(
    "",
    "."
  ),
  
  show_col_types = FALSE
)

unlink(
  c(
    query_tmp,
    query_err
  )
)


if (nrow(raw) == 0) {
  stop(
    "No variant records extracted from ",
    VCF_PATH,
    call. = FALSE
  )
}


cat(
  "Variant x sample rows extracted:",
  nrow(raw),
  "\n"
)

# 8. Validate Step 1 output

variant_definitions <- raw %>%
  distinct(
    chr,
    pos,
    ref,
    alt
  )


bad_variants <- variant_definitions %>%
  filter(
    is.na(ref) |
      is.na(alt) |
      nchar(ref) != 1 |
      nchar(alt) != 1 |
      !ref %in% c("A", "C", "G", "T") |
      !alt %in% c("A", "C", "G", "T") |
      ref == alt
  )


if (nrow(bad_variants) > 0) {
  
  print(
    head(
      bad_variants,
      20
    )
  )
  
  stop(
    "Step 1 output contains unexpected non-biallelic/non-SNV records.",
    call. = FALSE
  )
}


# 9. Parse AD

ad_lengths <- map_int(
  raw$AD,
  function(x) {
    
    if (is.na(x) || x == "" || x == ".") {
      return(NA_integer_)
    }
    
    length(
      strsplit(
        x,
        ",",
        fixed = TRUE
      )[[1]]
    )
  }
)


bad_ad <- which(
  !is.na(ad_lengths) &
    ad_lengths != 2L
)


if (length(bad_ad) > 0) {
  
  print(
    raw[
      head(bad_ad, 20),
      c(
        "chr",
        "pos",
        "ref",
        "alt",
        "vcf_sample_id",
        "AD"
      )
    ]
  )
  
  stop(
    "Unexpected FORMAT/AD structure. ",
    "Expected exactly REF,ALT counts.",
    call. = FALSE
  )
}


raw <- raw %>%
  mutate(
    
    ref_count = map_int(
      AD,
      parse_ad_component,
      index = 1L
    ),
    
    alt_count = map_int(
      AD,
      parse_ad_component,
      index = 2L
    )
  )

# 10. Attach sample metadata

long <- raw %>%
  left_join(
    meta %>%
      select(
        vcf_sample_id,
        sample_id,
        role
      ),
    by = "vcf_sample_id"
  )


if (any(is.na(long$sample_id))) {
  stop(
    "At least one VCF sample could not be assigned metadata.",
    call. = FALSE
  )
}

# 11. Calculate mutation IDs, informative depth and VAF

long <- long %>%
  mutate(
    
    patient = PATIENT,
    
    mutation_id = paste(
      chr,
      pos,
      ref,
      alt,
      sep = ":"
    ),
    
    # Number of reads specifically assigned to REF or the target ALT allele

    allelic_depth = if_else(
      !is.na(ref_count) &
        !is.na(alt_count),
      
      ref_count + alt_count,
      
      NA_integer_
    ),
    
    # VAF is based on reads assigned to REF or the target ALT
    vaf = if_else(
      !is.na(allelic_depth) &
        allelic_depth > 0 &
        !is.na(alt_count),
      
      alt_count / allelic_depth,
      
      NA_real_
    )
  )


# 12. Verify complete mutation x sample matrix

expected_samples <- nrow(meta)

matrix_check <- long %>%
  count(
    mutation_id,
    name = "n_samples"
  ) %>%
  filter(
    n_samples != expected_samples
  )


if (nrow(matrix_check) > 0) {
  
  print(
    head(
      matrix_check,
      20
    )
  )
  
  stop(
    "Some mutations do not contain exactly one row per patient sample.",
    call. = FALSE
  )
}


duplicate_check <- long %>%
  count(
    mutation_id,
    sample_id
  ) %>%
  filter(n != 1)


if (nrow(duplicate_check) > 0) {
  
  print(
    head(
      duplicate_check,
      20
    )
  )
  
  stop(
    "Duplicated mutation x sample rows detected.",
    call. = FALSE
  )
}


# 13. Define tissue mutation presence

long <- long %>%
  mutate(
    
    present =
      
      role %in% c(
        "tumor",
        "adjacent"
      ) &
      
      coalesce(
        total_depth >= TISSUE_DP_MIN,
        FALSE
      ) &
      
      coalesce(
        alt_count >= TISSUE_ALT_MIN,
        FALSE
      ) &
      
      coalesce(
        vaf >= TISSUE_VAF_MIN,
        FALSE
      )
  )
# 14. Matched-blood filter

blood <- long %>%
  filter(
    role == "blood"
  ) %>%
  mutate(
    
    pass_blood =
      
      coalesce(
        total_depth >= BLOOD_DP_MIN,
        FALSE
      ) &
      
      coalesce(
        alt_count <= BLOOD_ALT_MAX,
        FALSE
      ) &
      
      coalesce(
        vaf < BLOOD_VAF_MAX,
        FALSE
      ),
    
    blood_fail_reason = case_when(
  
  is.na(total_depth) ~
    "missing_blood_DP",
  
  total_depth < BLOOD_DP_MIN ~
    "low_blood_depth",
  
  is.na(alt_count) |
    is.na(vaf) ~
    "missing_blood_AD",
  
  alt_count > BLOOD_ALT_MAX &
    vaf >= BLOOD_VAF_MAX ~
    "blood_alt_and_vaf",
  
  alt_count > BLOOD_ALT_MAX ~
    "blood_alt_count",
  
  vaf >= BLOOD_VAF_MAX ~
    "blood_vaf",
  
  TRUE ~
    "pass")
  ) %>%
  select(
    mutation_id,
    
    blood_ref_count = ref_count,
    blood_alt_count = alt_count,
    blood_allelic_depth = allelic_depth,
    blood_total_depth = total_depth,
    blood_vaf = vaf,
    
    pass_blood,
    blood_fail_reason
  )

# 15. Summarise tissue evidence

tissue_summary <- long %>%
  filter(
    role %in% c(
      "tumor",
      "adjacent"
    )
  ) %>%
  group_by(
    patient,
    mutation_id,
    chr,
    pos,
    ref,
    alt
  ) %>%
  summarise(
    
    n_tissue_samples = n(),
    
    n_tissues_present =
      sum(
        present,
        na.rm = TRUE
      ),
    
    n_tumor_samples =
      sum(role == "tumor"),
    
    n_tumor_present =
      sum(
        role == "tumor" &
          present,
        na.rm = TRUE
      ),
    
    n_adjacent_samples =
      sum(role == "adjacent"),
    
    n_adjacent_present =
      sum(
        role == "adjacent" &
          present,
        na.rm = TRUE
      ),
    
    present_in_tumor =
      any(
        role == "tumor" &
          present
      ),
    
    present_in_adjacent =
      any(
        role == "adjacent" &
          present
      ),
    
    pass_tissue_signal =
      any(present),
    
    .groups = "drop"
  )


# 16. Final candidate-somatic classification

variant_summary <- tissue_summary %>%
  left_join(
    blood,
    by = "mutation_id"
  ) %>%
  mutate(
    
    is_candidate_somatic =
      pass_blood &
      pass_tissue_signal,
    
    filter_status = case_when(
      
      !pass_blood ~
        paste0(
          "excluded_",
          blood_fail_reason
        ),
      
      !pass_tissue_signal ~
        "excluded_no_tissue_signal",
      
      TRUE ~
        "candidate_somatic"
    )
  )

# 17. Add variant-level filtering information to long data

long_flagged <- long %>%
  left_join(
    variant_summary %>%
      select(
        mutation_id,
        pass_blood,
        blood_fail_reason,
        pass_tissue_signal,
        n_tissues_present,
        n_tumor_present,
        n_adjacent_present,
        present_in_tumor,
        present_in_adjacent,
        is_candidate_somatic,
        filter_status
      ),
    by = "mutation_id"
  )


candidate_long <- long_flagged %>%
  filter(
    is_candidate_somatic
  )


# 18. Create wide table

wide_counts <- long %>%
  select(
    mutation_id,
    sample_id,
    ref_count,
    alt_count,
    allelic_depth,
    total_depth,
    vaf,
    present
  ) %>%
  pivot_wider(
    
    names_from = sample_id,
    
    values_from = c(
      ref_count,
      alt_count,
      allelic_depth,
      total_depth,
      vaf,
      present
    ),
    
    names_glue =
      "{sample_id}_{.value}"
  )


candidate_wide <- variant_summary %>%
  filter(
    is_candidate_somatic
  ) %>%
  left_join(
    wide_counts,
    by = "mutation_id"
  )

# 19. QC summaries

sample_qc <- long_flagged %>%
  group_by(
    patient,
    sample_id,
    role
  ) %>%
  summarise(
    
    n_target_variants = n(),
    
    n_with_AD =
      sum(
        !is.na(ref_count) &
          !is.na(alt_count)
      ),
    
    n_with_total_DP =
      sum(
        !is.na(total_depth)
      ),
    
    median_allelic_depth =
      safe_median(
        allelic_depth
      ),
    
    median_total_depth =
      safe_median(
        total_depth
      ),
    
    n_with_alt_support =
      sum(
        alt_count > 0,
        na.rm = TRUE
      ),
    
    n_present =
      sum(
        present,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


filter_qc <- tibble(
  
  patient = PATIENT,
  
  n_input_snvs =
    nrow(variant_summary),
  
  n_pass_blood =
    sum(
      variant_summary$pass_blood,
      na.rm = TRUE
    ),
  
  n_with_tissue_signal =
    sum(
      variant_summary$pass_tissue_signal,
      na.rm = TRUE
    ),
  
  n_candidate_somatic =
    sum(
      variant_summary$is_candidate_somatic,
      na.rm = TRUE
    )
)


blood_filter_qc <- blood %>%
  count(
    blood_fail_reason,
    name = "n_variants"
  ) %>%
  arrange(
    desc(n_variants)
  )


thresholds <- tribble(
  
  ~parameter,
  ~value,
  
  "blood_vaf_max",
  BLOOD_VAF_MAX,
  
  "blood_alt_max",
  BLOOD_ALT_MAX,
  
  "blood_DP_min",
  BLOOD_DP_MIN,
  
  "tissue_vaf_min",
  TISSUE_VAF_MIN,
  
  "tissue_alt_min",
  TISSUE_ALT_MIN,
  
  "tissue_DP_min",
  TISSUE_DP_MIN
)

# 20. Write tables

write_tsv(
  long_flagged,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_all_variants_long.tsv"
    )
  )
)


write_tsv(
  variant_summary,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_variant_summary.tsv"
    )
  )
)


write_tsv(
  candidate_long,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_candidate_somatic_long.tsv"
    )
  )
)


write_tsv(
  candidate_wide,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_candidate_somatic_wide.tsv"
    )
  )
)


write_tsv(
  sample_qc,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_sample_qc.tsv"
    )
  )
)


write_tsv(
  filter_qc,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_filter_qc.tsv"
    )
  )
)


write_tsv(
  blood_filter_qc,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_blood_filter_qc.tsv"
    )
  )
)


write_tsv(
  thresholds,
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_filter_thresholds.tsv"
    )
  )
)

# 21. VAF QC plot

plot_data <- candidate_long %>%
  filter(
    role %in% c(
      "tumor",
      "adjacent"
    ),
    !is.na(vaf)
  )


if (nrow(plot_data) > 0) {
  
  p <- ggplot(
    plot_data,
    aes(x = vaf)
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    geom_vline(
      xintercept = TISSUE_VAF_MIN,
      linetype = "dashed"
    ) +
    
    facet_wrap(
      ~ sample_id,
      ncol = 1,
      scales = "free_y"
    ) +
    
    labs(
      title = paste(
        PATIENT,
        "candidate somatic SNV VAF distributions"
      ),
      
      subtitle = paste0(
      "Presence: VAF >= ",
      TISSUE_VAF_MIN,
      ", ALT >= ",
      TISSUE_ALT_MIN,
      ", DP >= ",
      TISSUE_DP_MIN
    ),
      
      x = "Variant allele frequency",
      y = "Number of SNVs"
    ) +
    
    theme_minimal()
  
  
  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        PATIENT,
        "_candidate_somatic_vaf_qc.pdf"
      )
    ),
    
    p,
    
    width = 8,
    
    height = max(
      5,
      2 * n_distinct(
        plot_data$sample_id
      )
    )
  )
}


# 22. Record software versions

bcftools_version <- system2(
  BCFTOOLS,
  "--version",
  stdout = TRUE,
  stderr = TRUE
)


writeLines(
  c(
    paste0(
      "Patient: ",
      PATIENT
    ),
    
    paste0(
      "Input VCF: ",
      VCF_PATH
    ),
    
    "",
    
    "bcftools:",
    bcftools_version,
    
    "",
    
    "R session:",
    capture.output(
      sessionInfo()
    )
  ),
  
  file.path(
    OUT_DIR,
    paste0(
      PATIENT,
      "_software_versions.txt"
    )
  )
)

# 23. Final report

cat("\n")

cat("Step 2 complete:", PATIENT, "\n")


cat(
  "Input SNVs:             ",
  nrow(variant_summary),
  "\n"
)

cat(
  "Pass matched blood:     ",
  sum(
    variant_summary$pass_blood,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Signal in >=1 tissue:   ",
  sum(
    variant_summary$pass_tissue_signal,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Candidate somatic SNVs: ",
  nrow(candidate_wide),
  "\n"
)

cat(
  "Output directory:       ",
  OUT_DIR,
  "\n"
)