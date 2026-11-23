# script4_prepare_pyclone_input.R

# Purpose:
#   Merge candidate somatic SNV allele counts from Step 2 with
#   ASCAT allele-specific copy number and purity estimates to
#   generate per-patient PyClone-VI input files.

# Strategy:
#   Tumor samples:
#     - require usable ASCAT purity
#     - require ASCAT segmentation
#     - use ASCAT major/minor copy number

#   Adjacent/NAT samples:
#     - tumour_content is set to 1.0
#     - use ASCAT major/minor copy number when available
#     - if no usable ASCAT segment covers a locus, assume diploid copy number: major_cn = 1, minor_cn = 1

#   No GOF/MAPD thresholds are used for sample exclusion.
#   These metrics are retained for QC only.

# PyClone-VI input columns:
#   mutation_id
#   sample_id
#   ref_counts
#   alt_counts
#   normal_cn
#   major_cn
#   minor_cn
#   tumour_content

# Intended use:
#   Open the script in RStudio, set paths and run

# Required metadata columns:
#   patient
#   sample_id
#   role

# Required sample ID recoding columns (sample_id_recoding.tsv):
#   source_id
#   canonical_sample_id


# 1. Packages

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(GenomicRanges)
  library(IRanges)
})


# 2. Paths

STEP2_DIR <- "path/to/step2/output/folder"
META_PATH <- "path/to/sample_metadata.tsv"
RECODE_PATH <- "path/to/sample_id_recoding.tsv"
ASCAT_SUMMARY_PATH <- "path/to/ASCAT/summary/file.csv"
ASCAT_SEG_DIR <- "path/to/ASCAT_sample_files/"
OUT_DIR <- "path/to/output/directory"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 3. Analysis settings

KEEP_CHRX <- TRUE
KEEP_CHRY <- FALSE

ADJACENT_TUMOUR_CONTENT <- 1.0

ADJACENT_FALLBACK_MAJOR_CN <- 1L
ADJACENT_FALLBACK_MINOR_CN <- 1L


# 4. Helper functions

normalize_id <- function(x) {
  
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all(
      "-",
      "_"
    ) %>%
    stringr::str_replace_all(
      "_+",
      "_"
    )
}


normalize_ascat_id <- function(x) {
  
  normalize_id(x) %>%
    stringr::str_replace(
      "_[CN]_",
      "_"
    ) %>%
    stringr::str_replace_all(
      "_+",
      "_"
    )
}


normalize_chr <- function(x) {
  
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_remove(
      "^chr"
    )
}


normalize_sex <- function(x) {
  
  x <- stringr::str_to_lower(
    stringr::str_trim(
      as.character(x)
    )
  )
  
  dplyr::case_when(
    
    x %in% c(
      "m",
      "male",
      "xy"
    ) ~
      "male",
    
    x %in% c(
      "f",
      "female",
      "xx"
    ) ~
      "female",
    
    TRUE ~
      NA_character_
  )
}


parse_logical_strict <- function(
    x,
    variable_name
) {
  
  y <- stringr::str_to_lower(
    stringr::str_trim(
      as.character(x)
    )
  )
  
  result <- rep(
    NA,
    length(y)
  )
  
  result[
    y %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  ] <- TRUE
  
  result[
    y %in% c(
      "false",
      "f",
      "0",
      "no",
      "n"
    )
  ] <- FALSE
  
  invalid <-
    is.na(result) &
    !is.na(y) &
    y != ""
  
  if (any(invalid)) {
    
    stop(
      "Unexpected value(s) in ",
      variable_name,
      ": ",
      paste(
        unique(y[invalid]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  result
}


read_table_auto <- function(path) {
  
  extension <- stringr::str_to_lower(
    tools::file_ext(path)
  )
  
  if (extension == "csv") {
    
    readr::read_csv(
      path,
      col_types = readr::cols(
        .default = "c"
      ),
      show_col_types = FALSE
    )
    
  } else {
    
    readr::read_tsv(
      path,
      col_types = readr::cols(
        .default = "c"
      ),
      show_col_types = FALSE
    )
  }
}


apply_sample_recoding <- function(
    x,
    recoding
) {
  
  raw_key <-
    normalize_id(x)
  
  normalized_key <-
    normalize_ascat_id(x)
  
  
  raw_match <- match(
    raw_key,
    recoding$source_key
  )
  
  
  normalized_match <- match(
    normalized_key,
    recoding$source_key
  )
  
  
  mapped_key <-
    recoding$canonical_key[
      raw_match
    ]
  
  
  use_normalized_match <-
    is.na(mapped_key) &
    !is.na(normalized_match)
  
  
  mapped_key[
    use_normalized_match
  ] <-
    recoding$canonical_key[
      normalized_match[
        use_normalized_match
      ]
    ]
  
  
  dplyr::if_else(
    !is.na(mapped_key),
    mapped_key,
    normalized_key
  )
}


# 5. Validate input paths

if (!dir.exists(STEP2_DIR)) {
  
  stop(
    "Step 2 directory not found:\n",
    STEP2_DIR,
    call. = FALSE
  )
}


if (!file.exists(META_PATH)) {
  
  stop(
    "Sample metadata file not found:\n",
    META_PATH,
    call. = FALSE
  )
}


if (!file.exists(RECODE_PATH)) {
  
  stop(
    "Sample ID recoding file not found:\n",
    RECODE_PATH,
    call. = FALSE
  )
}


if (!file.exists(ASCAT_SUMMARY_PATH)) {
  
  stop(
    "ASCAT summary file not found:\n",
    ASCAT_SUMMARY_PATH,
    call. = FALSE
  )
}


if (!dir.exists(ASCAT_SEG_DIR)) {
  
  stop(
    "ASCAT segment directory not found:\n",
    ASCAT_SEG_DIR,
    call. = FALSE
  )
}


# 6. Load sample metadata

meta <- readr::read_tsv(
  META_PATH,
  col_types = readr::cols(
    .default = "c"
  ),
  show_col_types = FALSE
)


required_metadata <- c(
  "patient",
  "sample_id",
  "role"
)


missing_metadata <- setdiff(
  required_metadata,
  names(meta)
)


if (length(missing_metadata) > 0) {
  
  stop(
    "Missing metadata column(s): ",
    paste(
      missing_metadata,
      collapse = ", "
    ),
    call. = FALSE
  )
}


meta <- meta %>%
  dplyr::mutate(
    
    patient =
      stringr::str_trim(patient),
    
    sample_id =
      stringr::str_trim(sample_id),
    
    sample_key =
      normalize_id(sample_id),
    
    role =
      stringr::str_to_lower(
        stringr::str_trim(role)
      )
  )


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
    "Unexpected metadata role(s): ",
    paste(
      bad_roles,
      collapse = ", "
    ),
    call. = FALSE
  )
}


if (anyDuplicated(meta$sample_key)) {
  
  duplicated_metadata <- meta %>%
    dplyr::count(
      sample_key
    ) %>%
    dplyr::filter(
      n > 1L
    )
  
  print(
    duplicated_metadata
  )
  
  stop(
    "Duplicated normalized sample IDs detected in metadata.",
    call. = FALSE
  )
}


# 7. Load sample ID recoding table

recoding <- readr::read_tsv(
  RECODE_PATH,
  col_types = readr::cols(
    .default = "c"
  ),
  show_col_types = FALSE
)


required_recoding <- c(
  "source_id",
  "canonical_sample_id"
)


missing_recoding <- setdiff(
  required_recoding,
  names(recoding)
)


if (length(missing_recoding) > 0) {
  
  stop(
    "Missing recoding column(s): ",
    paste(
      missing_recoding,
      collapse = ", "
    ),
    call. = FALSE
  )
}


recoding <- recoding %>%
  dplyr::mutate(
    
    source_id =
      stringr::str_trim(source_id),
    
    canonical_sample_id =
      stringr::str_trim(
        canonical_sample_id
      ),
    
    source_key =
      normalize_id(source_id),
    
    canonical_key =
      normalize_id(
        canonical_sample_id
      )
  ) %>%
  dplyr::filter(
    source_key != "",
    canonical_key != ""
  ) %>%
  dplyr::select(
    source_id,
    canonical_sample_id,
    source_key,
    canonical_key
  )


if (anyDuplicated(recoding$source_key)) {
  
  duplicated_recoding <- recoding %>%
    dplyr::count(
      source_key
    ) %>%
    dplyr::filter(
      n > 1L
    )
  
  print(
    duplicated_recoding
  )
  
  stop(
    "Duplicated source IDs detected in sample_id_recoding.tsv.",
    call. = FALSE
  )
}


unused_recoding <- recoding %>%
  dplyr::filter(
    !canonical_key %in%
      meta$sample_key
  )


if (nrow(unused_recoding) > 0) {
  
  message(
    nrow(unused_recoding),
    " recoding entries refer to samples that are ",
    "not part of the WES analysis."
  )
}


readr::write_tsv(
  unused_recoding,
  file.path(
    OUT_DIR,
    "sample_id_recoding_unused.tsv"
  )
)


# 8. Load Step 2 candidate somatic files

step2_files <- list.files(
  STEP2_DIR,
  pattern =
    "_candidate_somatic_long\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)


step2_files <- sort(
  step2_files
)


if (length(step2_files) == 0) {
  
  stop(
    "No *_candidate_somatic_long.tsv files found under:\n",
    STEP2_DIR,
    call. = FALSE
  )
}


read_step2_file <- function(path) {
  
  expected_patient <-
    basename(path) %>%
    stringr::str_remove(
      "_candidate_somatic_long\\.tsv$"
    )
  
  
  dat <- readr::read_tsv(
    path,
    col_types = readr::cols(
      .default = "c"
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  if (
    nrow(dat) > 0 &&
    "patient" %in% names(dat)
  ) {
    
    observed_patient <- unique(
      stringr::str_trim(
        dat$patient
      )
    )
    
    observed_patient <-
      observed_patient[
        !is.na(observed_patient) &
          observed_patient != ""
      ]
    
    
    if (
      length(observed_patient) != 1L ||
      observed_patient != expected_patient
    ) {
      
      stop(
        "Patient identifier in Step 2 file does not match filename:\n",
        path,
        call. = FALSE
      )
    }
  }
  
  
  dat
}


candidate_long <- purrr::map_dfr(
  step2_files,
  read_step2_file
)


required_step2_columns <- c(
  "patient",
  "sample_id",
  "role",
  "chr",
  "pos",
  "ref",
  "alt",
  "ref_count",
  "alt_count",
  "total_depth",
  "present",
  "is_candidate_somatic"
)


missing_step2_columns <- setdiff(
  required_step2_columns,
  names(candidate_long)
)


if (length(missing_step2_columns) > 0) {
  
  stop(
    "Step 2 output is missing required column(s): ",
    paste(
      missing_step2_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}


candidate_long <- candidate_long %>%
  dplyr::mutate(
    
    patient =
      stringr::str_trim(patient),
    
    sample_id =
      stringr::str_trim(sample_id),
    
    sample_key =
      normalize_id(sample_id),
    
    role =
      stringr::str_to_lower(
        stringr::str_trim(role)
      ),
    
    chr =
      stringr::str_trim(chr),
    
    chr_plain =
      normalize_chr(chr),
    
    pos =
      suppressWarnings(
        as.integer(pos)
      ),
    
    ref =
      stringr::str_to_upper(
        stringr::str_trim(ref)
      ),
    
    alt =
      stringr::str_to_upper(
        stringr::str_trim(alt)
      ),
    
    ref_count =
      suppressWarnings(
        as.integer(ref_count)
      ),
    
    alt_count =
      suppressWarnings(
        as.integer(alt_count)
      ),
    
    total_depth =
      suppressWarnings(
        as.integer(total_depth)
      ),
    
    present =
      parse_logical_strict(
        present,
        "present"
      ),
    
    is_candidate_somatic =
      parse_logical_strict(
        is_candidate_somatic,
        "is_candidate_somatic"
      )
  )


if (!"mutation_id" %in%
    names(candidate_long)) {
  
  candidate_long <- candidate_long %>%
    dplyr::mutate(
      
      mutation_id =
        paste(
          chr,
          pos,
          ref,
          alt,
          sep = ":"
        )
    )
  
} else {
  
  candidate_long <- candidate_long %>%
    dplyr::mutate(
      mutation_id =
        stringr::str_trim(
          mutation_id
        )
    )
}


# 9. Retain candidate somatic tissue rows

wes <- candidate_long %>%
  dplyr::filter(
    
    is_candidate_somatic,
    
    role %in% c(
      "tumor",
      "adjacent"
    )
  )


if (nrow(wes) == 0) {
  
  stop(
    "No candidate somatic tumor/adjacent rows were found.",
    call. = FALSE
  )
}


# 10. Validate SNVs

bad_snvs <- wes %>%
  dplyr::filter(
    
    is.na(chr_plain) |
      is.na(pos) |
      pos < 1L |
      
      is.na(ref) |
      is.na(alt) |
      
      nchar(ref) != 1L |
      nchar(alt) != 1L |
      
      !ref %in% c(
        "A",
        "C",
        "G",
        "T"
      ) |
      
      !alt %in% c(
        "A",
        "C",
        "G",
        "T"
      ) |
      
      ref == alt
  )


if (nrow(bad_snvs) > 0) {
  
  print(
    head(
      bad_snvs,
      20
    )
  )
  
  stop(
    "Malformed or non-SNV records detected in Step 2 input.",
    call. = FALSE
  )
}


# 11. Validate WES samples against metadata

wes_samples <- wes %>%
  dplyr::distinct(
    patient,
    sample_id,
    sample_key,
    role
  )


metadata_check <- wes_samples %>%
  dplyr::left_join(
    
    meta %>%
      dplyr::select(
        patient,
        sample_key,
        metadata_sample_id = sample_id,
        metadata_role = role
      ),
    
    by = c(
      "patient",
      "sample_key"
    )
  )


missing_wes_metadata <- metadata_check %>%
  dplyr::filter(
    is.na(metadata_sample_id)
  )


if (nrow(missing_wes_metadata) > 0) {
  
  print(
    missing_wes_metadata
  )
  
  stop(
    "At least one Step 2 sample is absent from sample metadata.",
    call. = FALSE
  )
}


role_mismatches <- metadata_check %>%
  dplyr::filter(
    role != metadata_role
  )


if (nrow(role_mismatches) > 0) {
  
  print(
    role_mismatches
  )
  
  stop(
    "Sample-role mismatch between Step 2 and metadata.",
    call. = FALSE
  )
}


# 12. Apply chromosome policy

keep_chr <- as.character(
  1:22
)


if (KEEP_CHRX) {
  
  keep_chr <- c(
    keep_chr,
    "X"
  )
}


if (KEEP_CHRY) {
  
  keep_chr <- c(
    keep_chr,
    "Y"
  )
}


rows_before_chr_filter <-
  nrow(wes)


wes <- wes %>%
  dplyr::filter(
    chr_plain %in% keep_chr
  )


cat(
  "Rows removed by chromosome policy:",
  rows_before_chr_filter -
    nrow(wes),
  "\n"
)


# 13. Verify complete mutation x sample matrices from Step 2

expected_samples <- wes %>%
  dplyr::distinct(
    patient,
    sample_id
  ) %>%
  dplyr::count(
    patient,
    name = "n_expected_samples"
  )


mutation_matrix_check <- wes %>%
  dplyr::count(
    patient,
    mutation_id,
    name = "n_observed_samples"
  ) %>%
  dplyr::left_join(
    expected_samples,
    by = "patient"
  ) %>%
  dplyr::filter(
    n_observed_samples !=
      n_expected_samples
  )


if (nrow(mutation_matrix_check) > 0) {
  
  print(
    head(
      mutation_matrix_check,
      20
    )
  )
  
  stop(
    "Step 2 candidate-somatic data do not form a complete ",
    "mutation x tissue-sample matrix.",
    call. = FALSE
  )
}


duplicate_rows <- wes %>%
  dplyr::count(
    patient,
    mutation_id,
    sample_id
  ) %>%
  dplyr::filter(
    n != 1L
  )


if (nrow(duplicate_rows) > 0) {
  
  print(
    head(
      duplicate_rows,
      20
    )
  )
  
  stop(
    "Duplicated mutation x sample rows detected.",
    call. = FALSE
  )
}


# 14. Load ASCAT summary

ascat_raw <- read_table_auto(
  ASCAT_SUMMARY_PATH
)


required_ascat_summary <- c(
  "sample_id",
  "sex",
  "purity",
  "ploidy",
  "goodness_of_fit",
  "tumour_mapd"
)


missing_ascat_summary <- setdiff(
  required_ascat_summary,
  names(ascat_raw)
)


if (length(missing_ascat_summary) > 0) {
  
  stop(
    "ASCAT summary is missing required column(s): ",
    paste(
      missing_ascat_summary,
      collapse = ", "
    ),
    call. = FALSE
  )
}


ascat_meta_all <- tibble::tibble(
  
  ascat_sample_raw =
    as.character(
      ascat_raw$sample_id
    ),
  
  source_key =
    normalize_id(
      ascat_raw$sample_id
    ),
  
  sample_key =
    apply_sample_recoding(
      ascat_raw$sample_id,
      recoding
    ),
  
  purity =
    suppressWarnings(
      as.numeric(
        ascat_raw$purity
      )
    ),
  
  ploidy =
    suppressWarnings(
      as.numeric(
        ascat_raw$ploidy
      )
    ),
  
  gof =
    suppressWarnings(
      as.numeric(
        ascat_raw$goodness_of_fit
      )
    ),
  
  mapd =
    suppressWarnings(
      as.numeric(
        ascat_raw$tumour_mapd
      )
    ),
  
  sex =
    normalize_sex(
      ascat_raw$sex
    )
)


# 15. Separate ASCAT summary samples not included in WES

ascat_meta_unused <- ascat_meta_all %>%
  dplyr::filter(
    !sample_key %in%
      meta$sample_key
  )


readr::write_tsv(
  ascat_meta_unused,
  file.path(
    OUT_DIR,
    "ascat_summary_samples_not_in_wes.tsv"
  )
)


ascat_meta <- ascat_meta_all %>%
  dplyr::filter(
    sample_key %in%
      meta$sample_key
  )


if (anyDuplicated(
  ascat_meta$sample_key
)) {
  
  duplicated_ascat <- ascat_meta %>%
    dplyr::count(
      sample_key
    ) %>%
    dplyr::filter(
      n > 1L
    )
  
  print(
    duplicated_ascat
  )
  
  stop(
    "Multiple ASCAT summary rows map to the same WES sample.",
    call. = FALSE
  )
}


# 16. Add WES metadata to ASCAT summary

ascat_meta <- ascat_meta %>%
  dplyr::left_join(
    
    meta %>%
      dplyr::select(
        sample_key,
        patient,
        sample_id,
        role
      ),
    
    by = "sample_key"
  )


# 17. Determine patient sex

patient_sex_check <- ascat_meta %>%
  dplyr::filter(
    !is.na(patient),
    !is.na(sex)
  ) %>%
  dplyr::distinct(
    patient,
    sex
  ) %>%
  dplyr::count(
    patient,
    name = "n_sex_values"
  ) %>%
  dplyr::filter(
    n_sex_values > 1L
  )


if (nrow(patient_sex_check) > 0) {
  
  print(
    patient_sex_check
  )
  
  stop(
    "Conflicting sex assignments found within at least one patient.",
    call. = FALSE
  )
}


patient_sex <- ascat_meta %>%
  dplyr::filter(
    !is.na(patient),
    !is.na(sex)
  ) %>%
  dplyr::distinct(
    patient,
    sex
  )


# 18. Find ASCAT segment files

seg_files <- list.files(
  ASCAT_SEG_DIR,
  pattern =
    "\\.segments\\.txt$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)


if (length(seg_files) == 0) {
  
  stop(
    "No *.segments.txt ASCAT files found under:\n",
    ASCAT_SEG_DIR,
    call. = FALSE
  )
}


# 19. Read ASCAT segment files

read_segment_file <- function(path) {
  
  raw_sample_id <-
    basename(path) %>%
    stringr::str_remove(
      "\\.segments\\.txt$"
    )
  
  
  df <- readr::read_tsv(
    path,
    col_types = readr::cols(
      .default = "c"
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  needed <- c(
    "chr",
    "startpos",
    "endpos",
    "nMajor",
    "nMinor"
  )
  
  
  missing <- setdiff(
    needed,
    names(df)
  )
  
  
  if (length(missing) > 0) {
    
    stop(
      "ASCAT segment file has unexpected columns: ",
      basename(path),
      "\nMissing: ",
      paste(
        missing,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  tibble::tibble(
    
    segment_file =
      basename(path),
    
    ascat_sample_raw =
      raw_sample_id,
    
    source_key =
      normalize_id(
        raw_sample_id
      ),
    
    sample_key =
      apply_sample_recoding(
        raw_sample_id,
        recoding
      ),
    
    chr_plain =
      normalize_chr(
        df$chr
      ),
    
    start =
      suppressWarnings(
        as.integer(
          df$startpos
        )
      ),
    
    end =
      suppressWarnings(
        as.integer(
          df$endpos
        )
      ),
    
    nMajor =
      suppressWarnings(
        as.numeric(
          df$nMajor
        )
      ),
    
    nMinor =
      suppressWarnings(
        as.numeric(
          df$nMinor
        )
      )
  )
}


segments_all <- purrr::map_dfr(
  seg_files,
  read_segment_file
)


if (nrow(segments_all) == 0) {
  
  stop(
    "No ASCAT segments could be loaded.",
    call. = FALSE
  )
}


# 20. Separate segment samples not included in WES

segments_unused <- segments_all %>%
  dplyr::filter(
    !sample_key %in%
      meta$sample_key
  )


readr::write_tsv(
  
  segments_unused %>%
    dplyr::distinct(
      ascat_sample_raw,
      sample_key,
      segment_file
    ),
  
  file.path(
    OUT_DIR,
    "ascat_segment_samples_not_in_wes.tsv"
  )
)


segments <- segments_all %>%
  dplyr::filter(
    sample_key %in%
      meta$sample_key
  )


# 21. Validate ASCAT segment coordinates

bad_segment_coordinates <- segments %>%
  dplyr::filter(
    
    is.na(chr_plain) |
      is.na(start) |
      is.na(end) |
      
      start < 1L |
      end < start
  )


if (nrow(bad_segment_coordinates) > 0) {
  
  print(
    head(
      bad_segment_coordinates,
      20
    )
  )
  
  stop(
    "Invalid ASCAT segment coordinates detected.",
    call. = FALSE
  )
}


# 22. Validate ASCAT copy numbers

non_integer_cn <- segments %>%
  dplyr::filter(
    
    (
      !is.na(nMajor) &
        abs(
          nMajor -
            round(nMajor)
        ) > 1e-8
    ) |
      
      (
        !is.na(nMinor) &
          abs(
            nMinor -
              round(nMinor)
          ) > 1e-8
      )
  )


if (nrow(non_integer_cn) > 0) {
  
  print(
    head(
      non_integer_cn,
      20
    )
  )
  
  stop(
    "Non-integer ASCAT major/minor copy-number values detected.",
    call. = FALSE
  )
}


invalid_cn <- segments %>%
  dplyr::filter(
    
    (
      !is.na(nMajor) &
        nMajor < 0
    ) |
      
      (
        !is.na(nMinor) &
          nMinor < 0
      ) |
      
      (
        !is.na(nMajor) &
          !is.na(nMinor) &
          nMinor > nMajor
      )
  )


if (nrow(invalid_cn) > 0) {
  
  print(
    head(
      invalid_cn,
      20
    )
  )
  
  stop(
    "Invalid ASCAT allele-specific copy-number values detected.",
    call. = FALSE
  )
}


segments <- segments %>%
  dplyr::mutate(
    
    nMajor =
      as.integer(
        round(nMajor)
      ),
    
    nMinor =
      as.integer(
        round(nMinor)
      )
  )


# 23. Sample-ID matching QC

wes_keys <- sort(
  unique(
    wes_samples$sample_key
  )
)


ascat_summary_keys <- sort(
  unique(
    ascat_meta$sample_key
  )
)


ascat_segment_keys <- sort(
  unique(
    segments$sample_key
  )
)


summary_matches <- intersect(
  wes_keys,
  ascat_summary_keys
)


segment_matches <- intersect(
  wes_keys,
  ascat_segment_keys
)


unmatched_summary <- setdiff(
  wes_keys,
  ascat_summary_keys
)


unmatched_segments <- setdiff(
  wes_keys,
  ascat_segment_keys
)


cat(
  "\nSample ID matching QC\n"
)


cat(
  "WES tissue samples:        ",
  length(wes_keys),
  "\n"
)


cat(
  "Matched ASCAT summary:     ",
  length(summary_matches),
  "\n"
)


cat(
  "Matched ASCAT segments:    ",
  length(segment_matches),
  "\n"
)


if (length(unmatched_summary) > 0) {
  
  cat(
    "\nWES samples without ASCAT summary:\n"
  )
  
  print(
    unmatched_summary
  )
}


if (length(unmatched_segments) > 0) {
  
  cat(
    "\nWES samples without ASCAT segments:\n"
  )
  
  print(
    unmatched_segments
  )
}


id_matching_qc <- wes_samples %>%
  dplyr::mutate(
    
    has_ascat_summary =
      sample_key %in%
      ascat_summary_keys,
    
    has_ascat_segments =
      sample_key %in%
      ascat_segment_keys
  )


readr::write_tsv(
  id_matching_qc,
  file.path(
    OUT_DIR,
    "ascat_wes_sample_matching_qc.tsv"
  )
)


# 24. Determine which samples can enter PyClone-VI

sample_status <- wes_samples %>%
  dplyr::left_join(
    
    ascat_meta %>%
      dplyr::select(
        sample_key,
        ascat_sample_raw,
        purity,
        ploidy,
        gof,
        mapd
      ),
    
    by = "sample_key"
  ) %>%
  dplyr::mutate(
    
    has_ascat_summary =
      !is.na(
        ascat_sample_raw
      ),
    
    has_segment_file =
      sample_key %in%
      ascat_segment_keys,
    
    purity_valid =
      !is.na(purity) &
      purity > 0 &
      purity <= 1,
    
    preliminary_included =
      dplyr::case_when(
        
        role == "adjacent" ~
          TRUE,
        
        role == "tumor" &
          purity_valid &
          has_segment_file ~
          TRUE,
        
        TRUE ~
          FALSE
      ),
    
    preliminary_reason =
      dplyr::case_when(
        
        role == "adjacent" &
          has_segment_file ~
          "included_adjacent_with_ascat",
        
        role == "adjacent" &
          !has_segment_file ~
          "included_adjacent_diploid_fallback",
        
        role == "tumor" &
          !has_ascat_summary ~
          "excluded_tumor_missing_ascat_summary",
        
        role == "tumor" &
          !purity_valid ~
          "excluded_tumor_missing_or_invalid_purity",
        
        role == "tumor" &
          !has_segment_file ~
          "excluded_tumor_missing_ascat_segments",
        
        role == "tumor" ~
          "included_tumor",
        
        TRUE ~
          "excluded"
      )
  )


# 25. Require at least one usable tumor per patient

patient_status <- sample_status %>%
  dplyr::group_by(
    patient
  ) %>%
  dplyr::summarise(
    
    n_tumor_samples =
      sum(
        role == "tumor"
      ),
    
    n_usable_tumor_samples =
      sum(
        role == "tumor" &
          preliminary_included
      ),
    
    patient_included =
      n_usable_tumor_samples > 0,
    
    .groups = "drop"
  )


sample_status <- sample_status %>%
  dplyr::left_join(
    patient_status,
    by = "patient"
  ) %>%
  dplyr::mutate(
    
    sample_included =
      preliminary_included &
      patient_included,
    
    final_reason =
      dplyr::case_when(
        
        !patient_included ~
          "excluded_patient_no_usable_tumor_ascat",
        
        TRUE ~
          preliminary_reason
      )
  )


readr::write_tsv(
  sample_status,
  file.path(
    OUT_DIR,
    "pyclone_sample_inclusion.tsv"
  )
)


included_samples <- sample_status %>%
  dplyr::filter(
    sample_included
  ) %>%
  dplyr::select(
    patient,
    sample_id,
    sample_key,
    role
  )


if (nrow(included_samples) == 0) {
  
  stop(
    "No samples remain eligible for PyClone-VI.",
    call. = FALSE
  )
}


# 26. Restrict WES data to included samples

wes_included <- wes %>%
  dplyr::semi_join(
    
    included_samples,
    
    by = c(
      "patient",
      "sample_id",
      "sample_key",
      "role"
    )
  )


# 27. Identify mutations with unusable allele counts

bad_count_mutations <- wes_included %>%
  dplyr::filter(
    
    is.na(ref_count) |
      is.na(alt_count) |
      
      ref_count < 0L |
      alt_count < 0L |
      
      ref_count +
      alt_count <= 0L
  ) %>%
  dplyr::distinct(
    patient,
    mutation_id
  ) %>%
  dplyr::mutate(
    exclusion_reason =
      "missing_or_zero_allele_counts"
  )


# 28. Build ASCAT segment GRanges

seg_gr <- GenomicRanges::GRanges(
  
  seqnames =
    segments$chr_plain,
  
  ranges =
    IRanges::IRanges(
      start = segments$start,
      end = segments$end
    ),
  
  sample_key =
    segments$sample_key,
  
  nMajor =
    segments$nMajor,
  
  nMinor =
    segments$nMinor
)


# 29. Build variant GRanges

var_gr <- GenomicRanges::GRanges(
  
  seqnames =
    wes_included$chr_plain,
  
  ranges =
    IRanges::IRanges(
      start = wes_included$pos,
      end = wes_included$pos
    )
)


# 30. Match variants to sample-specific ASCAT segments

hits <- GenomicRanges::findOverlaps(
  var_gr,
  seg_gr
)


hit_df <- tibble::tibble(
  
  query_index =
    S4Vectors::queryHits(
      hits
    ),
  
  segment_index =
    S4Vectors::subjectHits(
      hits
    ),
  
  segment_sample_key =
    S4Vectors::mcols(
      seg_gr
    )$sample_key[
      S4Vectors::subjectHits(
        hits
      )
    ],
  
  nMajor =
    S4Vectors::mcols(
      seg_gr
    )$nMajor[
      S4Vectors::subjectHits(
        hits
      )
    ],
  
  nMinor =
    S4Vectors::mcols(
      seg_gr
    )$nMinor[
      S4Vectors::subjectHits(
        hits
      )
    ]
) %>%
  
  dplyr::mutate(
    
    variant_sample_key =
      wes_included$sample_key[
        query_index
      ]
  ) %>%
  
  dplyr::filter(
    segment_sample_key ==
      variant_sample_key
  )


duplicate_segment_hits <- hit_df %>%
  dplyr::count(
    query_index
  ) %>%
  dplyr::filter(
    n > 1L
  )

if (nrow(duplicate_segment_hits) > 0) {
  
  print(
    head(
      duplicate_segment_hits,
      20
    )
  )
  
  stop(
    "At least one mutation overlaps more than one ASCAT ",
    "segment for the same sample.",
    call. = FALSE
  )
}


wes_included$nMajor <-
  NA_integer_

wes_included$nMinor <-
  NA_integer_


wes_included$nMajor[
  hit_df$query_index
] <- hit_df$nMajor


wes_included$nMinor[
  hit_df$query_index
] <- hit_df$nMinor


# 31. Add ASCAT summary information and patient sex

annot <- wes_included %>%
  dplyr::left_join(
    
    ascat_meta %>%
      dplyr::select(
        sample_key,
        purity,
        ploidy,
        gof,
        mapd
      ),
    
    by = "sample_key"
  ) %>%
  dplyr::left_join(
    patient_sex,
    by = "patient"
  )


# 32. Assign copy number and tumour content

annot <- annot %>%
  dplyr::mutate(
    
    has_usable_segment_cn =
      !is.na(nMajor) &
      !is.na(nMinor),
    
    cn_source =
      dplyr::case_when(
        
        has_usable_segment_cn ~
          "ascat",
        
        role == "adjacent" ~
          "diploid_fallback_adjacent",
        
        TRUE ~
          "missing_tumor_segment"
      ),
    
    major_cn =
      dplyr::case_when(
        
        has_usable_segment_cn ~
          nMajor,
        
        role == "adjacent" ~
          ADJACENT_FALLBACK_MAJOR_CN,
        
        TRUE ~
          NA_integer_
      ),
    
    minor_cn =
      dplyr::case_when(
        
        has_usable_segment_cn ~
          nMinor,
        
        role == "adjacent" ~
          ADJACENT_FALLBACK_MINOR_CN,
        
        TRUE ~
          NA_integer_
      ),
    
    tumour_content =
      dplyr::case_when(
        
        role == "tumor" ~
          purity,
        
        role == "adjacent" ~
          ADJACENT_TUMOUR_CONTENT,
        
        TRUE ~
          NA_real_
      ),
    
    purity_source =
      dplyr::case_when(
        
        role == "tumor" ~
          "ascat",
        
        role == "adjacent" ~
          "whole_adjacent_tissue",
        
        TRUE ~
          NA_character_
      )
  )


# 33. Assign normal copy number

annot <- annot %>%
  dplyr::mutate(
    
    sex_unknown_chrX =
      chr_plain == "X" &
      is.na(sex),
    
    normal_cn =
      dplyr::case_when(
        
        chr_plain == "X" &
          sex == "male" ~
          1L,
        
        chr_plain == "X" ~
          2L,
        
        chr_plain == "Y" ~
          1L,
        
        TRUE ~
          2L
      )
  )


# 34. Identify mutations that cannot be modeled

missing_tumor_cn <- annot %>%
  dplyr::filter(
    
    role == "tumor",
    
    is.na(major_cn) |
      is.na(minor_cn)
  ) %>%
  dplyr::distinct(
    patient,
    mutation_id
  ) %>%
  dplyr::mutate(
    exclusion_reason =
      "missing_tumor_copy_number"
  )


zero_major_cn <- annot %>%
  dplyr::filter(
    
    !is.na(major_cn),
    major_cn < 1L
  ) %>%
  dplyr::distinct(
    patient,
    mutation_id
  ) %>%
  dplyr::mutate(
    exclusion_reason =
      "major_cn_zero"
  )


invalid_final_cn <- annot %>%
  dplyr::filter(
    
    !is.na(major_cn) &
      !is.na(minor_cn) &
      
      (
        minor_cn < 0L |
          major_cn < minor_cn
      )
  )


if (nrow(invalid_final_cn) > 0) {
  
  print(
    head(
      invalid_final_cn,
      20
    )
  )
  
  stop(
    "Invalid final major/minor copy-number configuration detected.",
    call. = FALSE
  )
}


mutation_exclusions <- dplyr::bind_rows(
  bad_count_mutations,
  missing_tumor_cn,
  zero_major_cn
) %>%
  dplyr::distinct(
    patient,
    mutation_id,
    exclusion_reason
  )


readr::write_tsv(
  mutation_exclusions,
  file.path(
    OUT_DIR,
    "pyclone_mutation_exclusions.tsv"
  )
)


# 35. Remove excluded mutations from every included sample

annot_final <- annot %>%
  dplyr::anti_join(
    
    mutation_exclusions %>%
      dplyr::select(
        patient,
        mutation_id
      ) %>%
      dplyr::distinct(),
    
    by = c(
      "patient",
      "mutation_id"
    )
  )


# 36. Verify complete mutation x sample matrices after filtering

final_expected_samples <- annot_final %>%
  dplyr::distinct(
    patient,
    sample_id
  ) %>%
  dplyr::count(
    patient,
    name = "n_expected_samples"
  )


final_matrix_check <- annot_final %>%
  dplyr::count(
    patient,
    mutation_id,
    name = "n_observed_samples"
  ) %>%
  dplyr::left_join(
    final_expected_samples,
    by = "patient"
  ) %>%
  dplyr::filter(
    n_observed_samples !=
      n_expected_samples
  )


if (nrow(final_matrix_check) > 0) {
  
  print(
    head(
      final_matrix_check,
      20
    )
  )
  
  stop(
    "Final PyClone-VI dataset is not a complete ",
    "mutation x sample matrix.",
    call. = FALSE
  )
}


# 37. Final PyClone-VI validation

bad_final_rows <- annot_final %>%
  dplyr::filter(
    
    is.na(ref_count) |
      is.na(alt_count) |
      
      is.na(normal_cn) |
      is.na(major_cn) |
      is.na(minor_cn) |
      
      is.na(tumour_content) |
      
      ref_count < 0L |
      alt_count < 0L |
      
      major_cn < 1L |
      minor_cn < 0L |
      minor_cn > major_cn |
      
      tumour_content <= 0 |
      tumour_content > 1
  )


if (nrow(bad_final_rows) > 0) {
  
  print(
    head(
      bad_final_rows,
      20
    )
  )
  
  stop(
    "Invalid values remain in final PyClone-VI input.",
    call. = FALSE
  )
}


# 38. Identify patients with no retained mutations

patient_mutation_counts <- annot_final %>%
  dplyr::distinct(
    patient,
    mutation_id
  ) %>%
  dplyr::count(
    patient,
    name = "n_mutations"
  )


zero_mutation_patients <- setdiff(
  unique(
    included_samples$patient
  ),
  patient_mutation_counts$patient
)


if (length(zero_mutation_patients) > 0) {
  
  warning(
    "No mutations remained for patient(s): ",
    paste(
      zero_mutation_patients,
      collapse = ", "
    )
  )
}


annot_final <- annot_final %>%
  dplyr::filter(
    patient %in%
      patient_mutation_counts$patient
  )


# 39. Create exact PyClone-VI input

pyclone <- annot_final %>%
  dplyr::transmute(
    
    patient,
    
    mutation_id,
    
    sample_id,
    
    ref_counts =
      ref_count,
    
    alt_counts =
      alt_count,
    
    normal_cn,
    
    major_cn,
    
    minor_cn,
    
    tumour_content
  )


# 40. Write per-patient PyClone-VI files

patients_final <- sort(
  unique(
    pyclone$patient
  )
)


for (patient_id in patients_final) {
  
  patient_input <- pyclone %>%
    dplyr::filter(
      patient == patient_id
    ) %>%
    dplyr::select(
      -patient
    ) %>%
    dplyr::arrange(
      mutation_id,
      sample_id
    )
  
  
  readr::write_tsv(
    patient_input,
    file.path(
      OUT_DIR,
      paste0(
        patient_id,
        "_pyclone_input.tsv"
      )
    )
  )
}


# 41. Write combined cohort input for record keeping

readr::write_tsv(
  pyclone,
  file.path(
    OUT_DIR,
    "ALL_PATIENTS_pyclone_input_with_patient.tsv"
  )
)


# 42. Write full merge provenance

merge_provenance <- annot_final %>%
  dplyr::select(
    
    patient,
    mutation_id,
    sample_id,
    role,
    
    chr,
    pos,
    ref,
    alt,
    
    ref_count,
    alt_count,
    total_depth,
    present,
    
    nMajor,
    nMinor,
    
    major_cn,
    minor_cn,
    normal_cn,
    
    cn_source,
    
    purity,
    tumour_content,
    purity_source,
    
    ploidy,
    gof,
    mapd,
    
    sex,
    sex_unknown_chrX
  )


readr::write_tsv(
  merge_provenance,
  file.path(
    OUT_DIR,
    "pyclone_merge_provenance.tsv"
  )
)


# 43. Sample-level QC

sample_qc <- annot_final %>%
  dplyr::group_by(
    patient,
    sample_id,
    role
  ) %>%
  dplyr::summarise(
    
    n_mutations =
      dplyr::n_distinct(
        mutation_id
      ),
    
    tumour_content =
      dplyr::first(
        tumour_content
      ),
    
    purity_source =
      dplyr::first(
        purity_source
      ),
    
    ascat_purity =
      dplyr::first(
        purity
      ),
    
    ascat_ploidy =
      dplyr::first(
        ploidy
      ),
    
    ascat_gof =
      dplyr::first(
        gof
      ),
    
    ascat_mapd =
      dplyr::first(
        mapd
      ),
    
    pct_cn_from_ascat =
      round(
        100 *
          mean(
            cn_source ==
              "ascat"
          ),
        1
      ),
    
    pct_cn_diploid_fallback =
      round(
        100 *
          mean(
            cn_source ==
              "diploid_fallback_adjacent"
          ),
        1
      ),
    
    median_allelic_depth =
      median(
        ref_count +
          alt_count,
        na.rm = TRUE
      ),
    
    median_vaf =
      median(
        alt_count /
          (
            ref_count +
              alt_count
          ),
        na.rm = TRUE
      ),
    
    n_chrX =
      sum(
        chr_plain == "X"
      ),
    
    n_chrX_unknown_sex =
      sum(
        sex_unknown_chrX
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    patient,
    role,
    sample_id
  )

readr::write_tsv(sample_qc, file.path(OUT_DIR, "pyclone_input_sample_qc.tsv"))

# 44. Patient-level QC

patient_qc <- pyclone %>%
  dplyr::group_by(
    patient
  ) %>%
  dplyr::summarise(
    
    n_samples =
      dplyr::n_distinct(
        sample_id
      ),
    
    n_mutations =
      dplyr::n_distinct(
        mutation_id
      ),
    
    n_rows =
      dplyr::n(),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    patient
  )

readr::write_tsv(patient_qc, file.path(OUT_DIR, "pyclone_input_patient_qc.tsv"))

# 45. Mutation-exclusion QC

mutation_exclusion_qc <-
  mutation_exclusions %>%
  dplyr::count(
    exclusion_reason,
    name = "n_patient_mutations"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      n_patient_mutations
    )
  )

readr::write_tsv(mutation_exclusion_qc, file.path(OUT_DIR, "pyclone_mutation_exclusion_qc.tsv"))

# 46. Analysis settings

analysis_settings <- tibble::tribble(
  
  ~parameter,
  ~value,
  
  "input_variants",
  "Step 2 candidate somatic SNVs",
  
  "additional_minimum_depth_filter",
  "none",
  
  "systematic_ASCAT_ID_normalization",
  "remove _C_ or _N_ token",
  
  "exception_sample_ID_mapping",
  "sample_id_recoding.tsv",
  
  "tumor_copy_number",
  "ASCAT allele-specific copy number",
  
  "tumor_tumour_content",
  "ASCAT purity",
  
  "adjacent_copy_number",
  "ASCAT when available; otherwise diploid 1/1 fallback",
  
  "adjacent_tumour_content",
  as.character(
    ADJACENT_TUMOUR_CONTENT
  ),
  
  "adjacent_fallback_major_cn",
  as.character(
    ADJACENT_FALLBACK_MAJOR_CN
  ),
  
  "adjacent_fallback_minor_cn",
  as.character(
    ADJACENT_FALLBACK_MINOR_CN
  ),
  
  "ASCAT_GOF_exclusion",
  "none",
  
  "ASCAT_MAPD_exclusion",
  "none",
  
  "chrX_retained",
  
  as.character(KEEP_CHRX),
  "chrY_retained", as.character(KEEP_CHRY), "unknown_sex_chrX_normal_cn", "2")

readr::write_tsv(analysis_settings, file.path(OUT_DIR, "pyclone_input_analysis_settings.tsv"))

# 47. Record software versions

writeLines(c("Step 4: ASCAT + WES merge for PyClone-VI", "", paste0("ASCAT summary: ", ASCAT_SUMMARY_PATH),
    paste0("ASCAT segment directory: ", ASCAT_SEG_DIR), "", "R session:", capture.output(sessionInfo())), file.path(OUT_DIR, "pyclone_input_software_versions.txt"))


# 48. Final reports

cat("Step 2 patient files:        ", length(step2_files), "\n")
cat("WES tissue samples:          ", length(wes_keys), "\n")
cat("Matched ASCAT summary:       ", length(summary_matches), "\n")
cat("Matched ASCAT segments:      ", length(segment_matches), "\n")
cat("Patients with PyClone input: ", length(patients_final), "\n")
cat("Samples retained:            ", dplyr::n_distinct(pyclone$sample_id), "\n")
cat("Unique mutations retained:   ", dplyr::n_distinct(paste(pyclone$patient, pyclone$mutation_id)), "\n")
cat("Mutation x sample rows:      ", nrow(pyclone), "\n")
cat("Output directory:            ", OUT_DIR, "\n")

