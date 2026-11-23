# step6_analyse_pyclone_results.R

#   Summarizes PyClone-VI mutation-level output into cluster-level, sample-level and patient-level tables for downstream analyses.

# Cluster classification:
#   shared: CP >= 0.10 in at least one tumor and at least one NAT sample

#   tumor_specific: CP >= 0.10 in at least one tumor and < 0.10 in all NAT samples

#   NAT_specific: CP >= 0.10 in at least one NAT and < 0.10 in all tumor samples

#   below_detection: CP < 0.10 in all included samples

#   Tumor truncal: maximum tumor CP > 0.85

# High-confidence cluster:  >5 mutations

# Low-confidence clusters are retained and flagged.


# 1. Packages

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# 2. Paths

RESULTS_DIR <- "path/to/step5_pyclone_output/folder"
META_PATH <- "path/to/sample_metadata.tsv"
OUT_DIR <- "path/to/output/directory"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 3. Analysis settings

CP_DETECTION_THRESHOLD <- 0.10

CP_TRUNCAL_THRESHOLD <- 0.85

HIGH_CONFIDENCE_MIN_MUTATIONS <- 5L


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


safe_max <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      NA_real_
    )
  }
  
  max(x)
}


safe_min <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      NA_real_
    )
  }
  
  min(x)
}


safe_median <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    
    return(
      NA_real_
    )
  }
  
  median(x)
}


safe_fraction <- function(
    numerator,
    denominator
) {
  
  if (
    is.na(denominator) ||
    denominator <= 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  numerator / denominator
}


collapse_ids <- function(x) {
  
  x <- unique(
    x[
      !is.na(x) &
        x != ""
    ]
  )
  
  if (length(x) == 0) {
    
    return(
      NA_character_
    )
  }
  
  paste(
    sort(x),
    collapse = ";"
  )
}


# 5. Validate input paths

if (!dir.exists(RESULTS_DIR)) {
  
  stop(
    "PyClone-VI results directory not found:\n",
    RESULTS_DIR,
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
      stringr::str_trim(
        patient
      ),
    
    sample_id =
      stringr::str_trim(
        sample_id
      ),
    
    sample_key =
      normalize_id(
        sample_id
      ),
    
    role =
      stringr::str_to_lower(
        stringr::str_trim(
          role
        )
      )
  )


allowed_roles <- c(
  "blood",
  "tumor",
  "adjacent"
)


unexpected_roles <- setdiff(
  unique(meta$role),
  allowed_roles
)


if (length(unexpected_roles) > 0) {
  
  stop(
    "Unexpected role(s) in sample metadata: ",
    paste(
      unexpected_roles,
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
    "Duplicated sample IDs detected in metadata.",
    call. = FALSE
  )
}


# 7. Find PyClone-VI result files

result_files <- list.files(
  RESULTS_DIR,
  pattern =
    "_pyclone_results\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)


result_files <- sort(
  result_files
)


if (length(result_files) == 0) {
  
  stop(
    "No *_pyclone_results.tsv files found under:\n",
    RESULTS_DIR,
    call. = FALSE
  )
}


patients_from_files <- basename(
  result_files
) %>%
  stringr::str_remove(
    "_pyclone_results\\.tsv$"
  )


if (anyDuplicated(
  patients_from_files
)) {
  
  stop(
    "More than one PyClone-VI result file was found ",
    "for at least one patient.",
    call. = FALSE
  )
}


# 8. Read PyClone-VI results

read_pyclone_results <- function(path) {
  
  patient_id <-
    basename(path) %>%
    stringr::str_remove(
      "_pyclone_results\\.tsv$"
    )
  
  
  dat <- readr::read_tsv(
    path,
    col_types = readr::cols(
      mutation_id = "c",
      sample_id = "c",
      cluster_id = "i",
      cellular_prevalence = "d",
      cellular_prevalence_std = "d",
      cluster_assignment_prob = "d"
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  required_columns <- c(
    "mutation_id",
    "sample_id",
    "cluster_id",
    "cellular_prevalence",
    "cellular_prevalence_std",
    "cluster_assignment_prob"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      "PyClone-VI result file ",
      basename(path),
      " is missing column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  dat %>%
    dplyr::mutate(
      
      patient =
        patient_id,
      
      mutation_id =
        stringr::str_trim(
          mutation_id
        ),
      
      sample_id =
        stringr::str_trim(
          sample_id
        ),
      
      sample_key =
        normalize_id(
          sample_id
        )
    )
}


all_results <- purrr::map_dfr(
  result_files,
  read_pyclone_results
)


if (nrow(all_results) == 0) {
  
  stop(
    "PyClone-VI result files contain no rows.",
    call. = FALSE
  )
}


# 9. Attach sample roles

all_results <- all_results %>%
  dplyr::left_join(
    
    meta %>%
      dplyr::select(
        metadata_patient = patient,
        sample_key,
        metadata_sample_id = sample_id,
        role
      ),
    
    by = "sample_key"
  )


missing_annotations <- all_results %>%
  dplyr::filter(
    is.na(role)
  ) %>%
  dplyr::distinct(
    patient,
    sample_id,
    sample_key
  )


if (nrow(missing_annotations) > 0) {
  
  print(
    missing_annotations
  )
  
  stop(
    "At least one PyClone-VI sample could not be matched ",
    "to sample_metadata.tsv.",
    call. = FALSE
  )
}


patient_mismatches <- all_results %>%
  dplyr::filter(
    patient != metadata_patient
  ) %>%
  dplyr::distinct(
    patient,
    metadata_patient,
    sample_id
  )


if (nrow(patient_mismatches) > 0) {
  
  print(
    patient_mismatches
  )
  
  stop(
    "Patient identifier mismatch between PyClone-VI filename ",
    "and sample metadata.",
    call. = FALSE
  )
}


unexpected_pyclone_roles <- all_results %>%
  dplyr::filter(
    !role %in% c(
      "tumor",
      "adjacent"
    )
  ) %>%
  dplyr::distinct(
    patient,
    sample_id,
    role
  )


if (nrow(unexpected_pyclone_roles) > 0) {
  
  print(
    unexpected_pyclone_roles
  )
  
  stop(
    "PyClone-VI results contain samples other than tumor/adjacent.",
    call. = FALSE
  )
}


all_results <- all_results %>%
  dplyr::select(
    -metadata_patient,
    -metadata_sample_id
  )


# 10. Validate PyClone-VI numerical results

bad_cp <- all_results %>%
  dplyr::filter(
    
    is.na(cellular_prevalence) |
      cellular_prevalence < 0 |
      cellular_prevalence > 1
  )


if (nrow(bad_cp) > 0) {
  
  print(
    head(
      bad_cp,
      20
    )
  )
  
  stop(
    "Invalid cellular prevalence values detected.",
    call. = FALSE
  )
}


bad_cp_std <- all_results %>%
  dplyr::filter(
    
    !is.na(cellular_prevalence_std) &
      cellular_prevalence_std < 0
  )


if (nrow(bad_cp_std) > 0) {
  
  print(
    head(
      bad_cp_std,
      20
    )
  )
  
  stop(
    "Negative cellular prevalence standard deviations detected.",
    call. = FALSE
  )
}


bad_assignment_prob <- all_results %>%
  dplyr::filter(
    
    !is.na(cluster_assignment_prob) &
      (
        cluster_assignment_prob < 0 |
          cluster_assignment_prob > 1
      )
  )


if (nrow(bad_assignment_prob) > 0) {
  
  print(
    head(
      bad_assignment_prob,
      20
    )
  )
  
  stop(
    "Invalid cluster assignment probabilities detected.",
    call. = FALSE
  )
}


# 11. Verify that each mutation belongs to one cluster

mutation_cluster_check <- all_results %>%
  dplyr::distinct(
    patient,
    mutation_id,
    cluster_id
  ) %>%
  dplyr::count(
    patient,
    mutation_id,
    name = "n_clusters"
  ) %>%
  dplyr::filter(
    n_clusters != 1L
  )


if (nrow(mutation_cluster_check) > 0) {
  
  print(
    head(
      mutation_cluster_check,
      20
    )
  )
  
  stop(
    "At least one mutation is assigned to more than one cluster.",
    call. = FALSE
  )
}


# 12. Write combined mutation-level PyClone-VI results

all_results_out <- all_results %>%
  dplyr::select(
    patient,
    mutation_id,
    sample_id,
    role,
    cluster_id,
    cellular_prevalence,
    cellular_prevalence_std,
    cluster_assignment_prob
  ) %>%
  dplyr::arrange(
    patient,
    cluster_id,
    mutation_id,
    sample_id
  )


readr::write_tsv(all_results_out, file.path(OUT_DIR, "pyclone_all_results.tsv"))


# 13. Sample composition per patient

sample_composition <- all_results %>%
  dplyr::distinct(
    patient,
    sample_id,
    role
  ) %>%
  dplyr::group_by(
    patient
  ) %>%
  dplyr::summarise(
    
    n_tumor_samples =
      sum(
        role == "tumor"
      ),
    
    n_NAT_samples =
      sum(
        role == "adjacent"
      ),
    
    tumor_sample_ids =
      collapse_ids(
        sample_id[
          role == "tumor"
        ]
      ),
    
    NAT_sample_ids =
      collapse_ids(
        sample_id[
          role == "adjacent"
        ]
      ),
    
    multi_tumor =
      n_tumor_samples > 1,
    
    multi_NAT =
      n_NAT_samples > 1,
    
    .groups = "drop"
  )


readr::write_tsv(sample_composition, file.path(OUT_DIR, "pyclone_sample_composition.tsv"))


# 14. Summarise CP for each cluster x sample

cluster_sample_summary <- all_results %>%
  dplyr::group_by(
    patient,
    cluster_id,
    sample_id,
    role
  ) %>%
  dplyr::summarise(
    
    n_cluster_mutations =
      dplyr::n_distinct(
        mutation_id
      ),
    
    cellular_prevalence =
      safe_median(
        cellular_prevalence
      ),
    
    cellular_prevalence_std =
      safe_median(
        cellular_prevalence_std
      ),
    
    cellular_prevalence_min =
      safe_min(
        cellular_prevalence
      ),
    
    cellular_prevalence_max =
      safe_max(
        cellular_prevalence
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    
    detected =
      cellular_prevalence >=
      CP_DETECTION_THRESHOLD
  ) %>%
  dplyr::arrange(
    patient,
    cluster_id,
    role,
    sample_id
  )


readr::write_tsv(cluster_sample_summary, file.path(OUT_DIR, "pyclone_cluster_sample_summary.tsv"))


# 15. Cluster size

cluster_size <- all_results %>%
  dplyr::distinct(
    patient,
    cluster_id,
    mutation_id
  ) %>%
  dplyr::count(
    patient,
    cluster_id,
    name = "n_mutations"
  ) %>%
  dplyr::mutate(
    
    high_confidence_cluster =
      n_mutations >
      HIGH_CONFIDENCE_MIN_MUTATIONS
  )


# 16. Summarise mutation-level cluster assignment probabilities

mutation_assignment_prob <- all_results %>%
  dplyr::group_by(
    patient,
    cluster_id,
    mutation_id
  ) %>%
  dplyr::summarise(
    
    mutation_assignment_prob =
      safe_median(
        cluster_assignment_prob
      ),
    
    .groups = "drop"
  )


cluster_assignment_qc <- mutation_assignment_prob %>%
  dplyr::group_by(
    patient,
    cluster_id
  ) %>%
  dplyr::summarise(
    
    median_assignment_prob =
      safe_median(
        mutation_assignment_prob
      ),
    
    min_assignment_prob =
      safe_min(
        mutation_assignment_prob
      ),
    
    max_assignment_prob =
      safe_max(
        mutation_assignment_prob
      ),
    
    .groups = "drop"
  )


# 17. Summarise tumor CP for each cluster

tumor_cluster_summary <- cluster_sample_summary %>%
  dplyr::filter(
    role == "tumor"
  ) %>%
  dplyr::group_by(
    patient,
    cluster_id
  ) %>%
  dplyr::summarise(
    
    n_tumor_samples =
      dplyr::n(),
    
    max_tumor_cp =
      safe_max(
        cellular_prevalence
      ),
    
    median_tumor_cp =
      safe_median(
        cellular_prevalence
      ),
    
    n_tumor_detected =
      sum(
        cellular_prevalence >=
          CP_DETECTION_THRESHOLD,
        na.rm = TRUE
      ),
    
    tumor_samples_detected =
      collapse_ids(
        sample_id[
          cellular_prevalence >=
            CP_DETECTION_THRESHOLD
        ]
      ),
    
    .groups = "drop"
  )


# 18. Summarise NAT CP for each cluster

nat_cluster_summary <- cluster_sample_summary %>%
  dplyr::filter(
    role == "adjacent"
  ) %>%
  dplyr::group_by(
    patient,
    cluster_id
  ) %>%
  dplyr::summarise(
    
    n_NAT_samples =
      dplyr::n(),
    
    max_NAT_cp =
      safe_max(
        cellular_prevalence
      ),
    
    median_NAT_cp =
      safe_median(
        cellular_prevalence
      ),
    
    n_NAT_detected =
      sum(
        cellular_prevalence >=
          CP_DETECTION_THRESHOLD,
        na.rm = TRUE
      ),
    
    NAT_samples_detected =
      collapse_ids(
        sample_id[
          cellular_prevalence >=
            CP_DETECTION_THRESHOLD
        ]
      ),
    
    .groups = "drop"
  )


# 19. Create cluster-level summary and classification

cluster_summary <- cluster_size %>%
  dplyr::left_join(
    cluster_assignment_qc,
    by = c(
      "patient",
      "cluster_id"
    )
  ) %>%
  dplyr::left_join(
    tumor_cluster_summary,
    by = c(
      "patient",
      "cluster_id"
    )
  ) %>%
  dplyr::left_join(
    nat_cluster_summary,
    by = c(
      "patient",
      "cluster_id"
    )
  ) %>%
  dplyr::mutate(
    
    n_tumor_samples =
      dplyr::coalesce(
        n_tumor_samples,
        0L
      ),
    
    n_NAT_samples =
      dplyr::coalesce(
        n_NAT_samples,
        0L
      ),
    
    n_tumor_detected =
      dplyr::coalesce(
        n_tumor_detected,
        0L
      ),
    
    n_NAT_detected =
      dplyr::coalesce(
        n_NAT_detected,
        0L
      ),
    
    tumor_present =
      n_tumor_detected > 0L,
    
    NAT_present =
      n_NAT_detected > 0L,
    
    classification =
      dplyr::case_when(
        
        tumor_present &
          NAT_present ~
          "shared",
        
        tumor_present &
          !NAT_present ~
          "tumor_specific",
        
        !tumor_present &
          NAT_present ~
          "NAT_specific",
        
        TRUE ~
          "below_detection"
      ),
    
    tumor_truncal =
      !is.na(max_tumor_cp) &
      max_tumor_cp >
      CP_TRUNCAL_THRESHOLD,
    
    shared_truncal =
      classification ==
      "shared" &
      tumor_truncal
  ) %>%
  dplyr::arrange(
    patient,
    cluster_id
  )

readr::write_tsv(cluster_summary, file.path(OUT_DIR, "pyclone_cluster_summary.tsv"))

# 20. Write per-patient cluster tables

PER_PATIENT_DIR <- file.path(OUT_DIR, "per_patient")

dir.create(PER_PATIENT_DIR, recursive = TRUE, showWarnings = FALSE)

for (
  patient_id in
  sort(
    unique(
      cluster_summary$patient
    )
  )
) {
  
  patient_clusters <- cluster_summary %>%
    dplyr::filter(
      patient == patient_id
    )
  
  
  readr::write_tsv(
    patient_clusters,
    file.path(
      PER_PATIENT_DIR,
      paste0(
        patient_id,
        "_cluster_summary.tsv"
      )
    )
  )
}


# 21. Patient-level summary

summarise_patient <- function(
    cluster_df
) {
  
  patient_id <-
    unique(
      cluster_df$patient
    )
  
  
  if (length(patient_id) != 1L) {
    
    stop(
      "Expected one patient in summarise_patient().",
      call. = FALSE
    )
  }
  
  
  get_mutation_count <- function(
    classification_name,
    high_confidence_only = FALSE
  ) {
    
    x <- cluster_df %>%
      dplyr::filter(
        classification ==
          classification_name
      )
    
    
    if (high_confidence_only) {
      
      x <- x %>%
        dplyr::filter(
          high_confidence_cluster
        )
    }
    
    
    sum(
      x$n_mutations,
      na.rm = TRUE
    )
  }
  
  
  get_cluster_count <- function(
    classification_name,
    high_confidence_only = FALSE
  ) {
    
    x <- cluster_df %>%
      dplyr::filter(
        classification ==
          classification_name
      )
    
    
    if (high_confidence_only) {
      
      x <- x %>%
        dplyr::filter(
          high_confidence_cluster
        )
    }
    
    
    nrow(x)
  }
  
  
  n_mut_shared <-
    get_mutation_count(
      "shared"
    )
  
  
  n_mut_tumor_specific <-
    get_mutation_count(
      "tumor_specific"
    )
  
  
  n_mut_NAT_specific <-
    get_mutation_count(
      "NAT_specific"
    )
  
  
  n_mut_below_detection <-
    get_mutation_count(
      "below_detection"
    )
  
  
  n_mut_assigned <-
    n_mut_shared +
    n_mut_tumor_specific +
    n_mut_NAT_specific
  
  
  n_mut_shared_high_conf <-
    get_mutation_count(
      "shared",
      high_confidence_only = TRUE
    )
  
  
  n_mut_tumor_specific_high_conf <-
    get_mutation_count(
      "tumor_specific",
      high_confidence_only = TRUE
    )
  
  
  n_mut_NAT_specific_high_conf <-
    get_mutation_count(
      "NAT_specific",
      high_confidence_only = TRUE
    )
  
  
  n_mut_below_detection_high_conf <-
    get_mutation_count(
      "below_detection",
      high_confidence_only = TRUE
    )
  
  
  n_mut_assigned_high_conf <-
    n_mut_shared_high_conf +
    n_mut_tumor_specific_high_conf +
    n_mut_NAT_specific_high_conf
  
  
  total_mutations <-
    sum(
      cluster_df$n_mutations,
      na.rm = TRUE
    )
  
  
  total_mutations_high_conf <-
    sum(
      cluster_df$n_mutations[
        cluster_df$high_confidence_cluster
      ],
      na.rm = TRUE
    )
  
  
  tibble::tibble(
    
    patient =
      patient_id,
    
    n_tumor_samples =
      max(
        cluster_df$n_tumor_samples,
        na.rm = TRUE
      ),
    
    n_NAT_samples =
      max(
        cluster_df$n_NAT_samples,
        na.rm = TRUE
      ),
    
    n_mutations =
      total_mutations,
    
    n_clusters =
      nrow(
        cluster_df
      ),
    
    n_high_confidence_clusters =
      sum(
        cluster_df$high_confidence_cluster,
        na.rm = TRUE
      ),
    
    n_low_confidence_clusters =
      sum(
        !cluster_df$high_confidence_cluster,
        na.rm = TRUE
      ),
    
    n_shared_clusters =
      get_cluster_count(
        "shared"
      ),
    
    n_tumor_specific_clusters =
      get_cluster_count(
        "tumor_specific"
      ),
    
    n_NAT_specific_clusters =
      get_cluster_count(
        "NAT_specific"
      ),
    
    n_below_detection_clusters =
      get_cluster_count(
        "below_detection"
      ),
    
    n_mut_shared =
      n_mut_shared,
    
    n_mut_tumor_specific =
      n_mut_tumor_specific,
    
    n_mut_NAT_specific =
      n_mut_NAT_specific,
    
    n_mut_below_detection =
      n_mut_below_detection,
    
    n_mut_assigned =
      n_mut_assigned,
    
    frac_shared =
      safe_fraction(
        n_mut_shared,
        n_mut_assigned
      ),
    
    frac_tumor_specific =
      safe_fraction(
        n_mut_tumor_specific,
        n_mut_assigned
      ),
    
    frac_NAT_specific =
      safe_fraction(
        n_mut_NAT_specific,
        n_mut_assigned
      ),
    
    frac_below_detection_all =
      safe_fraction(
        n_mut_below_detection,
        total_mutations
      ),
    
    n_truncal_clusters =
      sum(
        cluster_df$tumor_truncal,
        na.rm = TRUE
      ),
    
    n_shared_truncal_clusters =
      sum(
        cluster_df$shared_truncal,
        na.rm = TRUE
      ),
    
    n_tumor_specific_truncal_clusters =
      sum(
        cluster_df$tumor_truncal &
          cluster_df$classification ==
          "tumor_specific",
        na.rm = TRUE
      ),
    
    n_mut_truncal =
      sum(
        cluster_df$n_mutations[
          cluster_df$tumor_truncal
        ],
        na.rm = TRUE
      ),
    
    n_mut_shared_truncal =
      sum(
        cluster_df$n_mutations[
          cluster_df$shared_truncal
        ],
        na.rm = TRUE
      ),
    
    n_mut_shared_high_conf =
      n_mut_shared_high_conf,
    
    n_mut_tumor_specific_high_conf =
      n_mut_tumor_specific_high_conf,
    
    n_mut_NAT_specific_high_conf =
      n_mut_NAT_specific_high_conf,
    
    n_mut_below_detection_high_conf =
      n_mut_below_detection_high_conf,
    
    n_mut_assigned_high_conf =
      n_mut_assigned_high_conf,
    
    frac_shared_high_conf =
      safe_fraction(
        n_mut_shared_high_conf,
        n_mut_assigned_high_conf
      ),
    
    frac_tumor_specific_high_conf =
      safe_fraction(
        n_mut_tumor_specific_high_conf,
        n_mut_assigned_high_conf
      ),
    
    frac_NAT_specific_high_conf =
      safe_fraction(
        n_mut_NAT_specific_high_conf,
        n_mut_assigned_high_conf
      ),
    
    frac_below_detection_high_conf_all =
      safe_fraction(
        n_mut_below_detection_high_conf,
        total_mutations_high_conf
      )
  )
}


patient_summary <- cluster_summary %>%
  split(
    .$patient
  ) %>%
  purrr::map_dfr(
    summarise_patient
  ) %>%
  dplyr::left_join(
    sample_composition,
    by = c(
      "patient",
      "n_tumor_samples",
      "n_NAT_samples"
    )
  ) %>%
  dplyr::arrange(
    patient
  )


# 22. Validate patient-level fractions

fraction_check <- patient_summary %>%
  dplyr::mutate(
    
    fraction_sum =
      frac_shared +
      frac_tumor_specific +
      frac_NAT_specific,
    
    high_conf_fraction_sum =
      frac_shared_high_conf +
      frac_tumor_specific_high_conf +
      frac_NAT_specific_high_conf
  )


bad_fraction_sum <- fraction_check %>%
  dplyr::filter(
    
    !is.na(fraction_sum) &
      abs(
        fraction_sum - 1
      ) > 1e-8
  )


if (nrow(bad_fraction_sum) > 0) {
  
  print(
    bad_fraction_sum
  )
  
  stop(
    "Patient-level mutation fractions do not sum to 1.",
    call. = FALSE
  )
}


bad_high_conf_fraction_sum <- fraction_check %>%
  dplyr::filter(
    
    !is.na(high_conf_fraction_sum) &
      abs(
        high_conf_fraction_sum - 1
      ) > 1e-8
  )


if (
  nrow(
    bad_high_conf_fraction_sum
  ) > 0
) {
  
  print(
    bad_high_conf_fraction_sum
  )
  
  stop(
    "High-confidence patient-level mutation fractions ",
    "do not sum to 1.",
    call. = FALSE
  )
}


readr::write_tsv(patient_summary, file.path(OUT_DIR, "pyclone_patient_summary.tsv"))


# 23. Cohort-level classification QC

classification_qc <- cluster_summary %>%
  dplyr::group_by(
    classification
  ) %>%
  dplyr::summarise(
    
    n_clusters =
      dplyr::n(),
    
    n_mutations =
      sum(
        n_mutations
      ),
    
    n_high_confidence_clusters =
      sum(
        high_confidence_cluster
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    classification
  )

readr::write_tsv(classification_qc, file.path(OUT_DIR, "pyclone_classification_qc.tsv"))

# 24. Analysis settings

analysis_settings <- tibble::tribble(
  
  ~parameter,
  ~value,
  
  "cluster_detection_threshold",
  paste0(
    "cellular_prevalence >= ",
    CP_DETECTION_THRESHOLD
  ),
  
  "shared_definition",
  "detected in >=1 tumor and >=1 NAT sample",
  
  "tumor_specific_definition",
  "detected in >=1 tumor and no NAT sample",
  
  "NAT_specific_definition",
  "detected in >=1 NAT and no tumor sample",
  
  "below_detection_definition",
  "cellular prevalence below detection threshold in all samples",
  
  "tumor_truncal_definition",
  paste0(
    "maximum tumor cellular prevalence > ",
    CP_TRUNCAL_THRESHOLD
  ),
  
  "high_confidence_cluster_definition",
  paste0(
    "n_mutations > ",
    HIGH_CONFIDENCE_MIN_MUTATIONS
  ),
  
  "low_confidence_clusters_removed",
  "FALSE",
  
  "cluster_assignment_probability_filter",
  "none",
  
  "multiple_tumor_samples",
  "maximum tumor CP used for aggregate classification",
  
  "multiple_NAT_samples",
  "maximum NAT CP used for aggregate classification",
  
  "fraction_denominator",
  "mutations in shared + tumor-specific + NAT-specific clusters"
)

readr::write_tsv(analysis_settings, file.path(OUT_DIR, "pyclone_analysis_settings.tsv"))

# 25. Software versions

writeLines(c("Step 6: PyClone-VI output analysis", "", paste0("PyClone-VI result directory: ", RESULTS_DIR), "", "R session:", capture.output(sessionInfo())), file.path(OUT_DIR, "pyclone_analysis_software_versions.txt"))

# 26. Final report

cat("Patients analysed:              ", dplyr::n_distinct(cluster_summary$patient), "\n")
cat("Samples analysed:               ", dplyr::n_distinct(all_results$sample_id), "\n")
cat("Unique mutations:               ", dplyr::n_distinct(paste(all_results$patient, all_results$mutation_id)), "\n")
cat("Clusters:                       ", nrow(cluster_summary), "\n")
cat("High-confidence clusters:       ", sum(cluster_summary$high_confidence_cluster), "\n")
cat("Shared clusters:                ", sum(cluster_summary$classification == "shared"), "\n")
cat("Tumor-specific clusters:        ", sum(cluster_summary$classification == "tumor_specific"), "\n")
cat("NAT-specific clusters:          ", sum(cluster_summary$classification == "NAT_specific"), "\n")
cat("Below-detection clusters:       ", sum(cluster_summary$classification == "below_detection"), "\n")
cat("Tumor-truncal clusters:         ", sum(cluster_summary$tumor_truncal), "\n")
cat("Output directory:               ", OUT_DIR, "\n")

