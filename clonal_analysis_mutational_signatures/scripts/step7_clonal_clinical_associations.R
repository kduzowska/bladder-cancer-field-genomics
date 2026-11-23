# step7_clonal_clinical_associations.R

# Purpose:
#   Test associations between patient-level clonal architecture measures derived from PyClone-VI and clinical characteristics.

# Primary outcomes:
#   - fraction of mutations assigned to shared clusters
#   - fraction of mutations assigned to NAT-specific clusters

# Primary analyses:
#   - Wilcoxon rank-sum test for two-level categorical variables
#   - Kruskal-Wallis test for categorical variables with >2 levels
#   - Spearman rank correlation for continuous variables
#   - Benjamini-Hochberg correction across all primary clinical tests

# Secondary threshold-based analyses:
#   - shared fraction: <5%, 5-25%, >25%
#   - comparison of low (<5%) versus high (>25%) shared groups
#   - NAT-specific fraction: >25% versus <=25%

# Sampling QC:
#   - association between number of NAT samples and clonal fractions

# 1. Packages

suppressPackageStartupMessages({
  library(tidyverse)
})

# 2. Paths

PYCLONE_SUMMARY <- "path/to/step6/output/pyclone_patient_summary.tsv"
CLINICAL_CSV <-  "path/to/clinical_data.csv"
OUT_DIR <- "path/to/output/directory"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 3. Analysis settings

SHARED_LOW_CUT <- 0.05
SHARED_HIGH_CUT <- 0.25
NAT_SPECIFIC_HIGH_CUT <- 0.25

PRIMARY_OUTCOMES <- c(
  "frac_shared",
  "frac_NAT_specific"
)

CATEGORICAL_VARS <- c(
  "sex",
  "grade",
  "muscle_invasive",
  "pt_stage",
  "smoking_status",
  "preoperative_treatment"
)

CONTINUOUS_VARS <- c(
  "age",
  "pack_years"
)


# 4. Helper functions

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}


safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}


safe_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  IQR(x)
}


group_summary <- function(values, groups) {
  dat <- tibble(value = values, group = groups) %>%
    filter(!is.na(value), !is.na(group))
  
  if (nrow(dat) == 0) return(NA_character_)
  
  dat %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median = median(value),
      iqr = IQR(value),
      .groups = "drop"
    ) %>%
    mutate(text = paste0(group, ": n=", n, ", median=", round(median, 3),
                         ", IQR=", round(iqr, 3))) %>%
    pull(text) %>%
    paste(collapse = "; ")
}


# 5. Validate inputs

if (!file.exists(PYCLONE_SUMMARY)) {
  stop("PyClone patient summary not found: ", PYCLONE_SUMMARY, call. = FALSE)
}

if (!file.exists(CLINICAL_CSV)) {
  stop("Clinical data file not found: ", CLINICAL_CSV, call. = FALSE)
}


# 6. Load PyClone patient summary

cohort <- readr::read_tsv(
  PYCLONE_SUMMARY,
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

required_cohort <- c(
  "patient",
  "frac_shared",
  "frac_NAT_specific",
  "n_NAT_samples",
  "n_shared_clusters"
)

missing_cohort <- setdiff(required_cohort, names(cohort))

if (length(missing_cohort) > 0) {
  stop(
    "PyClone patient summary is missing required column(s): ",
    paste(missing_cohort, collapse = ", "),
    call. = FALSE
  )
}

cohort <- cohort %>%
  mutate(
    patient = str_trim(patient),
    frac_shared = safe_numeric(frac_shared),
    frac_NAT_specific = safe_numeric(frac_NAT_specific),
    n_NAT_samples = as.integer(safe_numeric(n_NAT_samples)),
    n_shared_clusters = as.integer(safe_numeric(n_shared_clusters))
  )

if (anyDuplicated(cohort$patient)) {
  stop("Duplicated patient IDs detected in PyClone patient summary.", call. = FALSE)
}

bad_fractions <- cohort %>%
  filter(
    is.na(frac_shared) | frac_shared < 0 | frac_shared > 1 |
      is.na(frac_NAT_specific) | frac_NAT_specific < 0 | frac_NAT_specific > 1
  )

if (nrow(bad_fractions) > 0) {
  print(bad_fractions)
  stop("Invalid clonal fraction values detected.", call. = FALSE)
}


# 7. Load and clean clinical data

clin_raw <- readr::read_csv(
  CLINICAL_CSV,
  na = c("", "NA", "na", "N/A", "n/a", "-", "?"),
  show_col_types = FALSE
)

required_clinical <- c(
  "Donor",
  "Gender",
  "Age",
  "Grade",
  "pTNM_AJCC_2017",
  "Smoking status",
  "How many cigarettes a day (is/was)",
  "How long did they smoke (years)",
  "Neoadjuvant"
)

missing_clinical <- setdiff(required_clinical, names(clin_raw))

if (length(missing_clinical) > 0) {
  stop(
    "Clinical data are missing required column(s): ",
    paste(missing_clinical, collapse = ", "),
    call. = FALSE
  )
}

clin <- clin_raw %>%
  transmute(
    patient = str_trim(as.character(Donor)),
    sex_raw = str_trim(as.character(Gender)),
    age = safe_numeric(Age),
    grade_raw = str_trim(as.character(Grade)),
    ptnm_raw = str_trim(as.character(pTNM_AJCC_2017)),
    smoking_raw = str_trim(as.character(`Smoking status`)),
    cigs_per_day = safe_numeric(`How many cigarettes a day (is/was)`),
    years_smoked = safe_numeric(`How long did they smoke (years)`),
    preoperative_raw = str_trim(as.character(Neoadjuvant))
  ) %>%
  mutate(
    sex = case_when(
      str_to_lower(sex_raw) == "male" ~ "Male",
      str_to_lower(sex_raw) == "female" ~ "Female",
      TRUE ~ NA_character_
    ),
    
    grade = case_when(
      str_to_lower(grade_raw) == "low" ~ "Low",
      str_to_lower(grade_raw) == "high" ~ "High",
      TRUE ~ NA_character_
    ),
    
    pt_stage_raw = str_extract(
      str_to_lower(ptnm_raw),
      "pt(?:is|a|[0-4][a-c]?)"
    ),
    
    pt_stage = case_when(
      str_detect(pt_stage_raw, "^ptis") ~ "pTis",
      str_detect(pt_stage_raw, "^pta") ~ "pTa",
      str_detect(pt_stage_raw, "^pt1") ~ "pT1",
      str_detect(pt_stage_raw, "^pt2") ~ "pT2",
      str_detect(pt_stage_raw, "^pt3") ~ "pT3",
      str_detect(pt_stage_raw, "^pt4") ~ "pT4",
      TRUE ~ NA_character_
    ),
    
    muscle_invasive = case_when(
      pt_stage %in% c("pTis", "pTa", "pT1") ~ "NMIBC",
      pt_stage %in% c("pT2", "pT3", "pT4") ~ "MIBC",
      TRUE ~ NA_character_
    ),
    
    smoking_status = case_when(
      str_detect(str_to_lower(smoking_raw), "^non|never") ~ "Never",
      str_detect(str_to_lower(smoking_raw), "former|past") ~ "Former",
      str_detect(str_to_lower(smoking_raw), "current|present") ~ "Current",
      TRUE ~ NA_character_
    ),
    
    pack_years = case_when(
      smoking_status == "Never" ~ 0,
      smoking_status %in% c("Former", "Current") &
        !is.na(cigs_per_day) & !is.na(years_smoked) ~
        (cigs_per_day / 20) * years_smoked,
      TRUE ~ NA_real_
    ),
    
    preoperative_treatment = case_when(
      str_to_lower(preoperative_raw) == "yes" ~ "Yes",
      str_to_lower(preoperative_raw) == "no" ~ "No",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    grade = factor(grade, levels = c("Low", "High")),
    muscle_invasive = factor(muscle_invasive, levels = c("NMIBC", "MIBC")),
    pt_stage = factor(pt_stage, levels = c("pTis", "pTa", "pT1", "pT2", "pT3", "pT4")),
    smoking_status = factor(smoking_status, levels = c("Never", "Former", "Current")),
    preoperative_treatment = factor(preoperative_treatment, levels = c("No", "Yes"))
  ) %>%
  select(
    patient, sex, age, grade, pt_stage, muscle_invasive,
    smoking_status, pack_years, preoperative_treatment
  )

if (anyDuplicated(clin$patient)) {
  duplicated_clin <- clin %>% count(patient) %>% filter(n > 1)
  print(duplicated_clin)
  stop("Duplicated donor IDs detected in clinical_data.csv.", call. = FALSE)
}


# 8. Merge PyClone and clinical data

unmatched_pyclone <- setdiff(cohort$patient, clin$patient)

if (length(unmatched_pyclone) > 0) {
  warning(
    "PyClone patients without matching clinical data: ",
    paste(unmatched_pyclone, collapse = ", ")
  )
}

dat <- cohort %>%
  inner_join(clin, by = "patient") %>%
  mutate(
    shared_group = case_when(
      frac_shared < SHARED_LOW_CUT ~ "low",
      frac_shared > SHARED_HIGH_CUT ~ "high",
      TRUE ~ "intermediate"
    ),
    
    shared_group = factor(
      shared_group,
      levels = c("low", "intermediate", "high")
    ),
    
    shared_extreme_group = case_when(
      shared_group == "low" ~ "low",
      shared_group == "high" ~ "high",
      TRUE ~ NA_character_
    ),
    
    shared_extreme_group = factor(
      shared_extreme_group,
      levels = c("low", "high")
    ),
    
    NAT_specific_group = if_else(
      frac_NAT_specific > NAT_SPECIFIC_HIGH_CUT,
      "high",
      "not_high"
    ),
    
    NAT_specific_group = factor(
      NAT_specific_group,
      levels = c("not_high", "high")
    ),
    
    no_shared_clusters = n_shared_clusters == 0L
  )

write_tsv(dat, file.path(OUT_DIR, "clonal_clinical_merged.tsv"))

# 9. Clinical-variable QC

clinical_qc_cat <- map_dfr(CATEGORICAL_VARS, function(v) {
  x <- dat[[v]]
  counts <- table(x, useNA = "no")
  
  tibble(
    variable = v,
    type = "categorical",
    n_nonmissing = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    summary = paste(names(counts), as.integer(counts), sep = "=", collapse = "; ")
  )
})

clinical_qc_cont <- map_dfr(CONTINUOUS_VARS, function(v) {
  x <- dat[[v]]
  
  tibble(
    variable = v,
    type = "continuous",
    n_nonmissing = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    summary = if (sum(!is.na(x)) > 0) {
      paste0(
        "median=", round(median(x, na.rm = TRUE), 2),
        "; IQR=", round(IQR(x, na.rm = TRUE), 2),
        "; range=", round(min(x, na.rm = TRUE), 2),
        "-", round(max(x, na.rm = TRUE), 2)
      )
    } else {
      NA_character_
    }
  )
})

clinical_variable_qc <- bind_rows(clinical_qc_cat, clinical_qc_cont)

write_tsv(
  clinical_variable_qc,
  file.path(OUT_DIR, "clinical_variable_qc.tsv")
)


# 10. Primary analyses: continuous clonal fractions vs clinical variables

test_primary_categorical <- function(d, outcome, variable) {
  sub <- d %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[variable]]))
  
  grp <- droplevels(sub[[variable]])
  n_groups <- nlevels(grp)
  
  if (n_groups < 2) {
    return(tibble(
      outcome = outcome,
      clinical_variable = variable,
      test = "skipped",
      n = nrow(sub),
      statistic = NA_real_,
      p_value = NA_real_,
      detail = "fewer than two observed groups"
    ))
  }
  
  detail <- group_summary(sub[[outcome]], grp)
  
  if (n_groups == 2) {
    fit <- suppressWarnings(
      wilcox.test(sub[[outcome]] ~ grp, exact = FALSE)
    )
    
    tibble(
      outcome = outcome,
      clinical_variable = variable,
      test = "Wilcoxon rank-sum",
      n = nrow(sub),
      statistic = unname(fit$statistic),
      p_value = fit$p.value,
      detail = detail
    )
  } else {
    fit <- kruskal.test(sub[[outcome]] ~ grp)
    
    tibble(
      outcome = outcome,
      clinical_variable = variable,
      test = "Kruskal-Wallis",
      n = nrow(sub),
      statistic = unname(fit$statistic),
      p_value = fit$p.value,
      detail = detail
    )
  }
}


test_primary_continuous <- function(d, outcome, variable) {
  sub <- d %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[variable]]))
  
  if (
    nrow(sub) < 3 ||
    n_distinct(sub[[outcome]]) < 2 ||
    n_distinct(sub[[variable]]) < 2
  ) {
    return(tibble(
      outcome = outcome,
      clinical_variable = variable,
      test = "skipped",
      n = nrow(sub),
      statistic = NA_real_,
      p_value = NA_real_,
      detail = "insufficient variation or fewer than three complete observations"
    ))
  }
  
  fit <- suppressWarnings(
    cor.test(
      sub[[variable]],
      sub[[outcome]],
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    outcome = outcome,
    clinical_variable = variable,
    test = "Spearman",
    n = nrow(sub),
    statistic = unname(fit$estimate),
    p_value = fit$p.value,
    detail = paste0("rho=", round(unname(fit$estimate), 3))
  )
}


primary_results <- bind_rows(
  tidyr::crossing(
    outcome = PRIMARY_OUTCOMES,
    clinical_variable = CATEGORICAL_VARS
  ) %>%
    pmap_dfr(
      ~ test_primary_categorical(
        dat,
        outcome = ..1,
        variable = ..2
      )
    ),
  
  tidyr::crossing(
    outcome = PRIMARY_OUTCOMES,
    clinical_variable = CONTINUOUS_VARS
  ) %>%
    pmap_dfr(
      ~ test_primary_continuous(
        dat,
        outcome = ..1,
        variable = ..2
      )
    )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significant_BH = !is.na(p_adj_BH) & p_adj_BH < 0.05
  ) %>%
  arrange(p_adj_BH, p_value)

write_tsv(
  primary_results,
  file.path(OUT_DIR, "clonal_clinical_primary_tests.tsv")
)


# 11. Shared-fraction groups

shared_group_summary <- dat %>%
  count(shared_group, name = "n_patients") %>%
  mutate(
    percent_patients = 100 * n_patients / nrow(dat)
  ) %>%
  left_join(
    dat %>%
      group_by(shared_group) %>%
      summarise(
        n_no_shared_clusters = sum(no_shared_clusters, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "shared_group"
  )

write_tsv(
  shared_group_summary,
  file.path(OUT_DIR, "shared_fraction_group_summary.tsv")
)

# 12. Helper for threshold-based analyses

test_binary_group <- function(d, exposure, variable, variable_type) {
  sub <- d %>%
    filter(!is.na(.data[[exposure]]), !is.na(.data[[variable]]))
  
  exposure_factor <- droplevels(sub[[exposure]])
  
  if (nlevels(exposure_factor) < 2) {
    return(tibble(
      exposure = exposure,
      clinical_variable = variable,
      test = "skipped",
      n = nrow(sub),
      statistic = NA_real_,
      p_value = NA_real_,
      detail = "fewer than two exposure groups"
    ))
  }
  
  if (variable_type == "categorical") {
    clinical_factor <- droplevels(sub[[variable]])
    
    if (nlevels(clinical_factor) < 2) {
      return(tibble(
        exposure = exposure,
        clinical_variable = variable,
        test = "skipped",
        n = nrow(sub),
        statistic = NA_real_,
        p_value = NA_real_,
        detail = "fewer than two clinical groups"
      ))
    }
    
    tab <- table(exposure_factor, clinical_factor)
    fit <- fisher.test(tab)
    
    return(tibble(
      exposure = exposure,
      clinical_variable = variable,
      test = "Fisher exact",
      n = nrow(sub),
      statistic = NA_real_,
      p_value = fit$p.value,
      detail = paste(capture.output(print(tab)), collapse = " | ")
    ))
  }
  
  fit <- suppressWarnings(
    wilcox.test(sub[[variable]] ~ exposure_factor, exact = FALSE)
  )
  
  tibble(
    exposure = exposure,
    clinical_variable = variable,
    test = "Wilcoxon rank-sum",
    n = nrow(sub),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    detail = group_summary(sub[[variable]], exposure_factor)
  )
}


# 13. Shared fraction: low (<5%) versus high (>25%)

shared_extreme <- dat %>%
  filter(!is.na(shared_extreme_group))

shared_low_high_results <- bind_rows(
  map_dfr(
    CATEGORICAL_VARS,
    ~ test_binary_group(
      shared_extreme,
      exposure = "shared_extreme_group",
      variable = .x,
      variable_type = "categorical"
    )
  ),
  
  map_dfr(
    CONTINUOUS_VARS,
    ~ test_binary_group(
      shared_extreme,
      exposure = "shared_extreme_group",
      variable = .x,
      variable_type = "continuous"
    )
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

write_tsv(
  shared_low_high_results,
  file.path(OUT_DIR, "shared_low_vs_high_clinical_tests.tsv")
)


# 14. NAT-specific fraction: >25% versus <=25%

nat_specific_group_summary <- dat %>%
  count(NAT_specific_group, name = "n_patients") %>%
  mutate(
    percent_patients = 100 * n_patients / nrow(dat)
  )

write_tsv(
  nat_specific_group_summary,
  file.path(OUT_DIR, "NAT_specific_high_group_summary.tsv")
)

nat_specific_results <- bind_rows(
  map_dfr(
    CATEGORICAL_VARS,
    ~ test_binary_group(
      dat,
      exposure = "NAT_specific_group",
      variable = .x,
      variable_type = "categorical"
    )
  ),
  
  map_dfr(
    CONTINUOUS_VARS,
    ~ test_binary_group(
      dat,
      exposure = "NAT_specific_group",
      variable = .x,
      variable_type = "continuous"
    )
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

write_tsv(
  nat_specific_results,
  file.path(OUT_DIR, "NAT_specific_high_clinical_tests.tsv")
)


# 15. NAT sample-count QC

test_nat_count_qc <- function(outcome) {
  sub <- dat %>%
    filter(!is.na(n_NAT_samples), !is.na(.data[[outcome]]))
  
  if (
    nrow(sub) < 3 ||
    n_distinct(sub$n_NAT_samples) < 2 ||
    n_distinct(sub[[outcome]]) < 2
  ) {
    return(tibble(
      outcome = outcome,
      n = nrow(sub),
      spearman_rho = NA_real_,
      p_value = NA_real_,
      note = "insufficient variation"
    ))
  }
  
  fit <- suppressWarnings(
    cor.test(
      sub$n_NAT_samples,
      sub[[outcome]],
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    outcome = outcome,
    n = nrow(sub),
    spearman_rho = unname(fit$estimate),
    p_value = fit$p.value,
    note = NA_character_
  )
}

nat_sample_count_qc <- map_dfr(
  PRIMARY_OUTCOMES,
  test_nat_count_qc
)

write_tsv(
  nat_sample_count_qc,
  file.path(OUT_DIR, "NAT_sample_count_qc.tsv")
)


# 16. Analysis settings

analysis_settings <- tribble(
  ~parameter, ~value,
  
  "primary_clonal_outcomes",
  "frac_shared; frac_NAT_specific",
  
  "clonal_fraction_denominator",
  "mutations assigned to shared + tumor-specific + NAT-specific clusters; below-detection clusters excluded",
  
  "primary_two_group_test",
  "Wilcoxon rank-sum",
  
  "primary_multi_group_test",
  "Kruskal-Wallis",
  
  "primary_continuous_test",
  "Spearman rank correlation",
  
  "primary_multiple_testing",
  "Benjamini-Hochberg correction across all primary clinical association tests",
  
  "primary_significance_threshold",
  "BH-adjusted p < 0.05",
  
  "shared_low_definition",
  "frac_shared < 0.05",
  
  "shared_intermediate_definition",
  "0.05 <= frac_shared <= 0.25",
  
  "shared_high_definition",
  "frac_shared > 0.25",
  
  "NAT_specific_high_definition",
  "frac_NAT_specific > 0.25",
  
  "threshold_categorical_test",
  "Fisher exact",
  
  "threshold_continuous_test",
  "Wilcoxon rank-sum",
  
  "threshold_multiple_testing",
  "Benjamini-Hochberg correction within each secondary threshold-analysis family",
  
  "NAT_sample_count_QC",
  "Spearman correlation with frac_shared and frac_NAT_specific"
)

write_tsv(analysis_settings, file.path(OUT_DIR, "clonal_clinical_analysis_settings.tsv"))

# 17. Software versions

writeLines(
  c(
    "Step 7: Clonal architecture versus clinical characteristics",
    "",
    paste0("PyClone summary: ", PYCLONE_SUMMARY),
    paste0("Clinical data: ", CLINICAL_CSV),
    "",
    "R session:",
    capture.output(sessionInfo())
  ),
  file.path(OUT_DIR, "clonal_clinical_software_versions.txt")
)

# 18. Final report

cat("Patients analysed:", nrow(dat), "\n")
cat("Primary clinical tests:", nrow(primary_results), "\n")
cat("Shared low/high patients:", nrow(shared_extreme), "\n")
cat("Patients with >25% NAT-specific variants:", sum(dat$NAT_specific_group == "high"), "\n")
cat("Output directory:", OUT_DIR, "\n")

