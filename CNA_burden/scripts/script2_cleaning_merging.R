# script2_cleaning_merging.R

# Merge cell-fraction-weighted sCNA burden data with patient clinical data and prepare the clinical variables used in downstream analyses.

# Input:
#   scna_burden_per_sample.csv
#   Donor_clinical_data.csv

# Output:
#   scna_burden_with_clinical_data.csv


# Packages

library(tidyverse)

# Load data

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

scna_burden <- read_csv(file.path(INPUT_DIR, "scna_burden_per_sample.csv"), show_col_types = FALSE)

clinical <- read_csv(file.path(INPUT_DIR, "Donor_clinical_data.csv"), show_col_types = FALSE)


# Clean patient identifiers

scna_burden <- scna_burden %>%
  mutate(
    Donor = str_remove_all(Donor, "\\s")
  )

clinical <- clinical %>%
  mutate(
    Donor = str_remove_all(Donor, "\\s")
  )


# Check that the clinical table contains one row per patient

duplicated_donors <- clinical %>%
  count(Donor) %>%
  filter(n > 1)

if (nrow(duplicated_donors) > 0) {
  stop(
    "The clinical data contain duplicate Donor IDs. ",
    "Resolve these before merging with the sCNA burden data."
  )
}


# Merge burden and clinical data

data <- scna_burden %>%
  left_join(
    clinical,
    by = "Donor"
  )


# Check for samples without matching clinical data

unmatched_donors <- data %>%
  filter(
    if_all(
      -c(
        sample_id,
        Donor,
        Tissue,
        tissue_class,
        starts_with("burden_"),
        starts_with("cf_weighted_"),
        starts_with("n_"),
        starts_with("chr")
      ),
      is.na
    )
  ) %>%
  distinct(Donor)

if (nrow(unmatched_donors) > 0) {
  warning(
    nrow(unmatched_donors),
    " donor(s) did not have matching clinical information."
  )
}


# Define codes used for missing clinical information

missing_codes <- c(
  "",
  "no data",
  "No data",
  "N/A",
  "n.a.",
  "bd",
  "yes, no data",
  "-"
)

replace_missing_codes <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  
  x[x %in% missing_codes] <- NA_character_
  
  x
}


# Smoking status

# Smoking status is coded as never, former, and current smoker.

data <- data %>%
  mutate(
    Smoking_status = replace_missing_codes(Smoking_status),
    
    smoking_status = case_when(
      Smoking_status == "Non-smoker" ~ "never",
      Smoking_status == "Former-smoker" ~ "former",
      Smoking_status == "Current-Smoker" ~ "current",
      TRUE ~ NA_character_
    ),
    
    smoking_status = factor(
      smoking_status,
      levels = c("never", "former", "current")
    )
  )


# Sex

data <- data %>%
  mutate(
    Gender = replace_missing_codes(Gender),
    
    sex = case_when(
      str_to_lower(Gender) == "female" ~ "female",
      str_to_lower(Gender) == "male" ~ "male",
      TRUE ~ NA_character_
    ),
    
    sex = factor(
      sex,
      levels = c("female", "male")
    )
  )


# Clean age

data <- data %>%
  mutate(
    age_years = parse_number(as.character(Age))
  )


# Create family history of cancer variable
#
# "No family history" is coded as no.
# Any reported family history of cancer is coded as yes.
# Missing or unavailable information remains missing.

data <- data %>%
  mutate(
    Family_History = replace_missing_codes(Family_History),
    
    family_history_cancer = case_when(
      is.na(Family_History) ~ NA_character_,
      str_to_lower(Family_History) == "no family history" ~ "no",
      TRUE ~ "yes"
    ),
    
    family_history_cancer = factor(
      family_history_cancer,
      levels = c("no", "yes")
    )
  )


# Derive pathological T stage and muscle-invasiveness status
#
# According to the Methods:
#   pTa and pT1 -> NMIBC
#   pT2, pT3 and pT4 -> MIBC
#
# pTis is retained as a separate "Tis" category and is not classified as
# either NMIBC or MIBC.

data <- data %>%
  mutate(
    pt_stage = str_extract(
      pTNM_AJCC_2017,
      "pT(?:is|a|[0-4])"
    ),
    
    stage_group = case_when(
      pt_stage == "pTis" ~ "Tis",
      pt_stage %in% c("pTa", "pT1") ~ "NMIBC",
      pt_stage %in% c("pT2", "pT3", "pT4") ~ "MIBC",
      TRUE ~ NA_character_
    ),
    
    stage_group = factor(
      stage_group,
      levels = c("NMIBC", "Tis", "MIBC")
    ),
    
    muscle_invasive = case_when(
      stage_group == "NMIBC" ~ "NMIBC",
      stage_group == "MIBC" ~ "MIBC",
      TRUE ~ NA_character_
    ),
    
    muscle_invasive = factor(
      muscle_invasive,
      levels = c("NMIBC", "MIBC")
    )
  )


# Clean grade while retaining the original categories

data <- data %>%
  mutate(
    Grade = replace_missing_codes(Grade)
  )


# Check the derived clinical variables

clinical_summary <- data %>%
  distinct(
    Donor,
    age_years,
    sex,
    smoking_status,
    family_history_cancer,
    pt_stage,
    stage_group,
    muscle_invasive
  )

summary(clinical_summary$age_years)
table(clinical_summary$sex, useNA = "ifany")
table(clinical_summary$smoking_status, useNA = "ifany")
table(clinical_summary$family_history_cancer, useNA = "ifany")
table(clinical_summary$pt_stage, useNA = "ifany")
table(clinical_summary$stage_group, useNA = "ifany")


# Save

write_csv(data, file.path(OUT_DIR, "scna_burden_with_clinical_data.csv"))