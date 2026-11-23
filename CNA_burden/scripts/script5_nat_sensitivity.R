# 05_nat_patient_level_sensitivity.R

# Patient-level sensitivity analysis of associations between cell-fraction-weighted sCNA burden in normal-adjacent tissue (NAT) samples and clinical characteristics.

# For each patient, burden values are averaged across NAT samples.
# Patient-level mean burdens are converted to ranks and modeled using
# linear regression as a function of:
#   - age per 10-year increase
#   - sex
#   - smoking status (never, former, current)
#   - family history of cancer
#
# P-values are obtained from partial F-tests and adjusted across the 16 burden-clinical variable combinations using the Benjamini-Hochberg method.

# Input:
# scna_burden_with_clinical_data.csv

# Output:
# nat_patient_level_sensitivity_tests.csv
# nat_patient_level_sensitivity_coefficients.csv


# Packages

library(tidyverse)

# Load data

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

data <- read_csv(file.path(INPUT_DIR, "scna_burden_with_clinical_data.csv"), show_col_types = FALSE)


# Select NAT samples

nat <- data %>%
  filter(tissue_class == "non-cancerous")


# Define genome-wide burden measures

burden_variables <- c(
  "burden_total_pct",
  "burden_gain_pct",
  "burden_loss_pct",
  "burden_cnloh_pct"
)


# Check that burden values are available for all NAT samples

if (anyNA(nat[burden_variables])) {
  stop(
    "Missing sCNA burden values were detected in NAT samples."
  )
}


# Check that patient-level clinical information is consistent
# across NAT samples from the same patient

clinical_consistency <- nat %>%
  group_by(Donor) %>%
  summarise(
    n_age = n_distinct(age_years[!is.na(age_years)]),
    n_sex = n_distinct(sex[!is.na(sex)]),
    n_smoking = n_distinct(
      smoking_status[!is.na(smoking_status)]
    ),
    n_family_history = n_distinct(
      family_history_cancer[
        !is.na(family_history_cancer)
      ]
    ),
    .groups = "drop"
  )

inconsistent_patients <- clinical_consistency %>%
  filter(
    n_age > 1 |
      n_sex > 1 |
      n_smoking > 1 |
      n_family_history > 1
  )

if (nrow(inconsistent_patients) > 0) {
  stop(
    "Patient-level clinical variables are inconsistent across ",
    "NAT samples for at least one patient."
  )
}


# Average burden values across NAT samples for each patient

patient_level <- nat %>%
  group_by(Donor) %>%
  summarise(
    burden_total_pct = mean(burden_total_pct),
    burden_gain_pct = mean(burden_gain_pct),
    burden_loss_pct = mean(burden_loss_pct),
    burden_cnloh_pct = mean(burden_cnloh_pct),
    
    age_years = first(
      age_years[!is.na(age_years)],
      default = NA_real_
    ),
    
    sex = first(
      as.character(sex)[!is.na(sex)],
      default = NA_character_
    ),
    
    smoking_status = first(
      as.character(smoking_status)[!is.na(smoking_status)],
      default = NA_character_
    ),
    
    family_history_cancer = first(
      as.character(family_history_cancer)[
        !is.na(family_history_cancer)
      ],
      default = NA_character_
    ),
    
    n_NATs = n(),
    
    .groups = "drop"
  )


# Prepare clinical variables

patient_level <- patient_level %>%
  mutate(
    age_10y = age_years / 10,
    
    sex = factor(
      sex,
      levels = c("female", "male")
    ),
    
    smoking_status = factor(
      smoking_status,
      levels = c("never", "former", "current")
    ),
    
    family_history_cancer = factor(
      family_history_cancer,
      levels = c("no", "yes")
    )
  )


# Retain patients with complete clinical information

patient_level <- patient_level %>%
  filter(
    !is.na(age_10y),
    !is.na(sex),
    !is.na(smoking_status),
    !is.na(family_history_cancer)
  ) %>%
  droplevels()


cat(
  "Patients included in sensitivity analysis:",
  nrow(patient_level),
  "\n"
)

cat(
  "Median number of NAT samples per patient:",
  median(patient_level$n_NATs),
  "\n"
)


# Check that the clinical variables contain sufficient variation

if (n_distinct(patient_level$age_10y) < 2) {
  stop("Age does not vary among the included patients.")
}

if (nlevels(patient_level$sex) < 2) {
  stop("Sex has fewer than two observed groups.")
}

if (nlevels(patient_level$smoking_status) < 2) {
  stop("Smoking status has fewer than two observed groups.")
}

if (nlevels(patient_level$family_history_cancer) < 2) {
  stop(
    "Family history of cancer has fewer than two observed groups."
  )
}


# Convert patient-level mean burdens to ranks
#
# Average ranks are assigned to tied values, including tied zero burdens.

patient_level <- patient_level %>%
  mutate(
    rank_total = rank(
      burden_total_pct,
      ties.method = "average"
    ),
    
    rank_gain = rank(
      burden_gain_pct,
      ties.method = "average"
    ),
    
    rank_loss = rank(
      burden_loss_pct,
      ties.method = "average"
    ),
    
    rank_cnloh = rank(
      burden_cnloh_pct,
      ties.method = "average"
    )
  )


# Define ranked outcomes

rank_outcomes <- tribble(
  ~rank_variable, ~burden,
  "rank_total",   "burden_total_pct",
  "rank_gain",    "burden_gain_pct",
  "rank_loss",    "burden_loss_pct",
  "rank_cnloh",   "burden_cnloh_pct"
)


# Fit one linear regression model for each ranked burden

fit_patient_model <- function(rank_variable) {
  
  formula <- as.formula(
    paste0(
      rank_variable,
      " ~ age_10y + sex + smoking_status + ",
      "family_history_cancer"
    )
  )
  
  lm(
    formula,
    data = patient_level
  )
}


rank_outcomes <- rank_outcomes %>%
  mutate(
    model = map(
      rank_variable,
      fit_patient_model
    )
  )


# Extract one partial F-test for each clinical variable
#
# drop1() compares the full model with a model excluding one clinical
# variable at a time. Smoking therefore receives one overall test,
# regardless of the number of smoking-status coefficients.

extract_model_tests <- function(
    model,
    rank_variable,
    burden
) {
  
  test_table <- drop1(
    model,
    test = "F"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    filter(variable != "<none>")
  
  test_table %>%
    transmute(
      burden = burden,
      rank_variable = rank_variable,
      variable = variable,
      numerator_df = Df,
      denominator_df = df.residual(model),
      F_value = `F value`,
      p_value = `Pr(>F)`
    )
}


sensitivity_tests <- pmap_dfr(
  rank_outcomes,
  extract_model_tests
)


# Confirm that four burdens x four clinical variables were tested

if (nrow(sensitivity_tests) != 16) {
  stop(
    "Expected 16 burden-clinical variable tests, but obtained ",
    nrow(sensitivity_tests),
    "."
  )
}


# Adjust P-values across all 16 tests

sensitivity_tests <- sensitivity_tests %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significant = q_value < 0.05,
    
    clinical_characteristic = case_when(
      variable == "age_10y" ~
        "Age (per 10-year increase)",
      
      variable == "sex" ~
        "Sex",
      
      variable == "smoking_status" ~
        "Smoking status",
      
      variable == "family_history_cancer" ~
        "Family history of cancer",
      
      TRUE ~ variable
    )
  ) %>%
  select(
    burden,
    rank_variable,
    clinical_characteristic,
    variable,
    numerator_df,
    denominator_df,
    F_value,
    p_value,
    q_value,
    significant
  ) %>%
  arrange(q_value)


# Extract regression coefficients and 95% confidence intervals

extract_coefficients <- function(
    model,
    rank_variable,
    burden
) {
  
  estimates <- coef(model)
  
  confidence_intervals <- confint(
    model,
    level = 0.95
  )
  
  tibble(
    term = names(estimates),
    estimate = unname(estimates),
    CI_low = confidence_intervals[, 1],
    CI_high = confidence_intervals[, 2]
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        term == "age_10y" ~
          "age_10y",
        
        str_starts(term, "sex") ~
          "sex",
        
        str_starts(term, "smoking_status") ~
          "smoking_status",
        
        str_starts(
          term,
          "family_history_cancer"
        ) ~
          "family_history_cancer",
        
        TRUE ~ NA_character_
      ),
      
      burden = burden,
      rank_variable = rank_variable
    ) %>%
    select(
      burden,
      rank_variable,
      variable,
      term,
      estimate,
      CI_low,
      CI_high
    )
}


sensitivity_coefficients <- pmap_dfr(
  rank_outcomes,
  extract_coefficients
)


# Add the omnibus P-value and adjusted P-value for each clinical variable

sensitivity_coefficients <- sensitivity_coefficients %>%
  left_join(
    sensitivity_tests %>%
      select(
        burden,
        variable,
        p_value,
        q_value,
        significant
      ),
    by = c(
      "burden",
      "variable"
    )
  )


# Display results

sensitivity_tests %>%
  mutate(
    across(
      c(
        numerator_df,
        denominator_df,
        F_value
      ),
      ~ round(.x, 3)
    ),
    
    across(
      c(
        p_value,
        q_value
      ),
      ~ signif(.x, 3)
    )) %>%
  print(n = Inf)


# Save results

write_csv(sensitivity_tests, file.path(OUT_DIR, "nat_patient_level_sensitivity_tests.csv"))
write_csv(sensitivity_coefficients, file.path(OUT_DIR, "nat_patient_level_sensitivity_coefficients.csv"))