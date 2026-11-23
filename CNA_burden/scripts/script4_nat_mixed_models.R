# script4_nat_mixed_models.R

# Test associations between cell-fraction-weighted sCNA burden in normal-adjacent tissue (NAT) samples and clinical characteristics.

# For each genome-wide burden measure (total, gain, loss, and CN-LOH), log(1 + burden) is modeled as a function of:
#   - age per 10-year increase
#   - sex
#   - smoking status (never, former, current)
#   - family history of cancer

# A patient-specific random intercept accounts for multiple NAT samples collected from the same patient.

# P-values for clinical variables are obtained using F-tests with Satterthwaite's approximation. P-values are adjusted across the
# 16 burden-clinical variable combinations using the Benjamini-Hochberg method.

# Input:
# scna_burden_with_clinical_data.csv

# Output:
# nat_mixed_model_tests.csv
# nat_mixed_model_coefficients.csv


# Packages

library(tidyverse)
library(lme4)
library(lmerTest)


# Load data

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

data <- read_csv(file.path(INPUT_DIR,"scna_burden_with_clinical_data.csv"), show_col_types = FALSE)


# Select NAT samples and variables required for the analysis

nat_all <- data %>%
  filter(tissue_class == "non-cancerous") %>%
  transmute(
    sample_id = sample_id,
    Donor = Donor,
    burden_total_pct = burden_total_pct,
    burden_gain_pct = burden_gain_pct,
    burden_loss_pct = burden_loss_pct,
    burden_cnloh_pct = burden_cnloh_pct,
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


# Check that burden values are available for all NAT samples

burden_variables <- c(
  "burden_total_pct",
  "burden_gain_pct",
  "burden_loss_pct",
  "burden_cnloh_pct"
)

if (anyNA(nat_all[burden_variables])) {
  stop(
    "Missing sCNA burden values were detected in NAT samples."
  )
}


# Retain NAT samples with complete clinical information

nat <- nat_all %>%
  filter(
    !is.na(age_10y),
    !is.na(sex),
    !is.na(smoking_status),
    !is.na(family_history_cancer)
  ) %>%
  droplevels()


cat(
  "NAT samples before complete-case filtering:",
  nrow(nat_all),
  "\n"
)

cat(
  "NAT samples included in mixed-effects models:",
  nrow(nat),
  "\n"
)

cat(
  "Patients included in mixed-effects models:",
  n_distinct(nat$Donor),
  "\n"
)


# Check that patient-level clinical variables are consistent across NATs
# from the same patient

clinical_consistency <- nat %>%
  group_by(Donor) %>%
  summarise(
    n_age = n_distinct(age_10y),
    n_sex = n_distinct(sex),
    n_smoking = n_distinct(smoking_status),
    n_family_history = n_distinct(family_history_cancer),
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
    "Patient-level clinical variables are inconsistent across NAT samples ",
    "for at least one patient."
  )
}


# Check that all clinical variables contain sufficient variation

if (n_distinct(nat$age_10y) < 2) {
  stop("Age does not vary among the included patients.")
}

if (nlevels(nat$sex) < 2) {
  stop("Sex has fewer than two observed groups.")
}

if (nlevels(nat$smoking_status) < 3) {
  stop(
    "Smoking status does not contain all three groups: ",
    "never, former, and current."
  )
}

if (nlevels(nat$family_history_cancer) < 2) {
  stop("Family history of cancer has fewer than two observed groups.")
}


# Fit one mixed-effects model for each genome-wide burden measure

fit_nat_model <- function(burden) {
  
  formula <- as.formula(
    paste0(
      "log1p(",
      burden,
      ") ~ age_10y + sex + smoking_status + ",
      "family_history_cancer + (1 | Donor)"
    )
  )
  
  lmerTest::lmer(
    formula,
    data = nat
  )
}


models <- burden_variables %>%
  set_names() %>%
  map(fit_nat_model)


# Extract omnibus F-tests for each clinical variable
#
# These are the statistical tests described in the Methods.
# Smoking status receives one overall p-value despite having three groups.

extract_model_tests <- function(model, burden) {
  
  anova_table <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("variable")
  
  anova_table %>%
    transmute(
      burden = burden,
      variable = variable,
      numerator_df = NumDF,
      denominator_df = DenDF,
      F_value = `F value`,
      p_value = `Pr(>F)`
    )
}


nat_tests <- map2_dfr(
  models,
  names(models),
  extract_model_tests
)


# Confirm that the analysis contains the expected 16 tests

if (nrow(nat_tests) != 16) {
  stop(
    "Expected 16 burden-clinical variable tests, but obtained ",
    nrow(nat_tests),
    "."
  )
}


# Adjust P-values across all 16 tests

nat_tests <- nat_tests %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    significant = q_value < 0.05
  ) %>%
  mutate(
    clinical_characteristic = case_when(
      variable == "age_10y" ~ "Age (per 10-year increase)",
      variable == "sex" ~ "Sex",
      variable == "smoking_status" ~ "Smoking status",
      variable == "family_history_cancer" ~
        "Family history of cancer",
      TRUE ~ variable
    )
  ) %>%
  select(
    burden,
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


# Extract fixed-effect estimates and 95% confidence intervals
#
# These estimates provide the magnitude and direction of associations.
# Statistical significance is based on the omnibus F-tests above.

extract_coefficients <- function(model, burden) {
  
  estimates <- fixef(model)
  
  confidence_intervals <- confint(
    model,
    parm = "beta_",
    method = "Wald"
  )
  
  coefficient_table <- tibble(
    term = names(estimates),
    estimate = unname(estimates),
    CI_low = confidence_intervals[, 1],
    CI_high = confidence_intervals[, 2]
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        term == "age_10y" ~ "age_10y",
        str_starts(term, "sex") ~ "sex",
        str_starts(term, "smoking_status") ~ "smoking_status",
        str_starts(term, "family_history_cancer") ~
          "family_history_cancer",
        TRUE ~ NA_character_
      ),
      burden = burden
    ) %>%
    select(
      burden,
      variable,
      term,
      estimate,
      CI_low,
      CI_high
    )
  
  coefficient_table
}


nat_coefficients <- map2_dfr(
  models,
  names(models),
  extract_coefficients
)


# Add the omnibus adjusted p-value for each clinical variable
# to the corresponding coefficient estimates

nat_coefficients <- nat_coefficients %>%
  left_join(
    nat_tests %>%
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


# Check for singular model fits

singular_models <- map_lgl(
  models,
  ~ isSingular(.x, tol = 1e-4)
)

if (any(singular_models)) {
  
  warning(
    "Singular fit detected for: ",
    paste(
      names(singular_models)[singular_models],
      collapse = ", "
    )
  )
}


# Display primary results

nat_tests %>%
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
    )
  ) %>%
  print(n = Inf)


# Save results

write_csv(nat_tests, file.path(OUT_DIR,"nat_mixed_model_tests.csv"))

write_csv(nat_coefficients, file.path(OUT_DIR, "nat_mixed_model_coefficients.csv"))


# Secondary analysis: associations between NAT burden and tumor characteristics

# Each tumor characteristic is added separately to the primary mixed-effects model while adjusting for age, sex, smoking status, and family history.

# The tumor characteristics tested are:
#   - muscle invasiveness
#   - tumor grade
#   - neoadjuvant treatment
#   - prior cancer history

# P-values for the tumor characteristic are obtained using F-tests with Satterthwaite's approximation.

# Benjamini-Hochberg correction is performed separately across these secondary tests.


# Add tumor characteristics to the NAT analysis dataset

nat_tumor <- data %>%
  filter(tissue_class == "non-cancerous") %>%
  transmute(
    sample_id = sample_id,
    Donor = Donor,
    burden_total_pct = burden_total_pct,
    burden_gain_pct = burden_gain_pct,
    burden_loss_pct = burden_loss_pct,
    burden_cnloh_pct = burden_cnloh_pct,
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
    ),
    muscle_invasive = muscle_invasive,
    grade = grade,
    neoadjuvant = neoadjuvant,
    prior_cancer_status = prior_cancer_status
  )


# Define tumor characteristics

tumor_variables <- tribble(
  ~clinical_characteristic,     ~variable,
  "Muscle invasiveness",        "muscle_invasive",
  "Tumor grade",                "grade",
  "Neoadjuvant treatment",      "neoadjuvant",
  "Prior cancer",               "prior_cancer_status"
)


# Check that tumor characteristics are constant within each patient

tumor_variable_consistency <- nat_tumor %>%
  group_by(Donor) %>%
  summarise(
    across(
      all_of(tumor_variables$variable),
      ~ n_distinct(.x[!is.na(.x)])
    ),
    .groups = "drop"
  )

if (
  any(
    tumor_variable_consistency %>%
    select(all_of(tumor_variables$variable)) %>%
    as.matrix() > 1
  )
) {
  stop(
    "At least one tumor characteristic differs between NAT samples ",
    "from the same patient."
  )
}


# Fit one adjusted model for each burden and tumor characteristic

fit_tumor_characteristic_model <- function(
    burden,
    tumor_variable
) {
  
  model_data <- nat_tumor %>%
    filter(
      !is.na(age_10y),
      !is.na(sex),
      !is.na(smoking_status),
      !is.na(family_history_cancer),
      !is.na(.data[[tumor_variable]])
    ) %>%
    droplevels()
  
  if (
    n_distinct(model_data[[tumor_variable]]) < 2
  ) {
    stop(
      "The variable ",
      tumor_variable,
      " has fewer than two observed groups."
    )
  }
  
  formula <- as.formula(
    paste0(
      "log1p(",
      burden,
      ") ~ age_10y + sex + smoking_status + ",
      "family_history_cancer + ",
      tumor_variable,
      " + (1 | Donor)"
    )
  )
  
  lmerTest::lmer(
    formula,
    data = model_data
  )
}


# Run all burden-tumor characteristic combinations

tumor_test_grid <- crossing(
  burden = burden_variables,
  tumor_variables
)


tumor_models <- pmap(
  tumor_test_grid,
  function(
    burden,
    clinical_characteristic,
    variable
  ) {
    
    fit_tumor_characteristic_model(
      burden = burden,
      tumor_variable = variable
    )
  }
)


names(tumor_models) <- paste(
  tumor_test_grid$burden,
  tumor_test_grid$variable,
  sep = "__"
)


# Extract the omnibus F-test for each tumor characteristic

extract_tumor_test <- function(
    model,
    burden,
    clinical_characteristic,
    variable
) {
  
  anova_table <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  tumor_row <- anova_table %>%
    filter(term == variable)
  
  tibble(
    burden = burden,
    clinical_characteristic = clinical_characteristic,
    variable = variable,
    n_NATs = nobs(model),
    n_patients = n_distinct(model@frame$Donor),
    numerator_df = tumor_row$NumDF,
    denominator_df = tumor_row$DenDF,
    F_value = tumor_row$`F value`,
    p_value = tumor_row$`Pr(>F)`
  )
}


nat_tumor_tests <- pmap_dfr(
  list(
    model = tumor_models,
    burden = tumor_test_grid$burden,
    clinical_characteristic =
      tumor_test_grid$clinical_characteristic,
    variable = tumor_test_grid$variable
  ),
  extract_tumor_test
)


# Adjust P-values across the secondary testing family
#
# Four burdens x four tumor characteristics = 16 secondary tests.

nat_tumor_tests <- nat_tumor_tests %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    significant = q_value < 0.05
  ) %>%
  arrange(q_value)


# Extract model coefficients for the tumor characteristics

extract_tumor_coefficients <- function(
    model,
    burden,
    clinical_characteristic,
    variable
) {
  
  estimates <- fixef(model)
  
  confidence_intervals <- confint(
    model,
    parm = "beta_",
    method = "Wald"
  )
  
  coefficient_table <- tibble(
    term = names(estimates),
    estimate = unname(estimates),
    CI_low = confidence_intervals[, 1],
    CI_high = confidence_intervals[, 2]
  )
  
  coefficient_table %>%
    filter(
      str_starts(
        term,
        variable
      )
    ) %>%
    mutate(
      burden = burden,
      clinical_characteristic =
        clinical_characteristic,
      variable = variable,
      .before = 1
    )
}


nat_tumor_coefficients <- pmap_dfr(
  list(
    model = tumor_models,
    burden = tumor_test_grid$burden,
    clinical_characteristic =
      tumor_test_grid$clinical_characteristic,
    variable = tumor_test_grid$variable
  ),
  extract_tumor_coefficients
)


# Add omnibus P-values and adjusted P-values to coefficient table

nat_tumor_coefficients <- nat_tumor_coefficients %>%
  left_join(
    nat_tumor_tests %>%
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

nat_tumor_tests %>%
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
    )
  ) %>%
  print(n = Inf)


# Save secondary analysis

write_csv(nat_tumor_tests, file.path(OUT_DIR, "nat_tumor_characteristic_tests.csv"))
write_csv(nat_tumor_coefficients, file.path(OUT_DIR, "nat_tumor_characteristic_coefficients.csv"))