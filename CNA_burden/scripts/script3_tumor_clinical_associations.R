# script3_tumor_clinical_associations.R

# Test associations between cell-fraction-weighted sCNA burden in tumor
# samples and clinical characteristics.

# Continuous clinical variables are tested using Spearman's rank correlation.
# Binary clinical variables are tested using the Wilcoxon rank-sum test.

# P-values are adjusted across all burden-clinical variable combinations using the Benjamini-Hochberg method.

# Input:
# scna_burden_with_clinical_data.csv

# Output:
# tumor_clinical_associations.csv


# Packages

library(tidyverse)

# Load data

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

data <- read_csv(file.path(INPUT_DIR, "scna_burden_with_clinical_data.csv"), show_col_types = FALSE)


# Select tumor samples

tumor <- data %>%
  filter(tissue_class == "cancerous")


# Check that each patient contributes only one tumor sample

# The Wilcoxon and Spearman tests used here assume independent observations.

multiple_tumors <- tumor %>%
  count(Donor) %>%
  filter(n > 1)

if (nrow(multiple_tumors) > 0) {
  stop(
    "More than one tumor sample is present for at least one patient. ",
    "The planned Wilcoxon/Spearman analysis assumes one independent tumor observation per patient.")
}


# Prepare tumor-specific clinical variables

# Smoking is analyzed as ever versus never smoking because the Wilcoxon rank-sum test requires two groups.

# Former and current smokers are combined as ever smokers.
# The three-level never/former/current variable is retained for the NAT
# mixed-effects analysis in a later script.

tumor <- tumor %>%
  mutate(
    smoking_ever = case_when(
      smoking_status == "never" ~ "never",
      smoking_status %in% c("former", "current") ~ "ever",
      TRUE ~ NA_character_
    ),
    
    smoking_ever = factor(
      smoking_ever,
      levels = c("never", "ever")
    )
  )


# Clean additional categorical variables used only in the tumor analysis

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

clean_factor <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  
  x[x %in% missing_codes] <- NA_character_
  
  factor(x)
}


tumor <- tumor %>%
  mutate(
    grade = clean_factor(Grade),
    neoadjuvant = clean_factor(Neoadjuvant),
    prior_cancer_status = clean_factor(prior_cancer)
  )


# Define genome-wide burden measures

burden_variables <- c(
  "burden_total_pct",
  "burden_gain_pct",
  "burden_loss_pct",
  "burden_cnloh_pct"
)


# Define clinical variables and statistical tests

clinical_tests <- tribble(
  ~clinical_characteristic,      ~variable,                 ~test,
  "Age",                         "age_years",                "spearman",
  "Sex",                         "sex",                      "wilcoxon",
  "Smoking",                     "smoking_ever",             "wilcoxon",
  "Family history of cancer",    "family_history_cancer",   "wilcoxon",
  "Muscle invasiveness",         "muscle_invasive",         "wilcoxon",
  "Grade",                       "grade",                    "wilcoxon",
  "Neoadjuvant treatment",       "neoadjuvant",             "wilcoxon",
  "Prior cancer",                "prior_cancer_status",      "wilcoxon"
)


# Check that all required variables are present

required_variables <- c(
  burden_variables,
  clinical_tests$variable
)

missing_variables <- setdiff(
  required_variables,
  names(tumor)
)

if (length(missing_variables) > 0) {
  stop(
    "The following required variables are missing from the input data: ",
    paste(missing_variables, collapse = ", ")
  )
}


# Function for Spearman correlation

run_spearman <- function(data, burden, variable) {
  
  d <- data %>%
    select(
      burden = all_of(burden),
      clinical = all_of(variable)
    ) %>%
    filter(
      !is.na(burden),
      !is.na(clinical)
    )
  
  if (nrow(d) < 3) {
    stop(
      "Fewer than three complete observations are available for ",
      variable,
      "."
    )
  }
  
  if (n_distinct(d$clinical) < 2) {
    stop(
      "The clinical variable ",
      variable,
      " has no variation among the available tumor samples."
    )
  }
  
  if (n_distinct(d$burden) < 2) {
    
    return(
      tibble(
        n = nrow(d),
        test = "Spearman rank correlation",
        group_1 = NA_character_,
        n_group_1 = NA_integer_,
        median_group_1 = NA_real_,
        group_2 = NA_character_,
        n_group_2 = NA_integer_,
        median_group_2 = NA_real_,
        effect = NA_real_,
        effect_type = "Spearman rho",
        p_value = NA_real_
      )
    )
  }
  
  test_result <- cor.test(
    d$burden,
    d$clinical,
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    n = nrow(d),
    test = "Spearman rank correlation",
    group_1 = NA_character_,
    n_group_1 = NA_integer_,
    median_group_1 = NA_real_,
    group_2 = NA_character_,
    n_group_2 = NA_integer_,
    median_group_2 = NA_real_,
    effect = unname(test_result$estimate),
    effect_type = "Spearman rho",
    p_value = test_result$p.value
  )
}


# Function for Wilcoxon rank-sum test

run_wilcoxon <- function(data, burden, variable) {
  
  d <- data %>%
    select(
      burden = all_of(burden),
      clinical = all_of(variable)
    ) %>%
    filter(
      !is.na(burden),
      !is.na(clinical)
    ) %>%
    mutate(
      clinical = droplevels(factor(clinical))
    )
  
  if (nlevels(d$clinical) != 2) {
    stop(
      "The variable '",
      variable,
      "' has ",
      nlevels(d$clinical),
      " non-missing groups. ",
      "The Methods specify a Wilcoxon rank-sum test, which requires exactly two groups.")
  }
  
  groups <- levels(d$clinical)
  
  group_summary <- d %>%
    group_by(clinical) %>%
    summarise(
      n = n(),
      median = median(burden),
      .groups = "drop"
    )
  
  test_result <- wilcox.test(
    burden ~ clinical,
    data = d,
    exact = FALSE
  )
  
  median_1 <- group_summary %>%
    filter(clinical == groups[1]) %>%
    pull(median)
  
  median_2 <- group_summary %>%
    filter(clinical == groups[2]) %>%
    pull(median)
  
  n_1 <- group_summary %>%
    filter(clinical == groups[1]) %>%
    pull(n)
  
  n_2 <- group_summary %>%
    filter(clinical == groups[2]) %>%
    pull(n)
  
  tibble(
    n = nrow(d),
    test = "Wilcoxon rank-sum test",
    group_1 = groups[1],
    n_group_1 = n_1,
    median_group_1 = median_1,
    group_2 = groups[2],
    n_group_2 = n_2,
    median_group_2 = median_2,
    effect = median_2 - median_1,
    effect_type = "Difference in medians: group 2 - group 1",
    p_value = test_result$p.value
  )
}


# Function that selects the appropriate test

run_clinical_test <- function(
    data,
    burden,
    variable,
    test
) {
  
  if (test == "spearman") {
    
    run_spearman(
      data = data,
      burden = burden,
      variable = variable
    )
    
  } else if (test == "wilcoxon") {
    
    run_wilcoxon(
      data = data,
      burden = burden,
      variable = variable
    )
    
  } else {
    
    stop(
      "Unknown statistical test: ",
      test
    )
  }
}


# Run all burden-clinical characteristic combinations

test_grid <- crossing(
  burden = burden_variables,
  clinical_tests
)


tumor_results <- pmap_dfr(
  test_grid,
  function(
    burden,
    clinical_characteristic,
    variable,
    test
  ) {
    
    result <- run_clinical_test(
      data = tumor,
      burden = burden,
      variable = variable,
      test = test
    )
    
    result %>%
      mutate(
        burden = burden,
        clinical_characteristic = clinical_characteristic,
        variable = variable,
        .before = 1
      )
  }
)


# Adjust P-values for multiple testing
#
# Benjamini-Hochberg correction is performed across all tumor burden-clinical characteristic combinations tested in this script

tumor_results <- tumor_results %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    significant = q_value < 0.05
  ) %>%
  arrange(q_value)


# Inspect results

tumor_results %>%
  select(
    burden,
    clinical_characteristic,
    test,
    n,
    effect,
    effect_type,
    p_value,
    q_value
  ) %>%
  print(n = Inf)


# Save results

write_csv(tumor_results, file.path(OUT_DIR, "tumor_clinical_associations.csv"))