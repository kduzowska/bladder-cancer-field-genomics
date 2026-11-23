# script6_distance.R

# Test whether cell-fraction-weighted sCNA burden in normal-adjacent tissue (NAT) varies with distance from the tumor.

# Separate linear mixed-effects models are fitted for total, gain, loss, and CN-LOH burden:

# log(1 + burden) ~ distance from tumor + (1 | patient)

# Only patients with NAT samples collected at two or more different distances are included.

# P-values are adjusted across the four burden measures using the Benjamini-Hochberg method.

# Input:
# scna_burden_with_clinical_data.csv

# Output:
# nat_distance_models.csv


# Packages

library(tidyverse)
library(lme4)
library(lmerTest)


# Load data

INPUT_DIR <- "path/to/input/files"
OUT_DIR <- "path/to/your/directory"

data <- read_csv(file.path(INPUT_DIR, "scna_burden_with_clinical_data.csv"), show_col_types = FALSE)


# Select NAT samples

nat_distance <- data %>%
  filter(tissue_class == "non-cancerous") %>%
  transmute(
    sample_id = sample_id,
    Donor = Donor,
    distance_cm = Distance,
    burden_total_pct = burden_total_pct,
    burden_gain_pct = burden_gain_pct,
    burden_loss_pct = burden_loss_pct,
    burden_cnloh_pct = burden_cnloh_pct
  )


# Distance must be numeric

if (!is.numeric(nat_distance$distance_cm)) {
  stop(
    "Distance must be stored as a numeric variable measured in cm."
  )
}


# Remove NAT samples without distance information

nat_distance <- nat_distance %>%
  filter(
    !is.na(distance_cm)
  )


# Restrict the analysis to patients with NAT samples collected
# at two or more different distances

nat_distance <- nat_distance %>%
  group_by(Donor) %>%
  filter(
    n() >= 2,
    n_distinct(distance_cm) >= 2
  ) %>%
  ungroup()


# Check that eligible samples remain

if (nrow(nat_distance) == 0) {
  stop(
    "No patients have NAT samples collected at two or more ",
    "different distances."
  )
}


cat(
  "Patients included in distance analysis:",
  n_distinct(nat_distance$Donor),
  "\n"
)

cat(
  "NAT samples included in distance analysis:",
  nrow(nat_distance),
  "\n"
)

cat(
  "Distance range (cm):",
  min(nat_distance$distance_cm),
  "to",
  max(nat_distance$distance_cm),
  "\n"
)


# Define genome-wide burden measures

burden_variables <- c(
  "burden_total_pct",
  "burden_gain_pct",
  "burden_loss_pct",
  "burden_cnloh_pct"
)


# Check burden values

if (anyNA(nat_distance[burden_variables])) {
  stop(
    "Missing sCNA burden values were detected among NAT samples ",
    "included in the distance analysis."
  )
}


# Fit one mixed-effects model for each burden measure

fit_distance_model <- function(burden) {
  
  formula <- as.formula(
    paste0(
      "log1p(",
      burden,
      ") ~ distance_cm + (1 | Donor)"
    )
  )
  
  lmerTest::lmer(
    formula,
    data = nat_distance
  )
}


distance_models <- burden_variables %>%
  set_names() %>%
  map(
    fit_distance_model
  )


# Extract model results

extract_distance_results <- function(
    model,
    burden
) {
  
  anova_table <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("variable")
  
  distance_test <- anova_table %>%
    filter(variable == "distance_cm")
  
  confidence_interval <- confint(
    model,
    parm = "distance_cm",
    method = "Wald"
  )
  
  tibble(
    burden = burden,
    
    n_NATs = nobs(model),
    
    n_patients = n_distinct(
      model@frame$Donor
    ),
    
    estimate_per_cm = unname(
      fixef(model)["distance_cm"]
    ),
    
    CI_low = confidence_interval[1, 1],
    
    CI_high = confidence_interval[1, 2],
    
    numerator_df = distance_test$NumDF,
    
    denominator_df = distance_test$DenDF,
    
    F_value = distance_test$`F value`,
    
    p_value = distance_test$`Pr(>F)`
  )
}


distance_results <- map2_dfr(
  distance_models,
  names(distance_models),
  extract_distance_results
)


# Adjust P-values across the four burden measures

distance_results <- distance_results %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    
    significant = q_value < 0.05
  ) %>%
  arrange(q_value)


# Check for singular model fits

singular_models <- map_lgl(
  distance_models,
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


# Display results

distance_results %>%
  mutate(
    across(
      c(
        estimate_per_cm,
        CI_low,
        CI_high,
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

write_csv(distance_results, file.path(OUT_DIR, "nat_distance_models.csv"))