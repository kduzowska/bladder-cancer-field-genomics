# sCNA burden analysis

This repository contains the R code used to calculate and analyze cell-fraction-weighted somatic copy number alteration (sCNA) burden in bladder cancer tumors and normal-adjacent tissue (NAT) samples.

### Biological questions

1. Do NAT samples contain sCNA burden in bladder cancer patients?
2. Is the sCNA burden associated with clinical parameters (tumor grade, stage, muscle-invasiveness) and/or patient characteristics (sex, age, smoking status)?

Mosaic chromosomal alterations were called using the MoChA WDL pipeline and filtered as described in the Methods of the accompanying manuscript. The analysis in this repository begins with the resulting MoChA-derived call table (mocha_calls.csv). The MoChA calling and filtering pipeline itself is not reproduced here.

The scripts are numbered in the order in which they were run, with paths specified at the top of each script.

---

## script1_calculating_sCN_burden.R

For each retained autosomal gain, loss, or copy-neutral loss of heterozygosity (CN-LOH) event, the genomic length of the event was multiplied by its MoChA-estimated cell fraction. Cell-fraction-weighted sCNA burden was calculated as the sum of these weighted event lengths divided by the total length of autosomes 1-22. Burden was calculated genome-wide and separately for each autosome, both for all sCNAs combined and separately for gains, losses, and CN-LOH. Samples with no retained sCNA events were assigned a burden of zero.

Chromosome lengths correspond to the GRCh37/hg19 reference genome. The total autosomal denominator used for genome-wide burden calculations is 2881033286 bp. 

Input:
- mocha_calls.csv

Output:
- scna_burden_per_sample.csv

---

## script2_prepare_clinical_data.R

This script merges the per-sample cell-fraction-weighted sCNA burden data with patient clinical information and prepares the variables used in downstream analyses. Clinical variables including age, sex, smoking status, family history of cancer, pathological T stage, grade, preoperative treatment, and prior cancer history are cleaned and standardized.

Pathological T stage is derived from the pTNM classification. pTa and pT1 are classified as non-muscle-invasive bladder cancer (NMIBC), pT2-pT4 as muscle-invasive bladder cancer (MIBC), and pTis is retained as a separate Tis category. The resulting dataset is saved as `scna_burden_with_clinical_data.csv`, that is the input to subsequent scripts.

Input:
- scna_burden_per_sample.csv
- Donor_clinical_data.csv - a csv file with clinical information about the patients (age, sex, tumor stage, grade, neoadjuvant therapy, smoking history)

Output:
- scna_burden_with_clinical_data.csv

---

## script3_tumor_clinical_associations.R

This script tests associations between genome-wide cell-fraction-weighted sCNA burden in tumor samples and clinical characteristics. Total, gain, loss, and CN-LOH burdens are analyzed separately.

Spearman's rank correlation is used for continuous variables, and Wilcoxon rank-sum tests are used for binary clinical variables. Smoking is analyzed as ever versus never smoking for the tumor analysis. P-values are adjusted across all tested burden-clinical characteristic combinations using the Benjamini-Hochberg method, with q <0.05 considered statistically significant.

Input:
- scna_burden_with_clinical_data.csv

Output:
- tumor_clinical_associations.csv

---

## script4_nat_mixed_models.R

This script tests associations between genome-wide cell-fraction-weighted sCNA burden in normal-adjacent tissue (NAT) samples and patient clinical characteristics using linear mixed-effects models.

For each burden measure (total, gain, loss, and CN-LOH), log(1 + burden) is modeled as a function of age per 10-year increase, sex, smoking status, and family history of cancer, with a patient-specific random intercept to account for multiple NAT samples from the same patient. P-values are obtained using F-tests with Satterthwaite's approximation and adjusted across the 16 burden-clinical variable combinations using the Benjamini-Hochberg method.

The script also performs a secondary analysis testing whether NAT sCNA burden is associated with tumor-related characteristics, including muscle invasiveness, tumor grade, neoadjuvant treatment, and prior cancer history. Each tumor characteristic is added separately to the adjusted mixed-effects model, and these secondary tests are corrected for multiple testing separately from the primary analysis.

Results are saved as separate tables containing omnibus statistical tests and model coefficient estimates with 95% confidence intervals.

Input:
- scna_burden_with_clinical_data.csv

Output:
- nat_mixed_model_tests.csv
- nat_mixed_model_coefficients.csv

---

## 05_nat_patient_level_sensitivity.R

This script performs a patient-level sensitivity analysis of the associations between NAT cell-fraction-weighted sCNA burden and clinical characteristics. For each patient, total, gain, loss, and CN-LOH burdens are averaged across all available NAT samples and converted to ranks. Ranked burdens are modeled using linear regression with age, sex, smoking status, and family history of cancer as covariates.

P-values are obtained from partial F-tests and adjusted across the 16 burden-clinical variable combinations using the Benjamini-Hochberg method. This analysis assesses whether the findings from the primary mixed-effects analysis are robust to the highly skewed burden distribution and the presence of many zero-burden NAT samples.

Input:
- scna_burden_with_clinical_data.csv

Output:
- nat_patient_level_sensitivity_tests.csv
- nat_patient_level_sensitivity_coefficients.csv

---

## 06_nat_distance_models.R

This script tests whether cell-fraction-weighted sCNA burden in NAT samples varies with distance from the tumor. Separate linear mixed-effects models are fitted for total, gain, loss, and CN-LOH burden using log(1 + burden) as the outcome, distance from the tumor in centimeters as the predictor, and a patient-specific random intercept.

Only patients with NAT samples collected at two or more different distances are included. P-values are adjusted across the four burden measures using the Benjamini-Hochberg procedure. Model estimates represent the change in log(1 + burden) per 1 cm increase in distance from the tumor.

Input:
- scna_burden_with_clinical_data.csv

Output:
- nat_distance_models.csv

---

## Software versions

R version 4.5.2 (2025-10-31 ucrt)
tidyverse 2.0.0

---

