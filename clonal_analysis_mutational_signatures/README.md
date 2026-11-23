# Clonal analysis of bladder cancer and matched adjacent urothelium

This project analyzes clonal expansion in bladder cancer patients from whole-exome sequencing (WES) data and identifies mutational patterns present in bladder tumors and adjacent histologically normal tissues (NATs).

## Context of the project

For 38 patients, a tumor sample, 1–2 adjacent tissue samples, and a matched blood sample were selected for WES analysis.

Raw sequencing data were aligned to the hg38 (UCSC) reference genome and processed using the nf-core/sarek pipeline. Germline variants were called using Platypus and somatic variants using Octopus. Octopus somatic calling was performed as paired tissue-versus-blood comparisons for each tumor or adjacent tissue sample.

The downstream workflow described here starts from the Octopus somatic VCFs and corresponding BAM files. ANNOVAR annotations are incorporated later for driver-gene analysis.

## Biological questions

1. Which mutational signatures are detectable in NATs?
2. Which somatic variants are present in both NATs and matched tumors, and do they belong to the same clonal lineage?
3. Which mutational signatures are present in NAT-specific, shared, and tumor-specific clonal lineages?
4. Which BLCA driver-gene variants occur in shared, tumor-specific, or NAT-specific clonal lineages?

## Overview

Variant calling was performed upstream using Octopus in tissue-versus-blood comparisons.

A variant detected in one tissue may be present at low variant allele frequency (VAF) in another tissue without being independently called there. Therefore, the first part of this workflow creates a patient-level union of candidate SNVs and measures read support for the same REF/ALT allele across all available BAM files from that patient. The resulting force-called data are used for mutational signature analysis and clonal inference with PyClone-VI.

```text
Octopus somatic VCFs + patient BAMs
        |
        v
step1_snv_calling.sh (step1_loop_patients.sh)
        |
        |-- patient-level union of candidate SNVs
        |-- REF/ALT allele counts across all patient BAMs
        v
force-called VCFs
        |
        v
step2_filter_forcecalled_snvs.R (step2_loop_patients.sh)
        |
        |-- matched-blood filtering
        |-- tissue evidence filtering
        |-- VAF calculation
        |-- QC
        v
candidate somatic SNV dataset
        |
        |-------------------------------> step3_mutational_signatures.R
        |                                  sample-level COSMIC v3.2 SBS analysis
        |
        v
step4_prepare_pyclone_input.R + ASCAT
        |
        v
PyClone-VI input files
        |
        v
step5_run_pyclone.sh
        |
        v
PyClone-VI results
        |
        v
step6_analyse_pyclone_results.R
        |
        |---> step7_clonal_clinical_associations.R
        |---> step8_cluster_mutational_signatures.R
        |---> step9_driver_variants_clonal_distribution.R + ANNOVAR
```

## Input data

For each patient, the workflow expects:

- one matched blood BAM
- one or more tumor and/or adjacent urothelium BAMs
- one or more Octopus somatic VCFs generated from tissue-versus-blood comparisons
- a hg38 reference FASTA matching the reference used for alignment and variant calling
- a sample metadata table defining the biological role of each sample

Later steps additionally require ASCAT output, clinical data, ANNOVAR annotations, and a predefined BLCA driver-gene table.

## Sample metadata

Sample identity and biological role are defined in `sample_metadata.tsv`.

Required columns are:

- patient
- vcf_sample_id
- sample_id
- role (blood, tumor, or adjacent)
- signature_pair_selected

All tumor and adjacent samples are eligible for sample-level signature fitting in Step 3. `signature_pair_selected` identifies the single predefined tumor and NAT sample per patient used for the paired tumor–NAT statistical comparison.

Example:

```text
patient    vcf_sample_id             sample_id       role       signature_pair_selected
PAT001     XG-4050-PAT001-blood      PAT001_blood    blood      FALSE
PAT001     XG-4050-PAT001-PT1        PAT001_PT1      tumor      TRUE
PAT001     XG-4050-PAT001-UM1        PAT001_UM1      adjacent   TRUE
```

---

# Step 1 - Patient-level targeted SNV force-calling

Script: `step1_snv_calling.sh` and `step1_loop_patients.sh`

## Purpose

Step 1 creates the union of SNVs identified by Octopus for a patient and measures read support for those same REF/ALT alleles in every available BAM from that patient.

This allows downstream analyses to distinguish a variant not independently called in a sample from a variant that was directly evaluated and supported by a small number of reads.

## Variant selection

The source Octopus VCFs are restricted to biallelic SNVs before normalization and reference validation.

The workflow retains SNVs regardless of functional annotation - the union can contain coding, intronic, UTR, promoter, and intergenic SNVs, provided that the variant is present in one of the source Octopus VCFs and is a biallelic SNV.

Indels, multinucleotide variants, symbolic alleles, and other non-SNV records are not included in this workflow.

Because the sequencing data are WES, genomic regions that can be evaluated are limited by the original sequencing/capture coverage.

## Reference validation

Each retained SNV is checked against the hg38 reference FASTA using `bcftools norm`.

Reference mismatches cause the workflow to stop.

## Patient-level union

For each patient, SNVs from all available tissue-versus-blood Octopus VCFs are filtered to biallelic SNVs, checked against hg38, sorted and deduplicated by exact REF/ALT allele producing `<PATIENT>_union_snvs.vcf.gz`

Each variant is represented by its genomic coordinate and specific REF/ALT allele.

## Target allele list

The union VCF is converted into an explicit allele target file containing:

```text
CHROM    POS       REF,ALT
chr1     123456    C,T
chr5     1295113   G,A
```

This file constrains downstream calling to the same allele originally identified by Octopus.

Output:
- <PATIENT>_target_alleles.tsv.gz

## Allele counting across BAM files

`bcftools mpileup` is run across all BAM files belonging to the patient.

Current read-level settings are:
- minimum mapping quality: 20
- minimum base quality:    20
- maximum depth per BAM:   5000
- full BAQ is enabled

The analysis requests FORMAT/AD and FORMAT/DP for downstream calculation of allele counts, sequencing depth, and VAF.

The resulting pileup is passed to `bcftools call` using the explicit target allele list so that the same REF/ALT allele is evaluated across all samples.

## How to run

Single patient:

```bash
bash step1_snv_calling.sh PATIENT_ID
```

All patients using the shell wrapper:

```bash
bash step1_loop_patients.sh
```

## Outputs

For each patient:
- <PATIENT>_union_snvs.vcf.gz
- <PATIENT>_union_snvs.vcf.gz.tbi
- <PATIENT>_target_alleles.tsv.gz
- <PATIENT>_target_alleles.tsv.gz.tbi
- <PATIENT>_mpileup.vcf.gz - primary input to Step 2.
- <PATIENT>_mpileup.vcf.gz.tbi
- <PATIENT>_mpileup.err
- STATUS.txt

---

# Step 2 - Filtering force-called SNVs

Scripts: `step2_filter_forcecalled_snvs.R`and `step2_loop_patients.sh`

## How to run

```bash
bash /path/to/step2_loop_patients.sh /path/to/step1/output /path/to/sample_metadata.tsv /path/to/output/folder
```

## Purpose

Step 2 converts the force-called VCF into a clean per-mutation/per-sample dataset.

The script extracts REF and ALT read counts, calculates allelic depth, total depth, and VAF, evaluates matched blood and evidence in tumor/adjacent tissue and retains candidate somatic SNVs.

## VCF extraction

Variant information is extracted using `bcftools query`.

The following fields are retained:
CHROM
POS
REF
ALT
SAMPLE
FORMAT/AD
FORMAT/DP

The script verifies that Step 1 produced one biallelic SNV allele per variant record. Unexpected non-SNV or multiallelic records cause the script to stop.

### Allelic depth

allelic_depth = REF_count + ALT_count

### Total depth

`FORMAT/DP` is stored separately as `total_depth`.

## Variant allele frequency

VAF is calculated from allele-specific counts:

VAF = ALT_count / (REF_count + ALT_count)

## Tissue presence definition

A variant is considered present in a tumor or adjacent urothelium sample when all three criteria are met:

1. VAF >= 0.05
2. ALT-supporting reads >= 4
3. total depth (FORMAT/DP) >= 20

The same definition of variant presence is used throughout the workflow.

## Matched-blood filter

Candidate somatic SNVs must have sufficient matched-blood coverage and little evidence for the alternate allele:

blood total depth (FORMAT/DP) >= 10
blood ALT-supporting reads <= 2
blood VAF < 0.02

## Candidate somatic SNV

A variant is retained and considered a candidate somatic SNV when the matched-blood filter passes and the variant is present in >=1 tumor or NAT.

## Tumor and adjacent-tissue evidence

For each retained mutation, the script records:

- the number of tumor samples
- the number of adjacent samples
- the number of tumor samples in which the mutation is present
- the number of adjacent samples in which the mutation is present
- if present in any tumor: yes/no
- if present in any adjacent sample: yes/no

## Outputs
 - <PATIENT>_all_variants_long.tsv - Complete long-format table that contains every SNV + sample combination
 - <PATIENT>_variant_summary.tsv - Variant-level summary that contains one row per SNV with matched-blood and tissue-filtering results
 - <PATIENT>_candidate_somatic_long.tsv - Candidate somatic long-format table that contains one row per retained mutation + sample and is the input for further analyses
 - <PATIENT>_candidate_somatic_wide.tsv - Candidate somatic wide-format table that contains one row per retained SNV with sample-specific counts and VAFs in separate columns and is primarily intended for manual inspection and QC.

### QC outputs

- <PATIENT>_sample_qc.tsv
- <PATIENT>_filter_qc.tsv
- <PATIENT>_blood_filter_qc.tsv
- <PATIENT>_filter_thresholds.tsv
- <PATIENT>_candidate_somatic_vaf_qc.pdf
- <PATIENT>_software_versions.txt

---

# Step 3 - Mutational signatures

Script: `step3_mutational_signatures.R`

This script performs sample-level SBS96 mutational signature analysis using the somatic SNVs generated by Step 2.

The script:

1. loads candidate somatic SNVs from Step 2
2. retains mutations classified as present in each tumor/adjacent sample
3. constructs per-sample SBS96 mutation matrices
4. retrieves COSMIC v3.2 GRCh38 SBS signatures
5. selects the predefined 37-signature panel
6. performs strict signature refitting with `max_delta = 0.004`
7. calculates absolute and relative signature contributions
8. assesses reconstruction accuracy using cosine similarity
9. assigns mutation-count confidence categories (>100 SNVs high, 50-100 medium, 1-49 low)
10. compares one predefined tumor and one predefined NAT sample per patient using paired Wilcoxon signed-rank tests
11. applies Benjamini-Hochberg correction across signatures
12. writes analysis tables, QC outputs, plots, and software versions

All tumor and NAT samples with at least one SNV are fitted. `signature_pair_selected` is used only to select the predefined tumor–NAT pair for the paired statistical comparison.

Intended use: Specify paths at the top of the script and run the script in RStudio.

**Required metadata columns:**

1. patient
2. vcf_sample_id
3. sample_id
4. role
5. signature_pair_selected

---

# Step 4 - Preparation of PyClone-VI input files

Script: `step4_prepare_pyclone_input.R`

This step merges Step 2 allele counts with ASCAT allele-specific copy-number and purity estimates to generate one PyClone-VI input file per patient.

Required PyClone-VI input columns are:

1. `mutation_id` - unique mutation identifier, consistent across samples
2. `sample_id` - unique sample identifier
3. `ref_counts` - reads supporting the reference allele
4. `alt_counts` - reads supporting the alternate allele
5. `normal_cn` - normal copy number at the locus
6. `major_cn` - ASCAT major copy number
7. `minor_cn` - ASCAT minor copy number
8. `tumour_content` - sample tumor-content parameter used by PyClone-VI

Tumor samples require usable ASCAT purity and copy-number segments. Adjacent samples use ASCAT allele-specific copy number when available, otherwise a diploid major_cn = 1, minor_cn = 1 fallback is used. Adjacent samples are modeled with tumour_content = 1.0. ASCAT goodness-of-fit and MAPD are retained for QC.

The script also handles systematic ASCAT sample-name normalization and optional exceptional mappings from `sample_id_recoding.tsv`.

---

# Step 5 - Running PyClone-VI

Script: `step5_run_pyclone.sh`

This shell script loops over the per-patient input files generated in Step 4, fits PyClone-VI models, and writes per-patient result files.

The manuscript analysis uses PyClone-VI v0.2.0 with a beta-binomial emission distribution.

Run with:

```bash
bash step5_run_pyclone.sh
```

Set the input and output directories inside the shell script before running it.

---

# Step 6 - Analysis of PyClone-VI clusters

Script: `step6_analyse_pyclone_results.R`

This script processes per-patient PyClone-VI result files and generates cluster- and patient-level summaries for downstream analyses.

It uses PyClone-VI cellular prevalence estimates together with sample roles from `sample_metadata.tsv` to classify clusters as shared, tumor-specific, NAT-specific, or below the detection threshold. Tumor-truncal status and cluster confidence are recorded separately.

**Cluster classification:**

shared: CP >= 0.10 in at least one tumor and at least one NAT sample

tumor_specific: CP >= 0.10 in at least one tumor and < 0.10 in all NAT samples

NAT_specific: CP >= 0.10 in at least one NAT and < 0.10 in all tumor samples

below_detection: CP < 0.10 in all included samples

Tumor truncal: maximum tumor CP > 0.85

High-confidence cluster: >5 mutations

Low-confidence clusters are retained and flagged rather than removed.

Patient-level shared, tumor-specific, and NAT-specific mutation fractions are calculated using mutations assigned to those three detectable classes.
The variants in below_detection clusters are reported separately and excluded from this denominator.

Inputs:

- sample_metadata.tsv
- PyClone-VI result files

Main outputs:

- `pyclone_all_results.tsv` - combined mutation-level PyClone-VI results
- `pyclone_cluster_sample_summary.tsv` - cluster x sample cellular-prevalence table
- `pyclone_cluster_summary.tsv` - main cluster-classification table
- `pyclone_patient_summary.tsv` - patient-level clonal summary used for downstream analyses and Supplementary Table 12
- `pyclone_sample_composition.tsv` - tumor/NAT sample composition per patient
- `per_patient/` - per-patient cluster summaries
- QC, settings, and software-version files

---

# Step 7 - Associations between clonal architecture and clinical characteristics

Script: `step7_clonal_clinical_associations.R`

This script integrates the patient-level PyClone-VI summary with clinical data to evaluate associations between clonal architecture and patient or tumor characteristics.

The script:

- tests shared and NAT-specific mutation fractions against age, sex, tumor grade, pathological T stage, muscle invasiveness, smoking status, pack-years, and preoperative treatment using non-parametric statistical tests
- applies Benjamini-Hochberg correction across all primary clinical association tests
- performs the predefined <5%, 5-25%, and >25% shared-fraction analyses
- compares the low (<5%) and high (>25%) shared groups
- performs the predefined >25% NAT-specific analysis
- assesses the number of NAT samples per patient as a technical sampling-depth QC using Spearman correlation

Main outputs:

- clonal_clinical_merged.tsv
- clonal_clinical_primary_tests.tsv
- shared_fraction_group_summary.tsv
- shared_low_vs_high_clinical_tests.tsv
- NAT_specific_high_group_summary.tsv
- NAT_specific_high_clinical_tests.tsv
- NAT_sample_count_qc.tsv
- clinical_variable_qc.tsv
- clonal_clinical_analysis_settings.tsv
- clonal_clinical_software_versions.txt

---

# Step 8 - Mutational signatures in clonal compartments

Script: `step8_cluster_mutational_signatures.R`

This script characterizes mutational processes in shared, tumor-specific, and NAT-specific PyClone-VI clonal compartments.

Mutations are grouped by patient and clonal class to construct SBS96 spectra. Spectra are fitted against a focused COSMIC v3.2 panel comprising age-related (SBS1, SBS5, SBS40), APOBEC-associated (SBS2, SBS13), and tobacco-associated (SBS92) signatures using `fit_to_signatures_strict` with `max_delta = 0.004`.

The script generates patient-level and cohort-pooled signature contributions, reconstruction QC, predefined age-related (SBS1 + SBS5 + SBS40) and APOBEC (SBS2 + SBS13) process scores, and paired within-patient comparisons between clonal classes using Wilcoxon signed-rank tests with Benjamini-Hochberg correction.

No minimum mutation-count or cosine-similarity threshold is used to exclude fitted patient x class spectra. Mutation count and cosine similarity are retained as QC measures.

Main outputs:

- cluster_signature_mutations.tsv
- cluster_signature_96_matrix_patient_class.tsv
- cluster_signature_96_matrix_pooled.tsv
- cluster_signature_patient_class_contributions.tsv
- cluster_signature_pooled_contributions.tsv
- cluster_signature_process_contributions.tsv
- cluster_signature_individual_tests.tsv
- cluster_signature_process_tests.tsv
- cluster_signature_fit_qc.tsv
- cluster_signature_analysis_settings.tsv
- cluster_signature_software_versions.txt

---

# Step 9 - BLCA driver variants across clonal compartments

Script: `step9_driver_genes_clonal.R`

This script integrates PyClone-VI cluster classifications, Step 2 read support, ANNOVAR annotations, and `blca_driver_genes.tsv` to characterize variants in predefined BLCA driver genes across shared, tumor-specific, and NAT-specific clonal lineages.

Confident somatic variants must have:

- gnomAD popmax < 0.001 or absent from gnomAD
- ALT-supporting reads > 10
- VAF > 0.10

The ALT-read and VAF criteria must be met in the same tumor or NAT sample.

Protein-altering SNVs in genes listed in `blca_driver_genes.tsv` are retained as BLCA driver-gene variants. The recurrent TERT promoter C228T variant (chr5:1295113 G>A, GRCh38) is explicitly assigned to TERT, because promoter annotation may not map it to TERT in Gene.refGene.

The script generates variant-, gene-, and pathway-level clonal-distribution tables, a TERT-specific clonal summary, and the pathway-level stacked bar plot used in the manuscript as figure 4.

Main outputs:

- all_clustered_variants_annotated.tsv
- confident_somatic_clustered_variants.tsv
- BLCA_driver_gene_variants_by_clonal_class.tsv
- BLCA_driver_gene_clonal_distribution.tsv
- BLCA_driver_pathway_clonal_distribution.tsv
- TERT_C228T_clonal_summary.tsv
- driver_annotation_qc.tsv
- driver_annotation_settings.tsv
- driver_annotation_software_versions.txt
- Figure_driver_pathway_clonal_distribution.pdf
- Figure_driver_pathway_clonal_distribution.png

---

# Reproducibility

The workflow uses:

```text
R 4.5.2
bcftools 1.9
PyClone-VI 0.2.0
ASCAT 3.2
COSMIC SBS reference signatures v3.2 (GRCh38)
ANNOVAR
GRCh38 / hg38 reference genome
```

Key R packages include:

```text
tidyverse 2.0.0
MutationalPatterns 3.20.1
BSgenome.Hsapiens.UCSC.hg38 1.4.5
GenomicRanges 1.62.1
IRanges 2.44.0
GenomeInfoDb 1.46.2
S4Vectors 0.48.0
scales 1.4.0
```

The R scripts write `sessionInfo()` or equivalent software-version information to their output directories.

All stages of the workflow should use the same:

sample metadata
sample roles
variant-presence thresholds
genome build
mutation identifiers

To prevent inconsistent variant or lineage classifications between analysis stages.
