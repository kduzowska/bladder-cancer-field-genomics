# step9_driver_genes_clonal.R

# Purpose:
#   Annotate PyClone-VI clustered variants, identify confident somatic
#   variants affecting established BLCA driver genes, determine their
#   distribution across shared, tumor-specific and NAT-specific clonal
#   lineages, summarize the TERT promoter hotspot, and generate Figure 4.

# Confident somatic variant:
#   - gnomAD popmax < 0.001 or absent from gnomAD
#   - >10 alternate reads
#   - VAF >0.10
#   - ALT and VAF criteria must be satisfied in the same tissue sample

# Driver-gene variants:
#   - protein-altering SNVs in genes listed in config/blca_driver_genes.tsv
#   - TERT promoter chr5:1295113 G>A is explicitly assigned to TERT because
#     promoter annotation may not map it to TERT in Gene.refGene

# Outputs:
#   - annotated clustered variants
#   - confident somatic variants
#   - BLCA driver-gene variants by clonal class
#   - gene- and pathway-level summaries
#   - TERT C228T clonal summary
#   - Figure for manuscript


# 1. Packages

suppressPackageStartupMessages({
  library(tidyverse)
})

# 2. Paths

PYCLONE_RESULTS <- "path/to/step6/output/pyclone_all_results.tsv"
CLUSTER_SUMMARY <- "path/to/step6/output/pyclone_cluster_summary.tsv"

STEP2_DIR <- "path/to/step2/output/folder"
ANNOVAR_FILE <- "path/to/pyclone_annotated.hg38_multianno.txt"
DRIVER_FILE <- "path/to/blca_driver_genes.tsv"
OUT_DIR <- "path/to/output/directory"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 3. Analysis settings

POP_AF_MAX <- 0.001
ALT_READ_THRESHOLD <- 10
VAF_THRESHOLD <- 0.10

TERT_MUTATION_ID <- "chr5:1295113:G:A"

CLONAL_CLASSES <- c(
  "shared",
  "tumor_specific",
  "NAT_specific"
)

FIGURE_CLASS_LABELS <- c(
  shared = "Shared lineage",
  tumor_specific = "Tumor-specific lineage",
  NAT_specific = "NAT-specific lineage"
)

FIGURE_CLASS_COLORS <- c(
  shared = "#151D70",
  tumor_specific = "#C03830",
  NAT_specific = "#BBBBBC"
)

# 4. Helper functions

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_remove(x, "^chr")
  paste0("chr", x)
}

make_coord_key <- function(chr, pos, ref, alt) {
  paste(
    stringr::str_remove(normalize_chr(chr), "^chr"),
    as.integer(pos),
    stringr::str_to_upper(ref),
    stringr::str_to_upper(alt),
    sep = "|"
  )
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

get_optional_col <- function(df, candidates) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) == 0) return(rep(NA_character_, nrow(df)))
  as.character(df[[found[1]]])
}

parse_logical <- function(x) {
  x <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  dplyr::case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

# 5. Validate inputs

required_paths <- c(
  PYCLONE_RESULTS,
  CLUSTER_SUMMARY,
  ANNOVAR_FILE,
  DRIVER_FILE
)

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

if (!dir.exists(STEP2_DIR)) {
  stop("Step 2 directory not found: ", STEP2_DIR, call. = FALSE)
}


# 6. Load BLCA driver-gene panel

driver_panel <- readr::read_tsv(
  DRIVER_FILE,
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

required_driver_columns <- c("gene", "pathway")
missing_driver_columns <- setdiff(required_driver_columns, names(driver_panel))

if (length(missing_driver_columns) > 0) {
  stop(
    "blca_driver_genes.tsv is missing column(s): ",
    paste(missing_driver_columns, collapse = ", "),
    call. = FALSE
  )
}

driver_panel <- driver_panel %>%
  dplyr::transmute(
    gene = stringr::str_trim(gene),
    pathway = stringr::str_trim(pathway)
  )

if (any(is.na(driver_panel$gene) | driver_panel$gene == "" |
        is.na(driver_panel$pathway) | driver_panel$pathway == "")) {
  stop("Missing gene or pathway values detected in blca_driver_genes.tsv.", call. = FALSE)
}

if (anyDuplicated(driver_panel$gene)) {
  print(driver_panel %>% dplyr::count(gene) %>% dplyr::filter(n > 1))
  stop("Duplicated genes detected in blca_driver_genes.tsv.", call. = FALSE)
}

if (!"TERT" %in% driver_panel$gene) {
  stop("TERT must be present in blca_driver_genes.tsv.", call. = FALSE)
}

# 7. Load mutation-to-cluster assignments

pyclone <- readr::read_tsv(
  PYCLONE_RESULTS,
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

required_pyclone <- c("patient", "mutation_id", "cluster_id")
missing_pyclone <- setdiff(required_pyclone, names(pyclone))

if (length(missing_pyclone) > 0) {
  stop(
    "pyclone_all_results.tsv is missing column(s): ",
    paste(missing_pyclone, collapse = ", "),
    call. = FALSE
  )
}

mutation_clusters <- pyclone %>%
  dplyr::transmute(
    patient = stringr::str_trim(patient),
    mutation_id = stringr::str_trim(mutation_id),
    cluster_id = suppressWarnings(as.integer(cluster_id))
  ) %>%
  dplyr::distinct()

mutation_cluster_check <- mutation_clusters %>%
  dplyr::count(patient, mutation_id, name = "n_clusters") %>%
  dplyr::filter(n_clusters != 1L)

if (nrow(mutation_cluster_check) > 0) {
  print(mutation_cluster_check)
  stop("At least one mutation is assigned to more than one PyClone cluster.", call. = FALSE)
}

# 8. Load cluster classifications

clusters <- readr::read_tsv(
  CLUSTER_SUMMARY,
  col_types = readr::cols(.default = "c"),
  show_col_types = FALSE
)

required_cluster_columns <- c(
  "patient",
  "cluster_id",
  "classification",
  "n_mutations",
  "high_confidence_cluster",
  "tumor_truncal",
  "max_tumor_cp",
  "max_NAT_cp"
)

missing_cluster_columns <- setdiff(required_cluster_columns, names(clusters))

if (length(missing_cluster_columns) > 0) {
  stop(
    "pyclone_cluster_summary.tsv is missing column(s): ",
    paste(missing_cluster_columns, collapse = ", "),
    call. = FALSE
  )
}

clusters <- clusters %>%
  dplyr::transmute(
    patient = stringr::str_trim(patient),
    cluster_id = suppressWarnings(as.integer(cluster_id)),
    classification = stringr::str_trim(classification),
    n_cluster_mutations = suppressWarnings(as.integer(n_mutations)),
    high_confidence_cluster = parse_logical(high_confidence_cluster),
    tumor_truncal = parse_logical(tumor_truncal),
    max_tumor_cp = suppressWarnings(as.numeric(max_tumor_cp)),
    max_NAT_cp = suppressWarnings(as.numeric(max_NAT_cp)),
    median_assignment_prob = if ("median_assignment_prob" %in% names(clusters)) {
      suppressWarnings(as.numeric(clusters$median_assignment_prob))
    } else {
      NA_real_
    }
  )

if (anyDuplicated(clusters[c("patient", "cluster_id")])) {
  stop("Duplicated patient x cluster rows detected in cluster summary.", call. = FALSE)
}

# 9. Join mutations to cluster classifications

muts <- mutation_clusters %>%
  dplyr::left_join(clusters, by = c("patient", "cluster_id"))

if (any(is.na(muts$classification))) {
  print(muts %>% dplyr::filter(is.na(classification)) %>% head(20))
  stop("Some PyClone mutations could not be matched to a cluster classification.", call. = FALSE)
}

# 10. Parse mutation IDs

muts <- muts %>%
  tidyr::separate(
    mutation_id,
    into = c("chr", "pos", "ref", "alt"),
    sep = ":",
    remove = FALSE,
    convert = FALSE
  ) %>%
  dplyr::mutate(
    chr = normalize_chr(chr),
    pos = suppressWarnings(as.integer(pos)),
    ref = stringr::str_to_upper(ref),
    alt = stringr::str_to_upper(alt),
    coord_key = make_coord_key(chr, pos, ref, alt)
  )

bad_mutation_ids <- muts %>%
  dplyr::filter(
    is.na(pos) | pos < 1 |
      !ref %in% c("A", "C", "G", "T") |
      !alt %in% c("A", "C", "G", "T") |
      nchar(ref) != 1 |
      nchar(alt) != 1 |
      ref == alt
  )

if (nrow(bad_mutation_ids) > 0) {
  print(head(bad_mutation_ids, 20))
  stop("Malformed SNV mutation IDs detected.", call. = FALSE)
}


# 11. Load Step 2 tissue-level read support

step2_files <- sort(
  list.files(
    STEP2_DIR,
    pattern = "_candidate_somatic_long\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )
)

if (length(step2_files) == 0) {
  stop("No Step 2 *_candidate_somatic_long.tsv files found.", call. = FALSE)
}

read_step2_support <- function(path) {
  dat <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = "c"),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  required <- c(
    "patient", "sample_id", "role",
    "chr", "pos", "ref", "alt",
    "alt_count"
  )
  
  missing <- setdiff(required, names(dat))
  
  if (length(missing) > 0) {
    stop(
      basename(path),
      " is missing Step 2 column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!"mutation_id" %in% names(dat)) {
    dat <- dat %>%
      dplyr::mutate(
        mutation_id = paste(
          normalize_chr(chr),
          as.integer(pos),
          stringr::str_to_upper(ref),
          stringr::str_to_upper(alt),
          sep = ":"
        )
      )
  }
  
  if (!"vaf" %in% names(dat)) {
    if (!all(c("ref_count", "alt_count") %in% names(dat))) {
      stop("Step 2 file lacks both vaf and ref_count needed to calculate VAF.", call. = FALSE)
    }
    
    dat <- dat %>%
      dplyr::mutate(
        ref_count = suppressWarnings(as.numeric(ref_count)),
        alt_count_num = suppressWarnings(as.numeric(alt_count)),
        vaf = dplyr::if_else(
          ref_count + alt_count_num > 0,
          alt_count_num / (ref_count + alt_count_num),
          NA_real_
        )
      )
  }
  
  dat %>%
    dplyr::transmute(
      patient = stringr::str_trim(patient),
      mutation_id = stringr::str_trim(mutation_id),
      sample_id = stringr::str_trim(sample_id),
      role = stringr::str_to_lower(stringr::str_trim(role)),
      alt_count = suppressWarnings(as.numeric(alt_count)),
      vaf = suppressWarnings(as.numeric(vaf))
    ) %>%
    dplyr::filter(role %in% c("tumor", "adjacent"))
}

step2_support_long <- purrr::map_dfr(step2_files, read_step2_support)

read_support <- step2_support_long %>%
  dplyr::mutate(
    passes_sample_support =
      !is.na(alt_count) & alt_count > ALT_READ_THRESHOLD &
      !is.na(vaf) & vaf > VAF_THRESHOLD
  ) %>%
  dplyr::group_by(patient, mutation_id) %>%
  dplyr::summarise(
    max_alt_reads = safe_max(alt_count),
    max_vaf = safe_max(vaf),
    n_tissue_samples = dplyr::n_distinct(sample_id),
    n_samples_passing_support = sum(passes_sample_support, na.rm = TRUE),
    passes_read_support = any(passes_sample_support, na.rm = TRUE),
    .groups = "drop"
  )

# 12. Load ANNOVAR annotation

ann <- readr::read_tsv(
  ANNOVAR_FILE,
  col_types = readr::cols(.default = "c"),
  na = c("", ".", "NA"),
  show_col_types = FALSE,
  progress = FALSE
)

required_annovar <- c(
  "Chr", "Start", "Ref", "Alt",
  "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"
)

missing_annovar <- setdiff(required_annovar, names(ann))

if (length(missing_annovar) > 0) {
  stop(
    "ANNOVAR file is missing column(s): ",
    paste(missing_annovar, collapse = ", "),
    call. = FALSE
  )
}

popmax_candidates <- c(
  "gnomad312_AF_popmax",
  "gnomad312_genome_AF_popmax",
  "gnomAD_genome_AF_popmax",
  "AF_popmax"
)

popmax_col <- popmax_candidates[popmax_candidates %in% names(ann)][1]

if (is.na(popmax_col)) {
  stop(
    "No gnomAD popmax column was found in the ANNOVAR file. Expected one of: ",
    paste(popmax_candidates, collapse = ", "),
    call. = FALSE
  )
}

cat("Using population-frequency column:", popmax_col, "\n")

ann_sel <- tibble::tibble(
  coord_key = make_coord_key(ann$Chr, ann$Start, ann$Ref, ann$Alt),
  annovar_match = TRUE,
  Func.refGene = ann$Func.refGene,
  Gene.refGene = ann$Gene.refGene,
  ExonicFunc.refGene = ann$ExonicFunc.refGene,
  AAChange.refGene = ann$AAChange.refGene,
  popmax = suppressWarnings(as.numeric(ann[[popmax_col]])),
  CLNSIG = get_optional_col(ann, c("CLNSIG", "CLNSIG_ClinVar")),
  REVEL = get_optional_col(ann, c("REVEL", "REVEL_score")),
  cosmic = get_optional_col(
    ann,
    c("cosmic100_coding", "cosmic99_coding", "cosmic95_coding")
  )
)

if (anyDuplicated(ann_sel$coord_key)) {
  duplicated_ann <- ann_sel %>%
    dplyr::count(coord_key) %>%
    dplyr::filter(n > 1)
  
  print(head(duplicated_ann, 20))
  stop("Duplicated variant keys detected in ANNOVAR annotation.", call. = FALSE)
}

# 13. Combine cluster, read-support and annotation information

annotated <- muts %>%
  dplyr::left_join(read_support, by = c("patient", "mutation_id")) %>%
  dplyr::left_join(ann_sel, by = "coord_key") %>%
  dplyr::mutate(
    is_tert_promoter = mutation_id == TERT_MUTATION_ID,
    
    gene = stringr::str_trim(
      stringr::str_replace(Gene.refGene, "[;,].*$", "")
    ),
    
    gene = dplyr::if_else(is_tert_promoter, "TERT", gene),
    
    is_protein_altering =
      stringr::str_detect(
        stringr::str_to_lower(dplyr::coalesce(ExonicFunc.refGene, "")),
        "nonsynonymous|stopgain|stoploss|startloss"
      ) |
      stringr::str_detect(
        stringr::str_to_lower(dplyr::coalesce(Func.refGene, "")),
        "splicing"
      )
  ) %>%
  dplyr::left_join(driver_panel, by = "gene") %>%
  dplyr::mutate(
    is_driver_gene = !is.na(pathway),
    
    is_driver_gene_variant =
      is_driver_gene &
      (is_protein_altering | is_tert_promoter),
    
    passes_population_filter =
      is.na(popmax) | popmax < POP_AF_MAX,
    
    confident_somatic =
      passes_population_filter &
      dplyr::coalesce(passes_read_support, FALSE)
  )

# 14. Annotation QC

annotation_join_rate <- mean(annotated$annovar_match %in% TRUE)

cat(sprintf("ANNOVAR coordinate join rate: %.1f%%\n", 100 * annotation_join_rate))

if (annotation_join_rate < 1) {
  cat("Examples of mutations not found in the ANNOVAR table:\n")
  print(
    annotated %>%
      dplyr::filter(is.na(annovar_match)) %>%
      dplyr::select(patient, mutation_id) %>%
      head(20)
  )
}

# 15. Write annotated clustered variants

annotated_out <- annotated %>%
  dplyr::select(
    patient,
    mutation_id,
    chr,
    pos,
    ref,
    alt,
    gene,
    pathway,
    Func.refGene,
    ExonicFunc.refGene,
    AAChange.refGene,
    CLNSIG,
    REVEL,
    cosmic,
    popmax,
    max_alt_reads,
    max_vaf,
    n_samples_passing_support,
    passes_population_filter,
    passes_read_support,
    confident_somatic,
    is_protein_altering,
    is_tert_promoter,
    is_driver_gene,
    is_driver_gene_variant,
    cluster_id,
    classification,
    n_cluster_mutations,
    high_confidence_cluster,
    tumor_truncal,
    max_tumor_cp,
    max_NAT_cp,
    median_assignment_prob
  ) %>%
  dplyr::arrange(patient, chr, pos)

readr::write_tsv(
  annotated_out,
  file.path(OUT_DIR, "all_clustered_variants_annotated.tsv")
)

# 16. Confident somatic clustered variants

confident_variants <- annotated_out %>%
  dplyr::filter(confident_somatic)

readr::write_tsv(
  confident_variants,
  file.path(OUT_DIR, "confident_somatic_clustered_variants.tsv")
)

# 17. Confident BLCA driver-gene variants in interpretable clonal classes

driver_variants <- annotated_out %>%
  dplyr::filter(
    confident_somatic,
    is_driver_gene_variant,
    classification %in% CLONAL_CLASSES
  ) %>%
  dplyr::distinct(patient, mutation_id, gene, .keep_all = TRUE) %>%
  dplyr::arrange(pathway, gene, patient, mutation_id)

if (nrow(driver_variants) == 0) {
  stop("No confident BLCA driver-gene variants remain after filtering.", call. = FALSE)
}

readr::write_tsv(
  driver_variants,
  file.path(OUT_DIR, "BLCA_driver_gene_variants_by_clonal_class.tsv")
)

# 18. Gene-level clonal distribution

gene_variant_counts <- driver_variants %>%
  dplyr::count(gene, pathway, classification, name = "n_variants") %>%
  tidyr::complete(
    gene,
    pathway,
    classification = CLONAL_CLASSES,
    fill = list(n_variants = 0L)
  ) %>%
  tidyr::pivot_wider(
    names_from = classification,
    values_from = n_variants,
    names_glue = "n_{classification}_variants"
  )

gene_patient_counts <- driver_variants %>%
  dplyr::distinct(patient, gene, pathway, classification) %>%
  dplyr::count(gene, pathway, classification, name = "n_patients") %>%
  tidyr::complete(
    gene,
    pathway,
    classification = CLONAL_CLASSES,
    fill = list(n_patients = 0L)
  ) %>%
  tidyr::pivot_wider(
    names_from = classification,
    values_from = n_patients,
    names_glue = "n_{classification}_patients"
  )

gene_totals <- driver_variants %>%
  dplyr::group_by(gene, pathway) %>%
  dplyr::summarise(
    n_total_variants = dplyr::n(),
    n_total_patients = dplyr::n_distinct(patient),
    .groups = "drop"
  )

gene_summary <- gene_totals %>%
  dplyr::left_join(gene_variant_counts, by = c("gene", "pathway")) %>%
  dplyr::left_join(gene_patient_counts, by = c("gene", "pathway")) %>%
  dplyr::arrange(dplyr::desc(n_total_patients), dplyr::desc(n_total_variants), gene)

readr::write_tsv(
  gene_summary,
  file.path(OUT_DIR, "BLCA_driver_gene_clonal_distribution.tsv")
)

# 19. Pathway-level clonal distribution

pathway_long <- driver_variants %>%
  dplyr::count(pathway, classification, name = "n_variants") %>%
  tidyr::complete(
    pathway,
    classification = CLONAL_CLASSES,
    fill = list(n_variants = 0L)
  ) %>%
  dplyr::group_by(pathway) %>%
  dplyr::mutate(
    n_total_variants = sum(n_variants),
    proportion = n_variants / n_total_variants
  ) %>%
  dplyr::ungroup()

pathway_summary <- pathway_long %>%
  dplyr::select(pathway, classification, n_variants, n_total_variants, proportion) %>%
  tidyr::pivot_wider(
    names_from = classification,
    values_from = c(n_variants, proportion),
    names_glue = "{.value}_{classification}"
  ) %>%
  dplyr::arrange(dplyr::desc(n_total_variants), pathway)

readr::write_tsv(
  pathway_summary,
  file.path(OUT_DIR, "BLCA_driver_pathway_clonal_distribution.tsv")
)

# 20. TERT C228T clonal summary

tert_summary <- annotated_out %>%
  dplyr::filter(is_tert_promoter) %>%
  dplyr::select(
    patient,
    mutation_id,
    gene,
    pathway,
    max_alt_reads,
    max_vaf,
    n_samples_passing_support,
    passes_population_filter,
    passes_read_support,
    confident_somatic,
    cluster_id,
    classification,
    n_cluster_mutations,
    high_confidence_cluster,
    tumor_truncal,
    max_tumor_cp,
    max_NAT_cp,
    median_assignment_prob
  ) %>%
  dplyr::arrange(classification, patient)

if (nrow(tert_summary) == 0) {
  warning("TERT chr5:1295113 G>A was not found among PyClone-VI mutations.")
}

readr::write_tsv(
  tert_summary,
  file.path(OUT_DIR, "TERT_C228T_clonal_summary.tsv")
)

# 21. Figure

pathway_order <- pathway_summary %>%
  dplyr::arrange(n_total_variants) %>%
  dplyr::pull(pathway)

figure_data <- pathway_long %>%
  dplyr::mutate(
    pathway = factor(pathway, levels = pathway_order),
    classification = factor(classification, levels = CLONAL_CLASSES)
  )

figure4 <- ggplot(
  figure_data,
  aes(x = pathway, y = proportion, fill = classification)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(
    values = FIGURE_CLASS_COLORS,
    breaks = CLONAL_CLASSES,
    labels = FIGURE_CLASS_LABELS,
    name = "Clonal class"
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Proportion of driver variants"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    legend.position = "right"
  )

ggsave(
  file.path(OUT_DIR, "Figure_driver_pathway_clonal_distribution.pdf"),
  figure4,
  width = 8,
  height = 5
)

ggsave(
  file.path(OUT_DIR, "Figure_driver_pathway_clonal_distribution.png"),
  figure4,
  width = 8,
  height = 5,
  dpi = 300
)

# 22. QC summary

driver_variants_outside_classes <- annotated_out %>%
  dplyr::filter(
    confident_somatic,
    is_driver_gene_variant,
    !classification %in% CLONAL_CLASSES
  )

driver_annotation_qc <- tibble::tribble(
  ~metric, ~value,
  "clustered_mutations", nrow(annotated_out),
  "ANNOVAR_matched_mutations", sum(annotated$annovar_match %in% TRUE),
  "ANNOVAR_join_rate", annotation_join_rate,
  "mutations_with_read_support", sum(!is.na(annotated_out$max_alt_reads)),
  "confident_somatic_mutations", nrow(confident_variants),
  "confident_BLCA_driver_gene_variants", nrow(driver_variants),
  "shared_driver_variants", sum(driver_variants$classification == "shared"),
  "tumor_specific_driver_variants", sum(driver_variants$classification == "tumor_specific"),
  "NAT_specific_driver_variants", sum(driver_variants$classification == "NAT_specific"),
  "driver_variants_outside_three_classes", nrow(driver_variants_outside_classes),
  "TERT_C228T_patients", dplyr::n_distinct(tert_summary$patient)
)

readr::write_tsv(
  driver_annotation_qc,
  file.path(OUT_DIR, "driver_annotation_qc.tsv")
)

# 23. Analysis settings

analysis_settings <- tibble::tribble(
  ~parameter, ~value,
  
  "driver_gene_source",
  "config/blca_driver_genes.tsv",
  
  "population_frequency_filter",
  "gnomAD popmax < 0.001 or absent from gnomAD",
  
  "alternate_read_filter",
  ">10 alternate reads",
  
  "VAF_filter",
  ">0.10",
  
  "read_support_rule",
  "ALT >10 and VAF >0.10 must occur in the same tumor or NAT sample",
  
  "driver_variant_definition",
  "protein-altering SNV in predefined BLCA driver gene or TERT C228T promoter hotspot",
  
  "TERT_promoter_hotspot",
  TERT_MUTATION_ID,
  
  "clonal_classes",
  paste(CLONAL_CLASSES, collapse = ";"),
  
  "below_detection_clusters",
  "excluded from driver clonal-distribution summaries",
  
  "cluster_confidence_filter",
  "none",
  
  "Figure4_counting_unit",
  "unique patient x mutation driver variants"
)

readr::write_tsv(
  analysis_settings,
  file.path(OUT_DIR, "driver_annotation_settings.tsv")
)

# 24. Software versions

writeLines(
  c(
    "Step 9: BLCA driver-gene variants across PyClone-VI clonal compartments",
    "",
    paste0("PyClone mutation results: ", PYCLONE_RESULTS),
    paste0("PyClone cluster summary: ", CLUSTER_SUMMARY),
    paste0("Step 2 directory: ", STEP2_DIR),
    paste0("ANNOVAR file: ", ANNOVAR_FILE),
    paste0("Driver-gene file: ", DRIVER_FILE),
    paste0("gnomAD popmax column: ", popmax_col),
    "",
    "R session:",
    capture.output(sessionInfo())
  ),
  file.path(OUT_DIR, "driver_annotation_software_versions.txt")
)

# 25. Final report

cat("Clustered mutations:", nrow(annotated_out), "\n")
cat("Confident somatic mutations:", nrow(confident_variants), "\n")
cat("Confident BLCA driver-gene variants:", nrow(driver_variants), "\n")
cat("Shared driver variants:", sum(driver_variants$classification == "shared"), "\n")
cat("Tumor-specific driver variants:", sum(driver_variants$classification == "tumor_specific"), "\n")
cat("NAT-specific driver variants:", sum(driver_variants$classification == "NAT_specific"), "\n")
cat("TERT C228T patients:", dplyr::n_distinct(tert_summary$patient), "\n")
cat("Figure 4:", file.path(OUT_DIR, "Figure4_driver_pathway_clonal_distribution.pdf"), "\n")
cat("Output directory:", OUT_DIR, "\n")