#!/usr/bin/env Rscript

## ---- SETUP (NO INSTALLS) ----
suppressPackageStartupMessages({
  library(MungeSumstats)
  library(data.table)
  library(tidyGWAS)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
})

## Root where those cohort dirs/files live
base_dir  <- "/fh/fast/sinnott-armstrong_n/proj/TB_Meta_Mega_Folder/11_18_Munge_Files"

## Path to dbSNP155 directory used by tidyGWAS
dbsnp_dir <- "/fh/fast/sinnott-armstrong_n/proj/TB_Meta_Mega_Folder/Munge_Format_Testing/dbSNP155"

## ---- COHORT METADATA (NO AoU) ----
cohorts <- tibble::tribble(
  ~cohort,        ~in_path,                                                             ~ref_genome, ~CaseN, ~ControlN,
  # UK Biobank
  "UKBB",         file.path(base_dir, "UKBB",          "UK_Names_Fixed.tsv"),   "GRCh37",     3011,   394065,
  # China Kaduri
  "ChinaKaduri",  file.path(base_dir, "CKB",  "China_A16.tsv"),               "GRCh37",      488,    75529,
  # Estonian Biobank
  "EBB",          file.path(base_dir, "EBB",           "EBB_Raw.txt"),                 "GRCh38",     3181,   204942,
  # Finngen
  "Finngen",      file.path(base_dir, "Fingenn",       "finngen_wide.tsv"),            "GRCh38",     3063,   500348,
  # Genes & Health
  "GNH",          file.path(base_dir, "GNH",           "GNH_No_Id.tsv"),   "GRCh38",     1936,    44481,
  # Multi Ethnic Meta – Africa
  "GWAMA_Africa", file.path(base_dir, "GWAMA_Africa",  "GWAMA_Africa_NoEffects.txt"),      "GRCh37",     3104,     4316,
  # Multi Ethnic Meta – Asia
  "GWAMA_Asia",   file.path(base_dir, "GWAMA_Asia",    "GWAMA_Asia_NoEffects.txt"),        "GRCh37",     3929,     6763,
  # Multi Ethnic Meta – Europe
  "GWAMA_Europe", file.path(base_dir, "GWAMA_Europe",  "GWAMA_Europe_NoEffects.txt"),      "GRCh37",     6739,    13402,
  # Biobank Japan
  "BBJ",          file.path(base_dir, "BBJ",           "Japan_One_P.txt"),               "GRCh37",     7800,   170871,
  # Mass Gen Brigham
  "MGB",          file.path(base_dir, "MGB",           "MGB.tsv"),                     "GRCh38",      447,    41361,
  # Soumya
  "Soumya",       file.path(base_dir, "Soumya",        "Soumya_clean.tsv"),      "GRCh37",     2160,     1820,
  # Taiwan Biobank
  "Taiwan",       file.path(base_dir, "Taiwan",        "Taiwan_Fixed.txt"),              "GRCh38",     2555,   320283,
  # Zheng
  "Zheng",        file.path(base_dir, "Zheng",         "Zheng_Fixed.txt"),             "GRCh38",     2949,     5090
)

## ---- FUNCTION TO PROCESS ONE COHORT ----
process_cohort <- function(cohort, in_path, ref_genome, CaseN, ControlN, dbsnp_dir) {
  message("=== Processing ", cohort, " ===")

  if (!file.exists(in_path)) {
    warning("Input not found for ", cohort, ": ", in_path)
    return(NULL)
  }

  out_dir <- dirname(in_path)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_gz      <- file.path(out_dir, paste0(cohort, "_Munge.tsv.gz"))
  out_tidy    <- file.path(out_dir, paste0(cohort, "_Post_Munge_And_Tidy.tsv"))
  out_qq_pdf  <- file.path(out_dir, paste0(cohort, "_QQ_FRQbin.png"))

  ## 1. MungeSumstats
  Hmm <- format_sumstats(
    in_path,
    ref_genome            = ref_genome,
    dbSNP                 = 155,
    bi_allelic_filter     = FALSE,
    flip_frq_as_biallelic = TRUE,
    save_path             = out_gz
  )

  ## 2. Read back in
  Munge_Stats <- data.table::fread(Hmm)

  ## 3. tidyGWAS
  if (!is.na(CaseN) && !is.na(ControlN)) {
    Next_Step <- tidyGWAS(
      Munge_Stats,
      dbsnp_dir,
      CaseN    = CaseN,
      ControlN = ControlN
    )
  } else {
    message("  No CaseN/ControlN provided for ", cohort,
            " – running tidyGWAS without explicit Ns")
    Next_Step <- tidyGWAS(
      Munge_Stats,
      dbsnp_dir
    )
  }

  data.table::fwrite(
    Next_Step,
    file = out_tidy,
    sep  = "\t"
  )

  ## 4. QQ plot by FRQ decile
  if (!"FRQ" %in% names(Next_Step)) {
    warning("No FRQ column in tidyGWAS output for ", cohort, "; skipping QQ plot.")
    return(invisible(Next_Step))
  }

  Next_Step_Marked <- Next_Step %>%
    mutate(FRQ_bin = ggplot2::cut_number(FRQ, 10, dig.lab = 5)) %>%
    arrange(P) %>%
    group_by(FRQ_bin) %>%
    mutate(
      Expected = -log10(ppoints(n())),
      Observed = -log10(P)
    ) %>%
    ungroup()

  QQs <- ggplot(Next_Step_Marked, aes(x = Expected, y = Observed)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    facet_wrap(~ FRQ_bin, scales = "free") +
    labs(
      title = paste0(cohort, ": QQ Plots by FRQ Bin"),
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    )
ggsave(
    filename = out_qq_pdf,
    plot     = QQs,
    width    = 8,
    height   = 6
  )

  message("Finished ", cohort,
          "\n  Munge: ", out_gz,
          "\n  Tidy : ", out_tidy,
          "\n  QQ   : ", out_qq_pdf)

  invisible(Next_Step)
}

## ---- CLI LOGIC: ALL COHORTS OR ONE ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  ## No args -> run all cohorts sequentially
  message("No index or cohort name supplied; running ALL cohorts.")
  results <- purrr::pmap(
    cohorts,
    ~ process_cohort(
      cohort     = ..1,
      in_path    = ..2,
      ref_genome = ..3,
      CaseN      = ..4,
      ControlN   = ..5,
      dbsnp_dir  = dbsnp_dir
    )
  )
} else {
  key <- args[1]
  ## Allow numeric index OR cohort name
  idx_num <- suppressWarnings(as.integer(key))

  if (!is.na(idx_num)) {
    if (idx_num < 1 || idx_num > nrow(cohorts)) {
      stop("Index out of range: ", idx_num,
           " (must be between 1 and ", nrow(cohorts), ")")
    }
    row <- cohorts[idx_num, ]
  } else {
    if (!key %in% cohorts$cohort) {
      stop("Cohort name '", key, "' not found. Available: ",
           paste(cohorts$cohort, collapse = ", "))
    }
    row <- cohorts[cohorts$cohort == key, ]
  }

  process_cohort(
    cohort     = row$cohort,
    in_path    = row$in_path,
    ref_genome = row$ref_genome,
    CaseN      = row$CaseN,
    ControlN   = row$ControlN,
    dbsnp_dir  = dbsnp_dir
  )
}


