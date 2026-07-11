library(tidyverse)
library(ggforestplot)
library(patchwork)

# =====================================================================
# CONFIG  --  edit this block only
# =====================================================================

# --- 1. Per-phase SAIGE study sumstats (shared across ALL variants) ---
#   <-- EDIT the two paths to your SAIGE Phase 1 / Phase 2 sumstats.
study_manifest <- tribble(
  ~file,                                  ~study,      ~phase,
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/runSAIGE_and_AK_phase2B1-8/runSAIGE_phase1/runSaige_vcfs_ogCovar/plot/saigeFirstPass_phase1_harmonized_vcfs_og_covar_withBE_SE.gwaslab.tsv.gz", "Phase 1",   "P1",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/runGWASphase2B1-10/runSaige_withVCFs/plot/saigeFirstPass_phase2_harmonized_vcfs_withBE_SE.gwaslab.tsv.gz", "Phase 2",   "P2"
)

# --- 2. Meta-analysis result files, one per type ---------------------
#   (each file holds ALL variants; we filter by rs_number)
meta_files <- c(
  random = "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/runGWASphase2B1-10/meta-analysis_phase1_phase2B1-10_v4/runGWAMA_P1_OG_covar_vcfs_P2_vcfs/gwama_random.out",
  fixed  = "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/runGWASphase2B1-10/meta-analysis_phase1_phase2B1-10_v4/runGWAMA_P1_OG_covar_vcfs_P2_vcfs/gwama_fixed.out"
)

# --- 3. Which meta types to plot (subset of names(meta_files)) --------
#   Each requested type produces its OWN plot/file per variant.
meta_types <- c("random", "fixed")     

# --- 4. Variants of interest -----------------------------------------
#   label = short display name (rsID or protein change) used in title + filename.
variants <- tribble(
  ~gene,     ~variant_id,        ~label,
  "ITPKB", "1:226737538:G:A",     "rs117185933",
  "SNCA",  "4:89704662:A:G", "rs12512806",
  "SNCA",  "4:89704960:G:A", "rs356182",
  "SNCA",  "4:89827044:T:A", "rs920624",
  "LRRK2",   "12:40340400:G:A",  "G2019S",
  "RPS6KA2", "6:166745184:C:T",  "rs6456122",
  "HIP1R", "12:122842051:G:T", "rs10847864"
)

# --- 5. Global phase colour palette (shared by all variants) ---------
phase_colours <- c(
  "P1"   = "#200BDC",   
  "P2"   = "#b706e8",  
  "META" = "#000000"   
)

# --- 6. Output -------------------------------------------------------
out_dir <- "."                                  
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# =====================================================================
# END CONFIG  --  nothing below normally needs editing
# =====================================================================

stopifnot(all(meta_types %in% names(meta_files)))

# Pretty meta-type labels for titles
meta_label <- c(random = "Random Effects", fixed 
= "Fixed Effects")

all_ids <- unique(variants$variant_id)

# --- Read each study file ONCE, keep only variants of interest -------
study_hits <- study_manifest %>%
  transmute(
    study, phase,
    row = map(file, ~ read_tsv(.x, show_col_types = FALSE) %>%
                filter(SNPID %in% all_ids) %>%
                transmute(variant_id = SNPID,
                          estimate = BETA, se = SE, pvalue = P, n = N, AF = EAF))
  ) %>%
  unnest(row)

# --- Read each requested meta file ONCE, keep only variants of interest
#   OR_se in GWAMA is the SE of log(OR), i.e. SE of BETA.
meta_hits <- imap_dfr(meta_files[meta_types], function(path, type) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(rs_number %in% all_ids) %>%
    transmute(meta_type = type,
              variant_id = rs_number,
              estimate = log(OR), se = OR_se,
              pvalue = `p-value`, n = n_samples, AF = eaf)
})

# --- Per-variant assembly helpers ------------------------------------
# Studies: left_join onto the full manifest so a study MISSING the variant
# is kept as an NA row (drawn blank in the forest, "—" in the table).
build_studies <- function(vid) {
  study_manifest %>%
    select(study, phase) %>%
    left_join(filter(study_hits, variant_id == vid),
              by = c("study", "phase")) %>%
    select(study, phase, estimate, se, pvalue, n, AF)
}

build_meta <- function(vid, type) {
  meta_hits %>%
    filter(variant_id == vid, meta_type == type) %>%
    transmute(study = "Meta-analysis", phase = "META",
              estimate, se, pvalue, n, AF)
}

# --- Build + save one forest plot ------------------------------------
make_forest <- function(gene, variant_id, label, meta_type) {

  df_meta <- build_meta(variant_id, meta_type)
  if (nrow(df_meta) == 0)
    warning(sprintf("No %s meta result for %s (%s %s) - meta row omitted.",
                    meta_type, variant_id, gene, label))

  # meta row goes last
  df_all <- bind_rows(build_studies(variant_id), df_meta) %>%
    mutate(study = fct_inorder(study))

  title <- sprintf("%s %s %s Meta-analysis using SAIGE sumstats per Phase",
                   gene, label, meta_label[[meta_type]])

  # --- Forest panel ---
  p_forest <- ggforestplot::forestplot(
    df       = df_all,
    name     = study,
    estimate = estimate,
    se       = se,
    pvalue   = pvalue,
    logodds  = TRUE,
    colour   = phase,
    shape    = phase
  ) +
    labs(colour = "Phase", shape = "Phase") +
    xlab("OR (95% CI)") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = phase_colours) +
    theme_bw() +
    theme(legend.position = "bottom")

  # --- Right-hand table panel (AF / N / MAC / p-value) ---
  # ggforestplot orders rows top-to-bottom, so reverse the factor to match.
  df_table <- df_all %>%
    mutate(
      study = fct_rev(study),
      # MAC = minor-allele count; AF is the effect-allele freq, so use min(AF, 1-AF).
      # SAIGE N counts individuals (diploid genotypes), so multiply by 2.
      MAC = pmin(AF, 1 - AF) * n * 2,
      AF_lab     = ifelse(is.na(AF),     "—", sprintf("%.3f", AF)),
      N_lab      = ifelse(is.na(n),      "—", format(n, big.mark = ",", trim = TRUE)),
      MAC_lab    = ifelse(is.na(MAC),    "—", format(round(MAC), big.mark = ",", trim = TRUE)),
      pvalue_lab = ifelse(is.na(pvalue), "—", sprintf("%.2e", pvalue))
    )

  p_table <- ggplot(df_table, aes(y = study)) +
    ggforestplot::geom_stripes() +
    geom_text(aes(x = 0, label = AF_lab),     size = 3.2) +
    geom_text(aes(x = 1, label = N_lab),      size = 3.2) +
    geom_text(aes(x = 2, label = MAC_lab),    size = 3.2) +
    geom_text(aes(x = 3, label = pvalue_lab), size = 3.2) +
    scale_x_continuous(
      limits = c(-0.5, 3.5),
      breaks = c(0, 1, 2, 3),
      labels = c("AF", "N", "MAC", "p-value"),
      position = "top"
    ) +
    theme_void() +
    theme(
      axis.text.x.top = element_text(face = "bold", size = 9),
      plot.margin     = margin(t = 5, r = 5, b = 5, l = 5)
    )

  p_combined <- p_forest + p_table +
    plot_layout(widths = c(5, 2.8)) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
    )

  # --- Filename: gene_variant_metatype_timestamp ---
  safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  fname <- sprintf("forest_saige_%s_%s_%s_%s.png",
                   safe(gene), safe(label), meta_type, run_stamp)

  ggsave(file.path(out_dir, fname), plot = p_combined,
         width = 10, height = 5, dpi = 600)

  message("Wrote ", fname)
}

# --- Run: every variant x every requested meta type ------------------
plan <- tidyr::expand_grid(variants, meta_type = meta_types)
pwalk(plan, make_forest)

