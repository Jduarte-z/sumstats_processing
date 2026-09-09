library(tidyverse)
library(ggforestplot)
library(patchwork)

# Configuration

study_manifest <- tribble(
  ~file, ~study, ~ancestry,
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/EUR_P1_firstPass.tsv.gwaslab.tsv", "EUR P1", "EUR",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/EUR_P2_firstPass.tsv.gwaslab.tsv", "EUR P2", "EUR",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/AFR_P1_firstPass.tsv.gwaslab.tsv", "AFR P1", "AFR",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/AFR_P2_firstPass.tsv.gwaslab.tsv", "AFR P2", "AFR",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/NAT_P1_firstPass.tsv.gwaslab.tsv", "NAT P1", "NAT",
  "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/meta-analysis_v3/firstPass/NAT_P2_firstPass.tsv.gwaslab.tsv", "NAT P2", "NAT"
)

meta_files <- c(
  random = "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/MR-MEGA_3way_rsID/random_and_fixed/madma_random.out",
  fixed  = "/home/duartej3/isilon/0.JF/Backups/Brain_response_winter_2025-2026/tractor_manual_meta_brain/MR-MEGA_3way_rsID/random_and_fixed/madma_fixed.out"
)

meta_types <- c("random", "fixed")

variants <- tribble(
  ~gene, ~variant_id, ~label,
  "ITPKB",   "1:226737538:G:A",  "rs117185933",
  "SNCA",    "4:89704662:A:G",   "rs12512806",
  "SNCA",    "4:89704960:G:A",   "rs356182",
  "SNCA",    "4:89827044:T:A",   "rs920624",
  "LRRK2",   "12:40340400:G:A",  "G2019S",
  "RPS6KA2", "6:166745184:C:T",  "rs6456122",
  "HIP1R",   "12:122842051:G:T", "rs10847864"
)

ancestry_colours <- c(
  "EUR"  = "#EE2C2C",
  "AFR"  = "#1874CD",
  "NAT"  = "#008B45",
  "META" = "#000000"
)

out_dir <- "."
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

stopifnot(all(meta_types %in% names(meta_files)))

meta_label <- c(
  random = "Random Effects",
  fixed  = "Fixed Effects"
)

all_ids <- unique(variants$variant_id)

# Read ancestry-specific results

study_hits <- study_manifest %>%
  transmute(
    study,
    ancestry,
    row = map(
      file,
      ~ read_tsv(.x, show_col_types = FALSE) %>%
        filter(SNPID %in% all_ids) %>%
        transmute(
          variant_id = SNPID,
          estimate = BETA,
          se = SE,
          pvalue = P,
          n = N,
          AF = EAF
        )
    )
  ) %>%
  unnest(row)

# Read meta-analysis results

meta_hits <- imap_dfr(meta_files[meta_types], function(path, type) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(rs_number %in% all_ids) %>%
    transmute(
      meta_type = type,
      variant_id = rs_number,
      estimate = log(OR),
      se = OR_se,
      pvalue = `p-value`,
      n = n_samples,
      AF = eaf
    )
})

build_studies <- function(vid) {
  study_manifest %>%
    select(study, ancestry) %>%
    left_join(
      filter(study_hits, variant_id == vid),
      by = c("study", "ancestry")
    ) %>%
    select(study, ancestry, estimate, se, pvalue, n, AF)
}

build_meta <- function(vid, type) {
  meta_hits %>%
    filter(variant_id == vid, meta_type == type) %>%
    transmute(
      study = "Meta-analysis",
      ancestry = "META",
      estimate,
      se,
      pvalue,
      n,
      AF
    )
}

make_forest <- function(gene, variant_id, label, meta_type) {

  df_meta <- build_meta(variant_id, meta_type)

  if (nrow(df_meta) == 0) {
    warning(
      sprintf(
        "No %s meta result for %s (%s %s) - meta row omitted.",
        meta_type,
        variant_id,
        gene,
        label
      )
    )
  }

  df_all <- bind_rows(
    build_studies(variant_id),
    df_meta
  ) %>%
    mutate(study = fct_inorder(study))

  title <- sprintf(
    "%s %s %s Meta-analysis using TRACTOR sumstats per Ancestry",
    gene,
    label,
    meta_label[[meta_type]]
  )

  p_forest <- ggforestplot::forestplot(
    df = df_all,
    name = study,
    estimate = estimate,
    se = se,
    pvalue = pvalue,
    logodds = TRUE,
    colour = ancestry,
    shape = ancestry
  ) +
    labs(colour = "Ancestry", shape = "Ancestry") +
    xlab("OR (95% CI)") +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey60"
    ) +
    scale_colour_manual(values = ancestry_colours) +
    theme_bw() +
    theme(legend.position = "bottom")

  df_table <- df_all %>%
    mutate(
      study = fct_rev(study),

      # N represents ancestry-specific haplotypes in TRACTOR output
      MAC = pmin(AF, 1 - AF) * n,

      AF_lab = ifelse(
        is.na(AF),
        "—",
        sprintf("%.3f", AF)
      ),

      N_lab = ifelse(
        is.na(n),
        "—",
        format(n, big.mark = ",", trim = TRUE)
      ),

      MAC_lab = ifelse(
        is.na(MAC),
        "—",
        format(round(MAC), big.mark = ",", trim = TRUE)
      ),

      pvalue_lab = ifelse(
        is.na(pvalue),
        "—",
        sprintf("%.2e", pvalue)
      )
    )

  p_table <- ggplot(df_table, aes(y = study)) +
    ggforestplot::geom_stripes() +
    geom_text(aes(x = 0, label = AF_lab), size = 3.2) +
    geom_text(aes(x = 1, label = N_lab), size = 3.2) +
    geom_text(aes(x = 2, label = MAC_lab), size = 3.2) +
    geom_text(aes(x = 3, label = pvalue_lab), size = 3.2) +
    scale_x_continuous(
      limits = c(-0.5, 3.5),
      breaks = c(0, 1, 2, 3),
      labels = c("AF", "N", "MAC", "p-value"),
      position = "top"
    ) +
    theme_void() +
    theme(
      axis.text.x.top = element_text(
        face = "bold",
        size = 9
      ),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )

  p_combined <- p_forest + p_table +
    plot_layout(widths = c(5, 2.8)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 13
        )
      )
    )

  safe <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", x)
  }

  fname <- sprintf(
    "forest_tractor_%s_%s_%s_%s.png",
    safe(gene),
    safe(label),
    meta_type,
    run_stamp
  )

  ggsave(
    file.path(out_dir, fname),
    plot = p_combined,
    width = 10,
    height = 5,
    dpi = 600
  )

  message("Wrote ", fname)
}

plan <- tidyr::expand_grid(
  variants,
  meta_type = meta_types
)

pwalk(plan, make_forest)
