from adjustText import adjust_text
import gwaslab as gl

# =============================================================================
# Beta-beta replication plots 
#
# Two complementary, symmetric plots. They differ ONLY in which study's lead
# SNPs anchor the plot; everything else (axes meaning, filtering, labeling) is
# identical:
#   Plot 1 -> "Do lpd's hits replicate in amr?"  (lpd lead SNPs on the X axis)
#   Plot 2 -> "Do amr's hits replicate in lpd?"  (amr lead SNPs on the X axis)
#
# Convention: the ANCHOR study (whose hits are tested) is path1 = X axis;
#             the OTHER study (where we look for replication) is path2 = Y axis.
# All plotted lead variants are labeled (anno=True).
# =============================================================================

# ---- Input files ------------------------------------------------------------
amr_full = "../amr_genotools_regenie_harmonized.gwaslab.tsv.gz"
lpd_full = "../lpd_saige_meta_fixed_harmonized.gwaslab.tsv.gz"
amr_lead = "../lead_amr_genotools_regenie_harmonized.gwaslab.tsv.gz"
lpd_lead = "../lead_lpd_saige_meta_fixed_harmonized.gwaslab.tsv.gz"

amr_name = "amr_genotools_regenie"
lpd_name = "lpd_saige_meta_fixed"

# ---- Load the two full sumstats once (used for effect estimates on both axes)
amr = gl.Sumstats(amr_full, fmt="gwaslab", sep="\t", build="38")
lpd = gl.Sumstats(lpd_full, fmt="gwaslab", sep="\t", build="38")

# ---- Lead SNP lists (the "hits" of each study) ------------------------------
amr_lead_snps = gl.Sumstats(amr_lead, fmt="gwaslab", sep="\t", build="38")["SNPID"].tolist()
lpd_lead_snps = gl.Sumstats(lpd_lead, fmt="gwaslab", sep="\t", build="38")["SNPID"].tolist()
print(f"amr lead SNPs: {len(amr_lead_snps)}")
print(f"lpd lead SNPs: {len(lpd_lead_snps)}")


def replication_plot(anchor_full, other_full, lead_snps, anchor_name, other_name,
                     sig_level=0.05):
    """Beta-beta plot: do `anchor_name`'s lead SNPs replicate in `other_name`?

    anchor study -> path1 (X axis); other study -> path2 (Y axis).
    All lead SNPs that survive the both-datasets / allele-match merge are labeled.
    `sig_level` sets the replication/coloring threshold (P<sig_level). The default
    0.05 is nominal replication; pass 1e-6 to test at the discovery threshold.
    """
    # sig-level tag keeps the 0.05 filenames unchanged while making 1e-6 outputs distinct
    tag = "" if sig_level == 0.05 else f"_sig{sig_level:g}"

    # Diagnostic: how many anchor leads are found and replicate (P<sig_level) in the other study
    sub = other_full.data.loc[other_full.data["SNPID"].isin(lead_snps), ["SNPID", "P"]]
    n_present = sub["SNPID"].nunique()
    n_replicating = sub.loc[sub["P"] < sig_level, "SNPID"].nunique()
    print(f"{anchor_name}: {len(lead_snps)} lead SNPs; "
          f"{n_present} found in {other_name}; "
          f"{n_replicating} replicate at P<{sig_level:g}")

    return gl.compare_effect(
        path1=anchor_full,
        path2=other_full,
        snplist=lead_snps,
        label=[anchor_name, other_name, "Both", "None"],
        sig_level=sig_level,
        is_45_helper_line=True,
        anno=True,
        anno_kwargs={"bbox": {"facecolor": "lightblue", "alpha": 0.5,
                              "boxstyle": "round,pad=0.5"}, "fontsize": 8},
        adjust_text_kwargs_r={"expand_text": (4, 2.6)},
        adjust_text_kwargs_l={"expand_text": (2.2, 3.6)},
        allele_match=True,
        include_all=False,
        drop=True,
        save=f"./{anchor_name}_hits_in_{other_name}{tag}_beta_beta.png",
        save_kwargs={"dpi": 600, "facecolor": "white"},
    )


# --- Nominal replication (P<0.05) -------------------------------------------
# Plot 1: lpd hits -> replicate in amr?  (lpd on X)
replication_plot(lpd, amr, lpd_lead_snps, lpd_name, amr_name)

# Plot 2: amr hits -> replicate in lpd?  (amr on X)
replication_plot(amr, lpd, amr_lead_snps, amr_name, lpd_name)

# --- Replication tested at the discovery threshold (P<1e-6) ------------------
# Plot 3: lpd hits -> replicate in amr at P<1e-6?  (lpd on X)
replication_plot(lpd, amr, lpd_lead_snps, lpd_name, amr_name, sig_level=1e-6)

# Plot 4: amr hits -> replicate in lpd at P<1e-6?  (amr on X)
replication_plot(amr, lpd, amr_lead_snps, amr_name, lpd_name, sig_level=1e-6)
