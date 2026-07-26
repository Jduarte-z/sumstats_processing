#designed in version 4.1.6
import gwaslab as gl


gl.options.set_option("data_directory", "/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/")

threads=12



input_file="tractorNAT_phase2_harmonized.gwaslab.tsv.gz"
output_file="tractorNAT_phase2_harmonized_with_rsIDs"

cores=12

ss = gl.Sumstats(input_file, build="38",    
    snpid="SNPID",
    chrom="CHR",
    pos="POS",
    ea="EA",
    eaf="EAF",
    nea="NEA",
    n="N",
    beta="BETA",
    se="SE",
    p="P",
    OR="OR",
    OR_95L="OR_95L",
    OR_95U="OR_95U",
    sep="\t",
    na_values=["NA","."],
    verbose=True
)


#strict data clean

ss.fix_chr(remove=True)
ss.fix_pos(remove=True)
ss.fix_allele(remove=True)
ss.fix_id(fixchrpos=False,fixid=True,fixsep=False,forcefixid=True,overwrite=True)
ss.normalize_allele(threads=threads)
ss.sort_coordinate()



# "empty" (default): Only assign rsID for variants with missing/NA rsID values. This is the safest option and preserves existing rsID assignments.
# "invalid": Assign rsID for variants with invalid rsID format (not matching the pattern rs[0-9]+). Useful for fixing incorrectly formatted rsIDs.
# "all": Overwrite all rsIDs for eligible variants, regardless of existing values. Use with caution as this will replace all existing rsID assignments.

ss.assign_rsid2(
vcf_path="/home/duartej3/isilon/0.JF/Backups/LPD_DATA/june_2026/LPD_DATA/reference4ProjectedPCA/dbSNP/GCF_000001405.40_release157.gz",
threads=threads,
overwrite="empty",
)

# Replace invalid rsIDs with custom CHR:POS:NEA:EA format (instead of dropping)
bad_tokens = {"", "NA", "N/A", ".", "nan", "NaN", "None", "<NA>"}
s = ss.data["rsID"].astype("string").str.strip()
mask_token = s.isin(bad_tokens)
mask_na    = s.isna()
invalid    = mask_token | mask_na

# Build CHR:POS:NEA:EA IDs for invalid rows
custom_ids = (
    ss.data["CHR"].astype(str) + ":" +
    ss.data["POS"].astype(str) + ":" +
    ss.data["NEA"].astype(str) + ":" +
    ss.data["EA"].astype(str)
)

n_invalid = int(invalid.sum())
ss.data.loc[invalid, "rsID"] = custom_ids.loc[invalid]
print(f"Replaced {n_invalid} unmapped rsIDs "
      f"({int(mask_token.sum())} token NAs, {int(mask_na.sum())} true NAs) "
      f"with CHR:POS:NEA:EA IDs — total rows unchanged at {len(ss.data)}.")

# Sanity check: no bad tokens remain in rsID
assert not ss.data["rsID"].astype("string").str.strip().isin(bad_tokens).any(), \
    "Some rsIDs are still invalid after replacement"
assert ss.data["rsID"].notna().all(), "Some rsIDs are still NA after replacement"

#Tidy output 
ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=output_file, fmt="gwaslab")
lead_variants = ss.get_lead(
    windowsizekb=500,
    sig_level=5e-8,
    anno=True,
    build="38",
    source="ensembl")
lead_variants.to_csv("lead_vars.txt", index=False, sep='\t')
