#plinkGwasLPDcombined

import gwaslab as gl


#gl.options.set_option("data_directory", "/home/duartej3/beegfs/JF/programs/gwaslab/references_gwaslab416/")

threads=8

input_file="Mama_hg38_harmonized_withOR_SE.gwaslab.tsv.gz"
output_file="manhattan"
plot="manhattan.png"
title="MAMA GP2"

ss = gl.Sumstats(
    input_file, , build="38", fmt="gwaslab",
    verbose=True,
)


ss.plot_mqq(
    mode='mqq',
    #cut=14,
    #skip=5,
    title=title, 
    sig_level=1e-6,
    anno_sig_level=1e-6,
    anno="GENENAME",
    build="38",
    sig_line=True,
    additional_line=[5e-8],
    additional_line_color=['black'],
    font_family="DejaVu Sans",
    fontsize=11,
    anno_fontsize=12,
    colors=["#000000","#ABABAB"],
    save=plot, save_kwargs={"dpi":300, "facecolor":"white"}
)
# lead_variants = ss.get_lead(
#     windowsizekb=500,
#     sig_level=5e-8,
#     anno=True,
#     build="38",
#     source="ensembl")
# lead_variants.to_csv("lead_vars.txt", index=False, sep='\t')

