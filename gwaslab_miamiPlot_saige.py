# # plot miami gwama tractor + highlight novel hits v4.1.6
import os
from datetime import datetime
import pandas as pd
import gwaslab as gl
import matplotlib as mpl

HIGHLIGHT = True          # set to False to skip novel-hit highlighting
COLORS    = ["#000000", "#ABABAB"]

# -------------------------
# Tweakable config
# -------------------------
FIGSIZE = (9, 5)         
BASE_W  = 9              
SCALE   = FIGSIZE[0] / BASE_W

DPI = 600                 # DPI = sharpness/file size, not layout
PAD_IN = 0.25 * SCALE     # padding scales with figure size


# #mpl.rcParams["axes.titlesize"] = 20 * SCALE

IN_DIR  = "../../"
OUT_DIR = "./"
os.makedirs(OUT_DIR, exist_ok=True)


def load_pinpoints(path: str) -> list[str]:
    if not os.path.exists(path):
        return []
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return []
    if df.empty:
        return []
    col = "SNPID" if "SNPID" in df.columns else df.columns[0]
    snps = df[col].dropna().astype(str).str.strip()
    snps = [s for s in snps.tolist() if s and s.lower() != "nan"]
    return list(dict.fromkeys(snps))  # dedupe, keep order


script_start = datetime.now()
print(f"Script started: {script_start.strftime('%Y-%m-%d %H:%M:%S')}\n")

run_start = datetime.now()
print(f"[{run_start.strftime('%H:%M:%S')}] Processing miami plot...")

inRandom   = os.path.join(IN_DIR,  f"gwama_random_firstPass.gwaslab.tsv.gz")
inFixed    = os.path.join(IN_DIR,  f"gwama_fixed_firstPass.gwaslab.tsv.gz")
outPlot    = os.path.join(OUT_DIR, f"gwama_miami")

if HIGHLIGHT:
    novelhitsRandom_file = os.path.join(IN_DIR, "gwama_random_novelHits.tsv")
    pinpointsRandom = load_pinpoints(novelhitsRandom_file)
    novelhitsFixed_file = os.path.join(IN_DIR, "gwama_fixed_novelHits.tsv")
    pinpointsFixed = load_pinpoints(novelhitsFixed_file)
else:
    pinpointsRandom = []
    pinpointsFixed = []


gl1= gl.Sumstats(inRandom, build="38", fmt="gwaslab")
gl2= gl.Sumstats(inFixed, build="38", fmt="gwaslab")
fig,log = gl.plot_miami2(
    # Load dataframes
    gl1,
    gl2,
    id1                       ="SNPID",
    id2                       ="SNPID",
    suffixes                  =['_R', '_F'],  # R - Random, F - Fixed
    build                     ="38",
    mode                      ="m",           # Options "mqq", "qqm"
    cut                       =11,
    #skip                      =5,
    # Significance Lines
    sig_line                  =True,
    sig_level                 =5e-8,
    additional_line           =[1e-6],
    additional_line_color     =["gray"],
    # Titles
    titles                    = None,          # Set to none to add side titles later
    # Fonts
    font_family               ="DejaVu Sans",
    fontsize                  =8 * SCALE,
    # Annotation
    anno1                     ="GENENAME",
    anno2                     ="GENENAME",
    anno_style                ="right",       # Options: "tight", "expand"
    anno_fontsize             =10 * SCALE,
    #anno_sig_level            =ANNO_SIG_LEVEL.get(anc, ANNO_SIG_LEVEL_DEFAULT),
    # Colors
    colors=COLORS,
    # Highlight loci
    highlight1                =pinpointsRandom if pinpointsRandom else None, 
    highlight2                =pinpointsFixed if pinpointsFixed else None, 
    highlight_color1          ="#E70B0B",     # black for random-specific novel hits
    highlight_color2          ="#E70B0B",     # black for fixed-specific novel hits
    repel_force               =0.1,
    # Fix x and y axes
    xtight                    =False,
    same_ylim                 =True,
    # Figure size
    fig_kwargs={
        "figsize": FIGSIZE,
        "dpi": DPI,},
    # Save options
    save=None
)

# Get axes from fig
ax_top, ax_bottom = fig.axes

# Expand annotations
ymin_t, ymax_t = ax_top.get_ylim()
ax_top.set_ylim(ymin_t, ymax_t * 1.5) 

ymin_b, ymax_b = ax_bottom.get_ylim()
ax_bottom.set_ylim(ymin_b, ymax_b * 1.25)

# Tick label size
CHR_TICK = 7 * SCALE
ax_top.tick_params(axis="both", which="both", labelsize=CHR_TICK)
ax_bottom.tick_params(axis="both", which="both", labelsize=CHR_TICK)

# Add side labels
labels = ["Random effects", "Fixed effects"]
for ax, label in zip((ax_top, ax_bottom), labels):
    ax.set_ylabel("")
    bbox = ax.get_position()
    y_center = (bbox.y0 + bbox.y1) / 2
    fig.text(
        0.97, y_center, label,
        rotation=90, va="center", ha="center", fontsize=12,
    )

# Add padding to suptitle
fig.subplots_adjust(left=0.08, right=0.92, top=0.88, bottom=0.10)

# Add suptitle
fig.suptitle(
    f"SAIGE Phase1 & Phase2 meta-analysis",
    fontsize=14, 
    fontweight="bold",
    x=0.02, 
    y=0.97,
    ha="left"
)

fig.savefig(outPlot, 
            dpi              =DPI, 
            facecolor        ="white", 
            bbox_inches      ="tight", 
            pad_inches       =PAD_IN)

run_end = datetime.now()
elapsed = (run_end - run_start).seconds
print(f"[{run_end.strftime('%H:%M:%S')}] Finished miami plot. Elapsed: {elapsed}s\n")

print("Done process ancestry multi core function.")

print(f"\nScript finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
