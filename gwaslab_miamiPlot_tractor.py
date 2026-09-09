# # plot miami gwama tractor + highlight novel hits v4.1.6
import os
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import gwaslab as gl
import matplotlib as mpl

HIGHLIGHT = True          
# set to False to skip novel-hit highlighting

# Per-ancestry color pairs [dark, light] for alternating chromosomes
COLORS_DEFAULT = ["#000000", "#ABABAB"]


COLORS_BY_ANC = {
    "AFR": ["#1874CD","#78B3EE"],
    "NAT": ["#008B45","#A1DABD"],  
    "EUR": ["#EE2C2C", "#F38585"],
}

# Ancestries to process
ANCS = ["AFR", "EUR", "NAT"]  # adjust to your ancestry list
# ANCS = ['NAT']
# Config
FIGSIZE = (9, 5)         
BASE_W  = 9              
SCALE   = FIGSIZE[0] / BASE_W

DPI = 600                 # DPI = sharpness/file size, not layout
PAD_IN = 0.25 * SCALE     # padding scales with figure size


# #mpl.rcParams["axes.titlesize"] = 20 * SCALE

IN_DIR  = "../"
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


def process_ancestry(anc):
    anc_start = datetime.now()
    print(f"[{anc_start.strftime('%H:%M:%S')}] Processing ancestry: {anc}...")

    colors_anc = COLORS_BY_ANC.get(anc, COLORS_DEFAULT)

    inRandom   = os.path.join(IN_DIR,  f"gwama_random_{anc}_firstPass.gwaslab.tsv.gz")  
    inFixed    = os.path.join(IN_DIR,  f"gwama_fixed_{anc}_firstPass.gwaslab.tsv.gz")  
    outPlot    = os.path.join(OUT_DIR, f"gwama_{anc}_miami")

    if HIGHLIGHT:
        novelhitsRandom_file = os.path.join(IN_DIR, f"gwama_random_{anc}_novelHits.tsv")  
        pinpointsRandom = load_pinpoints(novelhitsRandom_file)
        novelhitsFixed_file = os.path.join(IN_DIR, f"gwama_fixed_{anc}_novelHits.tsv")    
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
        suffixes                  =['_R', '_F'],  
        # R - Random, F - Fixed
        build                     ="38",
        mode                      ="m",          
        cut                       =11,
        #skip                      =5,
        # Significance Lines
        sig_line                  =True,
        sig_level                 =1e-6,
        additional_line           =[5e-8],
        additional_line_color     =["gray"],
        # Titles
        # Set to none to add side titles later
        titles                    = None,          
        # Fonts
        font_family               ="DejaVu Sans",
        fontsize                  =8 * SCALE,
        # Annotation
        anno1                     ="GENENAME",
        anno2                     ="GENENAME",
        # Options: "tight", "expand"
        anno_style                ="right",       
        anno_fontsize             =10 * SCALE,
        anno_sig_level            =1e-6,
        # Colors
        colors=colors_anc,
        # Highlight loci
        highlight1                =pinpointsRandom if pinpointsRandom else None, 
        highlight2                =pinpointsFixed if pinpointsFixed else None, 

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
        f"TRACTOR-{anc} Meta-Analysis",
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

    anc_end = datetime.now()
    elapsed = (anc_end - anc_start).seconds
    print(f"[{anc_end.strftime('%H:%M:%S')}] Finished {anc}. Elapsed: {elapsed}s\n")


if __name__ == "__main__":
    print(f"Script started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    with ProcessPoolExecutor(max_workers=6) as executor:
        futures = [executor.submit(process_ancestry, anc) for anc in ANCS]
        for f in futures:
            f.result()

    print(f"\nScript finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
