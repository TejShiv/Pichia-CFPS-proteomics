# figure_control_proteome.py
#
# Shivakumar et al. — Mechanistic insights into microsomes in a
# Pichia pastoris cell-free system
#
# Generates QC and main figures for the control condition proteomic
# profiling (Figure 2 and associated supplementary panels).
#
# Outputs
# -------
# QC_figure.pdf / QC_figure.png
#   Six-panel quality-control figure:
#     A  Log2 intensity distributions across all six fractions
#     B  Protein detection counts per fraction
#     C  Protein amount loaded per fraction (BCA assay values)
#     D  Effect of normalisation on microsome fraction medians
#     E  Inter-fraction Pearson correlation heatmap (log2 intensities)
#     F  Protein rank-abundance plot
#
# Main_figure.pdf / Main_figure.png
#   Five-panel main results figure:
#     A  3-way Venn diagram — protein detection overlap (control fractions)
#     B  Stacked bar — GO Slim functional composition by fraction
#     C  Top 15 proteins — control lysate
#     D  Top 15 proteins — control microsomes
#     E  Secretory pathway detection table
#
# Input:  Proteomics_raw_data.xlsx
# Usage:  python figure_control_proteome.py
#
# Dependencies: pandas, numpy, matplotlib, seaborn, scipy, openpyxl
#               matplotlib-venn (optional; built-in fallback renderer provided)

# ── Imports ────────────────────────────────────────────────────────────────
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from scipy import stats

warnings.filterwarnings("ignore")
matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype":  42,
    "font.family":  "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        9,
    "axes.titlesize":  10,
    "axes.labelsize":   9,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
    "legend.fontsize":  7,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "axes.edgecolor":   "#333333",
    "axes.linewidth":   0.8,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})

try:
    from matplotlib_venn import venn3, venn3_circles
    HAS_VENN = True
except ImportError:
    HAS_VENN = False
    print("[INFO] matplotlib-venn not installed — using built-in Venn renderer.")

# =============================================================================
# 0.  CONFIG
# =============================================================================

DATA_PATH    = "Proteomics_raw_data.xlsx"
OUT_QC_PNG   = "QC_figure.png"
OUT_QC_PDF   = "QC_figure.pdf"
OUT_MAIN_PNG = "Main_figure.png"
OUT_MAIN_PDF = "Main_figure.pdf"

# ── Protein loaded per fraction (µg) ──────────────────────────────────────
# Measured by BCA/Bradford assay prior to MS injection.
# BfA microsomes: 250 µg vs control microsomes: 98 µg → 2.55× loading difference.
# Without this correction every BfA microsomal protein appears falsely enriched.
PROTEIN_LOADED_UG = {
    "BfA_lysate":       46,
    "BfA_supernatant":  23,
    "BfA_microsomes":  250,
    "lysate":           44,
    "supernatant":      31,
    "microsomes":       98,
}
REFERENCE_UG = 50   # scale to this notional amount; keeps intensities readable

QUANT_COLS   = ["BfA_lysate", "BfA_supernatant", "BfA_microsomes",
                "lysate",     "supernatant",      "microsomes"]
BFA_COLS     = ["BfA_lysate", "BfA_supernatant", "BfA_microsomes"]
CTRL_COLS    = ["lysate",     "supernatant",      "microsomes"]

FRACTION_LABELS = {
    "BfA_lysate":       "BfA\nLysate",
    "BfA_supernatant":  "BfA\nSupernatant",
    "BfA_microsomes":   "BfA\nMicrosomes",
    "lysate":           "Ctrl\nLysate",
    "supernatant":      "Ctrl\nSupernatant",
    "microsomes":       "Ctrl\nMicrosomes",
}

# Colour palette
PAL = {
    "Translation":   "#1B4F72",
    "Energy":        "#F0B27A",
    "Chaperone":     "#27AE60",
    "Redox":         "#6C3483",
    "Proteolysis":   "#1ABC9C",
    "AminoAcid":     "#922B21",
    "ERGolgi":       "#E74C3C",
    "Glycosylation": "#F48FB1",
    "NucleicAcid":   "#00BCD4",
    "Other":         "#95A5A6",
    "bfa":           "#E07B39",
    "ctrl":          "#4C78A8",
}

# Consistent panel-letter offset (axes-fraction coordinates)
PANEL_LETTER_KW = dict(fontsize=15, fontweight="bold", va="top", ha="left",
                        transform=None)  # set per-axis below

# =============================================================================
# 1.  GO SLIM ANNOTATION MAP
#     Curated map: UniProt accession / gene symbol → GO Slim category
#     Extend this as needed; unmapped proteins fall to "Other/Uncharacterised"
# =============================================================================

GO_SLIM_MAP = {
    # ── Translation / Ribosomal ──────────────────────────────────────────
    **{g: "Translation/Ribosomal" for g in [
        "RPS17B","RPS12","RPS13","RPS7","RPS7B","RPS19B","RPS1","RPS4A","RPS4B",
        "RPS8A","RPS9A","RPS9B","RPS10","RPS11","RPS14","RPS15","RPS16",
        "RPS18","RPS20","RPS21","RPS22","RPS23","RPS24","RPS25","RPS26",
        "RPS27","RPS28","RPS29","RPS30","RPS31",
        "RPL5","RPL7A","RPL10","RPL12","RPL14","RPL15","RPL16","RPL17",
        "RPL18","RPL19","RPL20A","RPL20B","RPL21","RPL22","RPL23","RPL24",
        "RPL25","RPL26","RPL27","RPL28","RPL29","RPL30","RPL31","RPL32",
        "RPL33","RPL34","RPL35","RPL36","RPL37","RPL38","RPL39","RPL40",
        "RPL41","RPL42","RPL43","RPP2B","RPP0","RPP1","RPP2",
        "EFT1","EFT2","EFB1","TEF1","TEF2","TEF4","YEF3","SUI1","SUI1-1",
        "TIF1","TIF2","STM1","ASC1","EGD1","EGD2",
    ]},
    # ── Energy metabolism ────────────────────────────────────────────────
    **{g: "Energy metabolism" for g in [
        "TDH1","TDH2","TDH3","ENO1","ENO2","PGK1","GPM1","FBA1","TPI1",
        "PGI1","GLK1","HXK1","HXK2","CDC19","PDC1","PDC5","PDC6",
        "ADH1","ADH2","ADH3","ADH4","ADH5","ADH6","GPD1","GPD2",
        "GDH1","GDH2","GDH3","LPD1","ACO1","CIT1","MDH1","MDH2","MDH3",
        "COR1","ATP1","ATP2","ATP3","ATP4","ATP5","ATP7","ATP14","ATP16",
        "TKL1","TKL2","TAL1","RKI1","SOL3","SOL4","ZWF1","GND1","GND2",
        "INO1","PYC1","PYC2","FUM1","SDH1","SDH2","IDH1","IDH2",
        "DAK1","DAK2","ERG10","ERG13",
    ]},
    # ── Chaperones / Protein folding ─────────────────────────────────────
    **{g: "Chaperones/Protein folding" for g in [
        "KAR2","PDI1","ERO1","EUG1","MPD1","MPD2","CPR1","CPR2","CPR3",
        "CPR5","CPR6","CPR7","FKB1","FKB2","SIL1","LHS1","SCJ1","JEM1",
        "SSA1","SSA2","SSA3","SSA4","SSB1","SSB2","SSC1","SSZ1","SSE1",
        "SSE2","STI1","HSP82","HSC82","AHA1","CNS1","CDC37","SBA1",
        "TCP1","CCT2","CCT3","CCT4","CCT5","CCT6","CCT7","CCT8",
        "HSP60","HSP10","MDJ1","GRP78","GRP94","HSP104","HSP26","HSP42",
        "HSP12","HSP150","BTT1","SHR3",
    ]},
    # ── ER / Golgi / Secretory pathway ───────────────────────────────────
    **{g: "ER/Golgi/Secretory" for g in [
        "SEC61","SEC62","SEC63","SEC11","SEC13","SEC16","SEC17","SEC18",
        "SEC23","SEC24","SEC31","SEC39","DSL1","DSL2","USE1",
        "SRP14","SRP21","SRP54","SRP68","SRP72","SRP101","SRP102",
        "SRB1","SRB2","SRB5",
        "OST1","OST2","OST3","OST4","OST5","OST6","WBP1","SWP1",
        "STT3","DAD1",
        "STE13","KEX1","KEX2","OPM1",
        "DPM1","SEC53","OCH1","MNS1","MNT1","MNT2","MNN9","VAN1",
        "ALG1","ALG2","ALG3","ALG5","ALG6","ALG8","ALG9","ALG10",
        "ALG11","ALG12","ALG13","ALG14",
        "PMI40","PMT1","PMT2","PMT3","PMT4","PMT5","PMT6",
        "IRE1","HAC1","PDI1",
        "CDC48","UFD1","NPL4","VCP",
        "ARL3","PIK1","SEC14","COP1","RET1","RET2","RET3",
        "ARF1","ARF2","GEA1","GEA2","SYS1","GOS1","SFT1","BET1",
        "LDB19","APM3","APL6","VPS21","VPS4",
    ]},
    # ── Amino acid metabolism ────────────────────────────────────────────
    **{g: "Amino acid metabolism" for g in [
        "MET6","CYS3","CYS4","ARG1","ARG2","ARG3","ARG4","ARG5",
        "ARG6","ARG7","ARG8","ARO1","ARO2","ARO3","ARO4","ARO7","ARO8",
        "LYS1","LYS2","LYS4","LYS9","LYS12","LYS14","LYS20","LYS21",
        "ILV1","ILV2","ILV3","ILV5","ILV6","LEU1","LEU2","LEU4","LEU5",
        "CAR1","CAR2","OAT1","PRO1","PRO2","PRO3","SER1","SER2","SER3",
        "SER33","TRP1","TRP2","TRP3","TRP4","TRP5","SAH1","MET2",
    ]},
    # ── Redox homeostasis ────────────────────────────────────────────────
    **{g: "Redox homeostasis" for g in [
        "TRR1","TRR2","TRX1","TRX2","TRX3","TSA1","TSA2","AHP1","PRX1",
        "GRX1","GRX2","GRX3","GRX4","GRX5","GLR1","GTT1","GTT2",
        "SOD1","SOD2","CCP1","CTT1","CTA1","TUM1",
    ]},
    # ── Proteolysis / Proteasome ─────────────────────────────────────────
    **{g: "Proteolysis/Proteasome" for g in [
        "PRE1","PRE2","PRE3","PRE4","PRE5","PRE6","PRE7","PRE8","PRE9",
        "PUP1","PUP2","PUP3","SCL1","RPN1","RPN2","RPN3","RPN5","RPN6",
        "RPN7","RPN8","RPN9","RPN10","RPN11","RPN12","RPN13",
        "RPT1","RPT2","RPT3","RPT4","RPT5","RPT6",
        "UBI4","RAD6","UBC1","UBC4","UBC5","UBC7","DOA4","RPN4",
        "APE1","APE3","CPY1","CPY2","PEP4","VPS10","PRB1","PRC1",
        "DUG1","DUG2","DUG3","GGT1",
    ]},
    # ── Nucleic acid metabolism ──────────────────────────────────────────
    **{g: "Nucleic acid metabolism" for g in [
        "ADO1","GUA1","ADE1","ADE2","ADE3","ADE4","ADE5","ADE6","ADE7",
        "ADE8","ADE12","ADE13","ADE16","ADE17","YNK1","GND1","GND2",
        "URA1","URA2","URA3","URA4","URA5","URA6","URA7","URA8","URA10",
        "CDD1","FMN1","RIB1","RIB2","RIB3","RIB4","RIB5","RIB7",
        "DBP6","RPO21","NOP1","NOP2","NOP56","NOP58",
    ]},
}

# Build reverse lookup: gene → GO Slim category
def get_go_slim(gene, desc):
    """Return GO Slim category for a gene symbol or fall back to description parsing."""
    g = str(gene).upper().strip()
    if g in GO_SLIM_MAP:
        return GO_SLIM_MAP[g]

    d = str(desc).lower()

    # Ribosomal / translation (catches numbered paralogs not in map)
    if re.search(r"ribosom|elongation factor|translation factor|eif\d|rpl\d|rps\d", d):
        return "Translation/Ribosomal"
    if re.search(r"glyceraldeh|enolase|phosphoglycerate|pyruvate kin|pyruvate dec|"
                 r"aldolase|hexokinase|dehydrogenase|transketol|citrate synth|"
                 r"oxidoreduct|reductase|alcohol dehyd", d):
        return "Energy metabolism"
    if re.search(r"chaperone|heat.shock|hsp\d|grp\d|cyclophilin|peptidyl.prolyl|"
                 r"fkbp|disulfide.isomerase|thioredoxin", d):
        return "Chaperones/Protein folding"
    if re.search(r"endoplasmic|golgi|sec\d|secretory|translocon|signal peptid|"
                 r"n.glycan|mannosyl|glycosyl|dolichol|vesicle|coatomer", d):
        return "ER/Golgi/Secretory"
    if re.search(r"aminotransfer|amino acid biosyn|methionine|cysteine synth|"
                 r"arginine|lysine|leucine|isoleucine|tryptophan", d):
        return "Amino acid metabolism"
    if re.search(r"peroxiredoxin|thioredoxin reduct|glutaredoxin|superoxide|"
                 r"catalase|glutathione", d):
        return "Redox homeostasis"
    if re.search(r"proteasome|ubiquitin|protease|peptidase", d):
        return "Proteolysis/Proteasome"
    if re.search(r"nucleotide|purine|pyrimidine|riboflavin|rna polymerase", d):
        return "Nucleic acid metabolism"
    if re.search(r"uncharacterized|hypothetical|unnamed|unknown", d):
        return "Uncharacterised"

    return "Other"

# Category display order and colours for stacked bars
CAT_ORDER = [
    "ER/Golgi/Secretory",
    "Nucleic acid metabolism",
    "Amino acid metabolism",
    "Proteolysis/Proteasome",
    "Redox homeostasis",
    "Chaperones/Protein folding",
    "Translation/Ribosomal",
    "Energy metabolism",
    "Other",
    "Uncharacterised",
]

CAT_COLOURS = {
    "Translation/Ribosomal":    PAL["Translation"],
    "Energy metabolism":         PAL["Energy"],
    "Chaperones/Protein folding": PAL["Chaperone"],
    "Redox homeostasis":         PAL["Redox"],
    "Proteolysis/Proteasome":    PAL["Proteolysis"],
    "Amino acid metabolism":     PAL["AminoAcid"],
    "ER/Golgi/Secretory":        PAL["ERGolgi"],
    "Nucleic acid metabolism":   PAL["NucleicAcid"],
    "Uncharacterised":           "#A9CCE3",
    "Other":                     PAL["Other"],
}

# =============================================================================
# 2.  PREPROCESSING PIPELINE
# =============================================================================

def load_and_preprocess(path):
    """
    Two-step normalisation pipeline:

    Step 1 - Protein loading correction (ug-based):
        intensity_corrected = raw_intensity / ug_loaded * REFERENCE_UG
        Corrects for unequal protein amounts loaded per fraction.
        BfA microsomes had 250 ug vs 98 ug for control microsomes (2.55x
        difference); without this step all BfA microsomal proteins appear
        falsely enriched purely due to loading.

    Step 2 - Median normalisation:
        Scales each fraction so its median equals the global median.
        Corrects for residual MS run-to-run variation and injection differences.

    Returns:
        df_norm  -- gene-level protein table (fully normalised)
        df_raw   -- gene-level protein table (raw, for QC rank plots)
        df_step1 -- gene-level protein table (after Step 1 only)
    """
    raw = pd.read_excel(path)
    raw.columns = raw.columns.str.strip()

    def _gene(desc):
        m = re.search(r"GN=(\S+)", str(desc))
        return m.group(1).upper() if m else ""

    def _prot_name(desc):
        m = re.match(r"^(.+?)\s+OS=", str(desc))
        return m.group(1).strip() if m else str(desc)[:80]

    raw["gene"]      = raw["protein.Description"].apply(_gene)
    raw["prot_name"] = raw["protein.Description"].apply(_prot_name)
    raw["go_slim"]   = raw.apply(
        lambda r: get_go_slim(r["gene"], r["protein.Description"]), axis=1)
    raw["gene_id"]   = raw.apply(
        lambda r: r["gene"] if r["gene"] else r["protein.Accession"], axis=1)

    agg_dict = {c: "sum" for c in QUANT_COLS}
    agg_dict.update({"prot_name": "first", "go_slim": "first",
                     "protein.Description": "first",
                     "protein.Accession":   "first"})
    df = raw.groupby("gene_id").agg(agg_dict).reset_index()
    df = df.rename(columns={"gene_id": "gene"})
    df_raw = df.copy()

    print(f"  Raw rows:            {len(raw)}")
    print(f"  After gene rollup:   {len(df)} unique gene-level entries")

    # === STEP 1: Protein loading correction ================================
    print(f"\n  Step 1: Protein loading correction (reference = {REFERENCE_UG} ug)")
    print(f"  {'Fraction':25s}  {'ug loaded':>10}  {'scale factor':>12}")
    print(f"  {'-'*52}")

    df_step1 = df.copy()
    for c in QUANT_COLS:
        ug = PROTEIN_LOADED_UG[c]
        sf = REFERENCE_UG / ug
        df_step1[c] = df[c].apply(lambda x, sf=sf: x * sf if x > 0 else 0)
        print(f"  {c:25s}  {ug:>10.0f}  {sf:>12.3f}x")

    # === STEP 2: Median normalisation ======================================
    print(f"\n  Step 2: Median normalisation (run-level technical correction)")
    global_median = np.median(
        [v for c in QUANT_COLS for v in df_step1[c].values if v > 0])
    print(f"  Global median after Step 1: {global_median:,.0f}")
    print(f"  {'Fraction':25s}  {'pre-median':>12}  {'post-median':>12}")
    print(f"  {'-'*55}")

    df_norm = df_step1.copy()
    for c in QUANT_COLS:
        col_med = df_step1[df_step1[c] > 0][c].median()
        if col_med > 0:
            df_norm[c] = df_step1[c].apply(
                lambda x, cm=col_med: x / cm * global_median if x > 0 else 0)
        m_after = df_norm[df_norm[c] > 0][c].median()
        print(f"  {c:25s}  {col_med:>12,.0f}  {m_after:>12,.0f}")

    bfa_med  = df_norm[df_norm["BfA_microsomes"] > 0]["BfA_microsomes"].median()
    ctrl_med = df_norm[df_norm["microsomes"]     > 0]["microsomes"].median()
    print(f"\n  Sanity check - normalised microsome medians:")
    print(f"    BfA  microsomes: {bfa_med:>10,.0f}")
    print(f"    Ctrl microsomes: {ctrl_med:>10,.0f}")
    print(f"    Ratio          : {bfa_med/ctrl_med:.3f}x  (target ~ 1.0)")

    return df_norm, df_raw, df_step1


# =============================================================================
# 3.  VENN HELPER
# =============================================================================

def draw_venn3_manual(ax, sets, labels, colours, title):
    """Pure matplotlib 3-way Venn (no external library required)."""
    sA, sB, sC = sets
    lA, lB, lC = labels
    cA, cB, cC = colours

    only_A  = sA - sB - sC
    only_B  = sB - sA - sC
    only_C  = sC - sA - sB
    AB      = (sA & sB) - sC
    AC      = (sA & sC) - sB
    BC      = (sB & sC) - sA
    ABC     = sA & sB & sC

    # Circle centres
    from matplotlib.patches import Circle
    cx = [0.38, 0.62, 0.50]
    cy = [0.55, 0.55, 0.32]
    r  = 0.28

    for (x, y, c) in zip(cx, cy, [cA, cB, cC]):
        ax.add_patch(Circle((x, y), r, color=c, alpha=0.30,
                             transform=ax.transAxes, clip_on=False))
        ax.add_patch(Circle((x, y), r, fill=False, edgecolor=c,
                             linewidth=1.6, transform=ax.transAxes, clip_on=False))

    tw = dict(ha="center", va="center", fontsize=9, fontweight="bold",
              transform=ax.transAxes)
    ax.text(0.26, 0.64, str(len(only_A)),  **tw)
    ax.text(0.74, 0.64, str(len(only_B)),  **tw)
    ax.text(0.50, 0.18, str(len(only_C)),  **tw)
    ax.text(0.50, 0.64, str(len(AB)),      **tw)
    ax.text(0.33, 0.36, str(len(AC)),      **tw)
    ax.text(0.67, 0.36, str(len(BC)),      **tw)
    ax.text(0.50, 0.45, str(len(ABC)),     **tw)

    for (x, y, c, lbl, tot) in [
        (0.14, 0.62, cA, lA, len(sA)),
        (0.86, 0.62, cB, lB, len(sB)),
        (0.50, 0.03, cC, lC, len(sC)),
    ]:
        ax.text(x, y, f"{lbl}\n(n={tot})", ha="center", va="center",
                fontsize=8, fontweight="bold", color=c,
                transform=ax.transAxes)

    ax.set_title(title, fontsize=9.5, fontweight="bold", pad=4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


# =============================================================================
# 4.  SECRETORY PATHWAY TABLE DATA
#     Expanded: includes Sec61 translocon, SRP, OST complex, COPII, Golgi
# =============================================================================

SECRETORY_PROTEINS = [
    # (gene,  display_name,           group)
    # ── ER translocation ──────────────────
    ("SEC61",  "SEC61 (translocon α)",    "ER translocon"),
    ("SEC62",  "SEC62 (translocon)",      "ER translocon"),
    ("SEC63",  "SEC63 (translocon)",      "ER translocon"),
    ("SEC11",  "SEC11 (signal peptidase)","ER translocon"),
    # ── SRP pathway ───────────────────────
    ("SRP54",  "SRP54",                   "SRP pathway"),
    ("SRP14",  "SRP14",                   "SRP pathway"),
    ("SRP68",  "SRP68",                   "SRP pathway"),
    # ── ER chaperones ─────────────────────
    ("KAR2",   "KAR2 (BiP/Hsp70)",       "ER chaperone"),
    ("PDI1",   "PDI1 (Protein disulfide isomerase)", "ER chaperone"),
    ("ERO1",   "ERO1",                    "ER chaperone"),
    ("CPR1",   "CPR1 (Cyclophilin)",      "ER chaperone"),
    # ── OST complex ───────────────────────
    ("STT3",   "STT3 (OST catalytic)",    "N-glycosylation"),
    ("OST1",   "OST1",                    "N-glycosylation"),
    ("WBP1",   "WBP1",                    "N-glycosylation"),
    # ── Glycan biosynthesis ───────────────
    ("DPM1",   "DPM1 (Dol-P-Man synth)", "N-glycosylation"),
    ("SEC53",  "SEC53 (Phosphomannomut)", "N-glycosylation"),
    ("OCH1",   "OCH1 (α-mannosyltransf)","N-glycosylation"),
    ("MNS1",   "MNS1 (α-mannosidase)",   "N-glycosylation"),
    # ── ERAD ──────────────────────────────
    ("CDC48",  "CDC48 (p97/VCP)",         "ERAD"),
    ("IRE1",   "IRE1 (UPR sensor)",       "ERAD"),
    # ── COPII vesicles ────────────────────
    ("SEC23",  "SEC23 (COPII inner coat)","COPII vesicles"),
    ("SEC24",  "SEC24 (COPII inner coat)","COPII vesicles"),
    ("SEC13",  "SEC13 (COPII outer coat)","COPII vesicles"),
    # ── Golgi / vesicle trafficking ───────
    ("ARL3",   "ARL3 (Golgi GTPase)",     "Golgi/trafficking"),
    ("PIK1",   "PIK1 (PI4-kinase)",       "Golgi/trafficking"),
    ("STE13",  "STE13 (DPAP-A)",          "Golgi/trafficking"),
    ("DSL1",   "DSL1 (ER tether)",        "Golgi/trafficking"),
    ("SEC39",  "SEC39 (DSL complex)",     "Golgi/trafficking"),
    # ── NAC complex ───────────────────────
    ("EGD1",   "EGD1 (NACβ)",            "Co-translational"),
    # ── Cytosolic chaperones ──────────────
    ("SSZ1",   "SSZ1 (Hsp70/RAC)",       "Cytosolic chaperone"),
    ("SSB2",   "SSB2 (Hsp70/RAC)",       "Cytosolic chaperone"),
    ("SSA1",   "SSA1 (Hsp70)",           "Cytosolic chaperone"),
    ("SSA3",   "SSA3 (Hsp70)",           "Cytosolic chaperone"),
    ("SSC1",   "SSC1 (mito Hsp70)",      "Cytosolic chaperone"),
    ("SSE1",   "SSE1 (Hsp110)",          "Cytosolic chaperone"),
    ("HSP82",  "HSP82 (Hsp90)",          "Cytosolic chaperone"),
    ("STI1",   "STI1 (Hop co-chaperone)","Cytosolic chaperone"),
    ("TCP1",   "TCP1 (TRiC/CCT)",        "Cytosolic chaperone"),
]

GROUP_COLOURS = {
    "ER translocon":      "#C0392B",
    "SRP pathway":        "#E74C3C",
    "ER chaperone":       "#E67E22",
    "N-glycosylation":    "#F48FB1",
    "ERAD":               "#9B59B6",
    "COPII vesicles":     "#2980B9",
    "Golgi/trafficking":  "#1ABC9C",
    "Co-translational":   "#3498DB",
    "Cytosolic chaperone":"#27AE60",
}

GROUP_ORDER = list(GROUP_COLOURS.keys())


# =============================================================================
# 5.  PANEL LETTER HELPER
# =============================================================================

def add_panel_letter(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="top", ha="left")


# =============================================================================
# 6.  QC FIGURE  (4 panels)
# =============================================================================

def make_qc_figure(df, df_raw, df_step1):

    fig = plt.figure(figsize=(16, 13))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(3, 2, figure=fig,
                           hspace=0.50, wspace=0.35,
                           left=0.09, right=0.96,
                           top=0.93, bottom=0.08)

    fig.suptitle("Quality Control – Untargeted Proteomics\n"
                 "K. phaffii Cell-Free Lysate Fractions  "
                 "(two-step normalisation: µg loading + median)",
                 fontsize=11, fontweight="bold", y=0.97)

    # ── QC-A: Log2 intensity distributions ──────────────────────────────
    ax_qa = fig.add_subplot(gs[0, 0])
    add_panel_letter(ax_qa, "A")

    cmap_fracs = plt.cm.tab10(np.linspace(0, 0.6, 6))
    for i, col in enumerate(QUANT_COLS):
        vals = np.log2(df[df[col] > 0][col].values)
        ax_qa.hist(vals, bins=30, alpha=0.55, color=cmap_fracs[i],
                   label=FRACTION_LABELS[col].replace("\n", " "),
                   density=True, linewidth=0)

    ax_qa.set_xlabel("log₂(Normalised Intensity)", fontsize=9)
    ax_qa.set_ylabel("Density", fontsize=9)
    ax_qa.set_title("Intensity Distribution\n(Median-Normalised)", fontsize=9.5,
                     fontweight="bold")
    ax_qa.legend(fontsize=6.5, ncol=2, framealpha=0.7)

    # ── QC-B: Protein detection per fraction (bar chart) ────────────────
    ax_qb = fig.add_subplot(gs[0, 1])
    add_panel_letter(ax_qb, "B")

    counts = [(FRACTION_LABELS[c].replace("\n", " "),
               (df[c] > 0).sum(),
               PAL["bfa"] if "BfA" in c else PAL["ctrl"])
              for c in QUANT_COLS]
    labels_b, heights_b, colours_b = zip(*counts)

    bars = ax_qb.bar(range(6), heights_b, color=colours_b,
                      edgecolor="white", linewidth=0.5, width=0.7)
    for bar, h in zip(bars, heights_b):
        ax_qb.text(bar.get_x() + bar.get_width() / 2, h + 3, str(h),
                   ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax_qb.set_xticks(range(6))
    ax_qb.set_xticklabels(labels_b, fontsize=7.5, rotation=25, ha="right")
    ax_qb.set_ylabel("Proteins detected (intensity > 0)", fontsize=9)
    ax_qb.set_title("Proteins Detected per Fraction", fontsize=9.5,
                     fontweight="bold")

    handles_b = [mpatches.Patch(color=PAL["bfa"],  label="BfA-treated"),
                 mpatches.Patch(color=PAL["ctrl"], label="Control")]
    ax_qb.legend(handles=handles_b, fontsize=7.5, framealpha=0.7)
    ax_qb.set_ylim(0, max(heights_b) * 1.15)

    # ── QC-C: Inter-fraction Pearson correlation (log2) ─────────────────
    # ── QC-C: Protein loading per fraction (bar + annotation) ───────────
    ax_qc_load = fig.add_subplot(gs[1, 0])
    add_panel_letter(ax_qc_load, "C")

    fracs_load  = list(PROTEIN_LOADED_UG.keys())
    ug_vals     = list(PROTEIN_LOADED_UG.values())
    colours_load = [PAL["bfa"] if "BfA" in f else PAL["ctrl"] for f in fracs_load]
    xlabels_load = [FRACTION_LABELS[f].replace("\n", " ") for f in fracs_load]

    bars_l = ax_qc_load.bar(range(6), ug_vals, color=colours_load,
                              edgecolor="white", linewidth=0.5, width=0.7)
    for bar, ug in zip(bars_l, ug_vals):
        ax_qc_load.text(bar.get_x() + bar.get_width()/2, ug + 2,
                         f"{ug} µg", ha="center", va="bottom",
                         fontsize=8.5, fontweight="bold")

    # Annotate the key ratio
    bfa_mic_ug  = PROTEIN_LOADED_UG["BfA_microsomes"]
    ctrl_mic_ug = PROTEIN_LOADED_UG["microsomes"]
    ax_qc_load.annotate(
        f"2.55× more protein\nloaded (BfA vs ctrl)",
        xy=(2, bfa_mic_ug), xytext=(3.5, bfa_mic_ug - 30),
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.4),
        fontsize=8, color="#C0392B", fontweight="bold",
        ha="center")

    ax_qc_load.set_xticks(range(6))
    ax_qc_load.set_xticklabels(xlabels_load, fontsize=7.5, rotation=25, ha="right")
    ax_qc_load.set_ylabel("Protein loaded (µg)", fontsize=9)
    ax_qc_load.set_title("Protein Amount Loaded per Fraction\n"
                          "(Step 1: used to derive loading scale factors)",
                          fontsize=9.5, fontweight="bold")
    ax_qc_load.set_ylim(0, max(ug_vals) * 1.25)
    handles_l = [mpatches.Patch(color=PAL["bfa"],  label="BfA-treated"),
                 mpatches.Patch(color=PAL["ctrl"], label="Control")]
    ax_qc_load.legend(handles=handles_l, fontsize=7.5, framealpha=0.7)

    # ── QC-D: Before vs after normalisation — microsome medians ──────────
    ax_qc_ba = fig.add_subplot(gs[1, 1])
    add_panel_letter(ax_qc_ba, "D")

    mic_pairs = [("BfA_microsomes", "BfA"), ("microsomes", "Ctrl")]
    stages    = ["Raw", "After\nStep 1\n(µg norm)", "After\nStep 2\n(+ median)"]
    x_pos_ba  = np.arange(len(stages))

    for label, (col, grp_label) in enumerate(mic_pairs):
        raw_med    = df_raw[df_raw[col] > 0][col].median()
        step1_med  = df_step1[df_step1[col] > 0][col].median()
        norm_med   = df[df[col] > 0][col].median()
        medians    = [raw_med, step1_med, norm_med]
        colour     = PAL["bfa"] if grp_label == "BfA" else PAL["ctrl"]
        offset     = -0.18 if grp_label == "BfA" else 0.18
        ax_qc_ba.bar(x_pos_ba + offset, medians, width=0.34,
                      color=colour, alpha=0.85,
                      edgecolor="white", linewidth=0.4,
                      label=f"{grp_label} Microsomes")
        for xi, m in zip(x_pos_ba + offset, medians):
            ax_qc_ba.text(xi, m + max(medians)*0.01,
                           f"{m/1e6:.1f}M", ha="center", va="bottom",
                           fontsize=7, fontweight="bold",
                           color=colour)

    ax_qc_ba.set_xticks(x_pos_ba)
    ax_qc_ba.set_xticklabels(stages, fontsize=8)
    ax_qc_ba.set_ylabel("Median protein intensity", fontsize=9)
    ax_qc_ba.set_title("Effect of Normalisation on Microsome Fractions\n"
                        "(median intensity at each step)",
                        fontsize=9.5, fontweight="bold")
    ax_qc_ba.legend(fontsize=7.5, framealpha=0.7)

    # Add expected-equal annotation after Step 2
    ax_qc_ba.annotate("Medians equalised\n(target of Step 2)",
                       xy=(2, ax_qc_ba.get_ylim()[1] * 0.6),
                       fontsize=7.5, color="#555555", ha="center",
                       style="italic")

    ax_qc = fig.add_subplot(gs[2, 0])
    add_panel_letter(ax_qc, "E")

    log_mat = np.log2(df[QUANT_COLS].replace(0, np.nan))
    corr = log_mat.corr(method="pearson", min_periods=20)
    corr.index   = [FRACTION_LABELS[c].replace("\n", " ") for c in QUANT_COLS]
    corr.columns = [FRACTION_LABELS[c].replace("\n", " ") for c in QUANT_COLS]

    mask = np.triu(np.ones_like(corr, dtype=bool))
    cmap_corr = LinearSegmentedColormap.from_list(
        "corr", ["#4C78A8", "#FAFAFA", "#E45756"], N=256)

    sns.heatmap(corr, ax=ax_qc, cmap=cmap_corr,
                vmin=0.3, vmax=1.0,
                annot=True, fmt=".2f", annot_kws={"size": 7},
                linewidths=0.5, linecolor="#DDDDDD",
                mask=mask,
                cbar_kws={"label": "Pearson r", "shrink": 0.75})
    ax_qc.set_title("Inter-fraction Correlation\n(Pearson r, log₂ intensity)",
                     fontsize=9.5, fontweight="bold")
    ax_qc.tick_params(axis="x", rotation=30, labelsize=7.5)
    ax_qc.tick_params(axis="y", rotation=0,  labelsize=7.5)

    # ── QC-D: Protein rank abundance (all 6 fractions) ──────────────────
    ax_qd = fig.add_subplot(gs[2, 1])
    add_panel_letter(ax_qd, "F")

    for i, col in enumerate(QUANT_COLS):
        vals = np.sort(df[df[col] > 0][col].values)[::-1]
        ax_qd.plot(np.arange(1, len(vals) + 1),
                   np.log10(vals),
                   color=cmap_fracs[i], linewidth=1.5, alpha=0.85,
                   label=FRACTION_LABELS[col].replace("\n", " "))

    ax_qd.set_xlabel("Protein rank (by abundance)", fontsize=9)
    ax_qd.set_ylabel("log₁₀(Intensity)", fontsize=9)
    ax_qd.set_title("Protein Rank–Abundance Plot\n(Dynamic range assessment)",
                     fontsize=9.5, fontweight="bold")
    ax_qd.legend(fontsize=6.5, ncol=2, framealpha=0.7)

    # Single-replicate note
    fig.text(0.5, 0.005,
             "Note: Data represent single biological replicates per condition. "
             "No CV or statistical testing between conditions is possible without biological replication.",
             ha="center", va="bottom", fontsize=7, color="#888888",
             style="italic")

    return fig


# =============================================================================
# 7.  MAIN FIGURE  (5 panels A–E)
# =============================================================================

def make_main_figure(df, df_raw):

    # ── Pre-compute sets ──────────────────────────────────────────────────
    ctrl_sets = {c: set(df.loc[df[c] > 0, "gene"]) for c in CTRL_COLS}
    sL = ctrl_sets["lysate"]
    sS = ctrl_sets["supernatant"]
    sM = ctrl_sets["microsomes"]
    total_ctrl = len(sL | sS | sM)

    # ── Functional composition (% normalised intensity) ────────────────
    frac_pct = {}
    for label, col in [("Lysate",    "lysate"),
                       ("Supernatant", "supernatant"),
                       ("Microsomes",  "microsomes")]:
        detected = df[df[col] > 0]
        total    = detected[col].sum()
        pct = {}
        for cat in CAT_ORDER:
            s = detected.loc[detected["go_slim"] == cat, col].sum()
            pct[cat] = s / total * 100 if total > 0 else 0
        frac_pct[label] = pct

    # ── Top 15 proteins ───────────────────────────────────────────────
    def top15(col):
        sub = df[df[col] > 0].sort_values(col, ascending=False).head(15)
        return sub[["gene", "go_slim", col]].copy()

    top_lys  = top15("lysate")
    top_mic  = top15("microsomes")

    # ── Secretory pathway presence ────────────────────────────────────
    sec_rows = []
    for gene, display, group in SECRETORY_PROTEINS:
        g_upper = gene.upper()
        match = df[df["gene"].str.upper() == g_upper]
        if len(match) > 0:
            row = match.iloc[0]
            in_lys = row["lysate"] > 0
            in_mic = row["microsomes"] > 0
            in_sup = row["supernatant"] > 0
        else:
            in_lys = in_mic = in_sup = False
        sec_rows.append({
            "gene": gene, "display": display, "group": group,
            "in_lysate": in_lys, "in_microsomes": in_mic,
            "in_supernatant": in_sup,
        })
    sec_df = pd.DataFrame(sec_rows)
    sec_df["group_order"] = sec_df["group"].map(
        {g: i for i, g in enumerate(GROUP_ORDER)})
    sec_df = sec_df.sort_values("group_order").reset_index(drop=True)

    # ── Figure layout ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 20))
    fig.patch.set_facecolor("white")

    outer = gridspec.GridSpec(
        3, 2, figure=fig,
        height_ratios=[1, 1, 1.55],
        hspace=0.38, wspace=0.30,
        left=0.08, right=0.97,
        top=0.95, bottom=0.03,
    )

    fig.suptitle(
        "Proteome Characterisation of Control K. phaffii Cell-Free Lysate Fractions",
        fontsize=11, fontweight="bold", y=0.975, color="#2B2B2B")

    # =======================================================================
    # PANEL A – 3-way Venn
    # =======================================================================
    ax_a = fig.add_subplot(outer[0, 0])
    add_panel_letter(ax_a, "A")

    gL = sL - sS - sM
    gS = sS - sL - sM
    gM = sM - sL - sS

    if HAS_VENN:
        gLS  = (sL & sS) - sM
        gLM  = (sL & sM) - sS
        gSM  = (sS & sM) - sL
        gLSM = sL & sS & sM
        v = venn3(
            subsets=(len(gL), len(gS), len(gLS),
                     len(gM), len(gLM), len(gSM), len(gLSM)),
            set_labels=("", "", ""), ax=ax_a)
        patch_colours = {
            "100": "#E74C3C", "010": "#3498DB", "001": "#27AE60",
            "110": "#AF7AC5", "101": "#F39C12", "011": "#17A589",
            "111": "#ABB2B9"}
        for pid, c in patch_colours.items():
            p = v.get_patch_by_id(pid)
            if p: p.set_color(c); p.set_alpha(0.38)
        for t in v.subset_labels:
            if t: t.set_fontsize(9); t.set_fontweight("bold")
        ax_a.text(-0.48, 0.50, f"Lysate\n(n={len(sL)})",
                  fontsize=9, fontweight="bold", color="#E74C3C",
                  ha="center", transform=ax_a.transAxes)
        ax_a.text(1.48, 0.50, f"Supernatant\n(n={len(sS)})",
                  fontsize=9, fontweight="bold", color="#3498DB",
                  ha="center", transform=ax_a.transAxes)
        ax_a.text(0.50, -0.10, f"Microsomes\n(n={len(sM)})",
                  fontsize=9, fontweight="bold", color="#27AE60",
                  ha="center", transform=ax_a.transAxes)
    else:
        draw_venn3_manual(
            ax_a, (sL, sS, sM),
            (f"Lysate", f"Supernatant", f"Microsomes"),
            ("#E74C3C", "#3498DB", "#27AE60"),
            "")

    ax_a.set_title(f"Protein Detection Overlap (Control)\n"
                   f"Total detected: {total_ctrl} proteins",
                   fontsize=9.5, fontweight="bold", pad=5)

    # =======================================================================
    # PANEL B – Stacked bar (% intensity by GO Slim)
    # =======================================================================
    ax_b = fig.add_subplot(outer[0, 1])
    add_panel_letter(ax_b, "B")

    fracs_b = ["Lysate", "Supernatant", "Microsomes"]
    x_pos   = np.arange(3)
    bottoms = {f: 0 for f in fracs_b}

    for cat in CAT_ORDER:
        heights = [frac_pct[f].get(cat, 0) for f in fracs_b]
        bot     = [bottoms[f] for f in fracs_b]
        c       = CAT_COLOURS.get(cat, PAL["Other"])
        ax_b.bar(x_pos, heights, 0.60, bottom=bot,
                 color=c, edgecolor="white", linewidth=0.3, label=cat)
        for i, f in enumerate(fracs_b):
            if heights[i] >= 5:
                ax_b.text(i, bottoms[f] + heights[i] / 2,
                          f"{heights[i]:.0f}%",
                          ha="center", va="center",
                          fontsize=7, color="white", fontweight="bold")
            bottoms[f] += heights[i]

    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(fracs_b, fontsize=9.5)
    ax_b.set_ylabel("% of Total Fraction Intensity\n(Median-normalised)", fontsize=8.5)
    ax_b.set_ylim(0, 103)
    ax_b.set_title("Functional Composition by Fraction\n(GO Slim annotation)",
                   fontsize=9.5, fontweight="bold", pad=5)

    handles_b = [mpatches.Patch(color=CAT_COLOURS.get(c, PAL["Other"]), label=c)
                 for c in reversed(CAT_ORDER)
                 if any(frac_pct[f].get(c, 0) > 0.1 for f in fracs_b)]
    ax_b.legend(handles=handles_b, fontsize=6.2, loc="upper left",
                bbox_to_anchor=(1.01, 1.0), frameon=True,
                fancybox=False, edgecolor="#CCCCCC",
                title="GO Slim category", title_fontsize=6.5)

    # =======================================================================
    # PANEL C – Top 15 proteins in Control Lysate
    # =======================================================================
    ax_c = fig.add_subplot(outer[1, 0])
    add_panel_letter(ax_c, "C")

    top_lys_s = top_lys.sort_values("lysate", ascending=True)
    cols_c    = [CAT_COLOURS.get(r["go_slim"], PAL["Other"])
                 for _, r in top_lys_s.iterrows()]
    ax_c.barh(range(15), top_lys_s["lysate"].values / 1e6,
              color=cols_c, edgecolor="white", linewidth=0.3, height=0.72)
    ax_c.set_yticks(range(15))
    ax_c.set_yticklabels(top_lys_s["gene"].values,
                          fontsize=8.5, fontstyle="italic")
    ax_c.set_xlabel("Intensity (×10⁶, normalised)", fontsize=8.5)
    ax_c.set_title("Top 15 Proteins — Control Lysate",
                   fontsize=9.5, fontweight="bold", pad=5)

    # Colour legend (only categories present)
    cats_c = top_lys_s["go_slim"].unique()
    leg_c  = [mpatches.Patch(color=CAT_COLOURS.get(c, PAL["Other"]), label=c)
              for c in cats_c]
    ax_c.legend(handles=leg_c, fontsize=6, loc="lower right",
                framealpha=0.8, edgecolor="#CCCCCC")

    # =======================================================================
    # PANEL D – Top 15 proteins in Control Microsomes
    # =======================================================================
    ax_d = fig.add_subplot(outer[1, 1])
    add_panel_letter(ax_d, "D")

    top_mic_s = top_mic.sort_values("microsomes", ascending=True)
    cols_d    = [CAT_COLOURS.get(r["go_slim"], PAL["Other"])
                 for _, r in top_mic_s.iterrows()]
    ax_d.barh(range(15), top_mic_s["microsomes"].values / 1e6,
              color=cols_d, edgecolor="white", linewidth=0.3, height=0.72)
    ax_d.set_yticks(range(15))
    ax_d.set_yticklabels(top_mic_s["gene"].values,
                          fontsize=8.5, fontstyle="italic")
    ax_d.set_xlabel("Intensity (×10⁶, normalised)", fontsize=8.5)
    ax_d.set_title("Top 15 Proteins — Control Microsomes",
                   fontsize=9.5, fontweight="bold", pad=5)

    cats_d = top_mic_s["go_slim"].unique()
    leg_d  = [mpatches.Patch(color=CAT_COLOURS.get(c, PAL["Other"]), label=c)
              for c in cats_d]
    ax_d.legend(handles=leg_d, fontsize=6, loc="lower right",
                framealpha=0.8, edgecolor="#CCCCCC")

    # =======================================================================
    # PANEL E – Secretory pathway detection table
    #           Built as a proper subplot grid (not manual coord drawing)
    # =======================================================================
    ax_e_outer = outer[2, :]
    n_rows = len(sec_df)
    n_extra = 4   # header + 2 legend rows + padding

    # inner grid: n_rows data rows + n_extra extra rows, 4 columns
    gs_e = gridspec.GridSpecFromSubplotSpec(
        n_rows + n_extra, 4,
        subplot_spec=ax_e_outer,
        hspace=0.0, wspace=0.0,
        width_ratios=[0.04, 0.40, 0.32, 0.12, 0.12][1:],  # 3 data cols
    )

    # Column widths as ratios (colour strip excluded)
    # We draw colour strip as a thin column via the leftmost subplot
    gs_e2 = gridspec.GridSpecFromSubplotSpec(
        n_rows + n_extra, 4,
        subplot_spec=ax_e_outer,
        hspace=0.0, wspace=0.02,
    )

    # Use a single hidden axis for the whole panel
    ax_e = fig.add_subplot(outer[2, :])
    ax_e.set_xlim(0, 1)
    ax_e.set_ylim(0, 1)
    ax_e.axis("off")
    add_panel_letter(ax_e, "E", x=-0.03, y=1.02)
    ax_e.set_title("Secretory Pathway Protein Detection — Control Condition\n"
                   "(expanded panel: translocon, SRP, OST complex, COPII, Golgi)",
                   fontsize=9.5, fontweight="bold", pad=6)

    # Layout parameters (in axes [0,1] × [0,1])
    row_h   = 0.88 / (n_rows + 3)      # leave room for header + legend
    header_y = 0.96
    col_x   = [0.04, 0.44, 0.76, 0.88]  # left edges: strip|protein|group|Lys|Mic
    col_w   = [0.40, 0.32, 0.12, 0.12]

    # Header
    header_labels = ["Protein", "Functional Group", "Lys", "Mic"]
    for x, w, lbl in zip(col_x, col_w, header_labels):
        ax_e.add_patch(mpatches.FancyBboxPatch(
            (x, header_y - row_h), w, row_h,
            boxstyle="square,pad=0", fc="#2C3E50", ec="white", lw=0.4,
            transform=ax_e.transAxes, clip_on=False))
        ax_e.text(x + w / 2, header_y - row_h / 2, lbl,
                  ha="center", va="center", fontsize=8.5,
                  fontweight="bold", color="white",
                  transform=ax_e.transAxes)

    prev_group = None
    for idx, row in sec_df.iterrows():
        y_top = header_y - (idx + 2) * row_h

        # Background
        if row["in_lysate"] and row["in_microsomes"]:
            bg = "#E8F5E9"
        elif row["in_lysate"] or row["in_microsomes"]:
            bg = "#FFF8E1"
        elif row["in_supernatant"]:
            bg = "#FFF3E0"
        else:
            bg = "#FDECEA"

        for x, w in zip(col_x, col_w):
            ax_e.add_patch(mpatches.FancyBboxPatch(
                (x, y_top), w, row_h,
                boxstyle="square,pad=0", fc=bg, ec="#E0E0E0", lw=0.25,
                transform=ax_e.transAxes, clip_on=False))

        # Group separator
        if row["group"] != prev_group and prev_group is not None:
            ax_e.plot([col_x[0], col_x[-1] + col_w[-1]],
                      [y_top + row_h, y_top + row_h],
                      color="#888888", lw=0.7,
                      transform=ax_e.transAxes, clip_on=False)

        # Colour strip (left of protein name)
        gc = GROUP_COLOURS.get(row["group"], "#AAAAAA")
        ax_e.add_patch(mpatches.FancyBboxPatch(
            (0.005, y_top + 0.01 * row_h), 0.028, row_h * 0.98,
            boxstyle="square,pad=0", fc=gc, ec="none",
            transform=ax_e.transAxes, clip_on=False))

        # Protein name
        ax_e.text(col_x[0] + 0.012, y_top + row_h / 2,
                  row["display"], ha="left", va="center",
                  fontsize=7.5, fontstyle="italic", fontweight="bold",
                  transform=ax_e.transAxes)

        # Functional group
        ax_e.text(col_x[1] + col_w[1] / 2, y_top + row_h / 2,
                  row["group"], ha="center", va="center",
                  fontsize=7, transform=ax_e.transAxes)

        # Detected symbols
        for col_name, x_sym in [("in_lysate", col_x[2]),
                                  ("in_microsomes", col_x[3])]:
            detected = row[col_name]
            in_sup   = row["in_supernatant"]
            cx_sym   = x_sym + col_w[2] / 2

            if detected:
                ax_e.text(cx_sym, y_top + row_h / 2, "●",
                          ha="center", va="center", fontsize=11,
                          color="#27AE60", transform=ax_e.transAxes)
            elif in_sup:
                ax_e.text(cx_sym, y_top + row_h / 2, "○",
                          ha="center", va="center", fontsize=10,
                          color="#F39C12", transform=ax_e.transAxes)
            else:
                ax_e.text(cx_sym, y_top + row_h / 2, "✗",
                          ha="center", va="center", fontsize=9,
                          color="#E74C3C", transform=ax_e.transAxes)

        prev_group = row["group"]

    # Detection legend
    leg_y = header_y - (n_rows + 2.6) * row_h
    det_items = [
        ("●", "#27AE60", "Detected"),
        ("○", "#F39C12", "Supernatant only"),
        ("✗", "#E74C3C", "Not detected"),
    ]
    for i, (sym, c, lbl) in enumerate(det_items):
        x0 = 0.04 + i * 0.28
        ax_e.text(x0, leg_y, sym, ha="left", va="center",
                  fontsize=11, color=c, transform=ax_e.transAxes)
        ax_e.text(x0 + 0.04, leg_y, lbl, ha="left", va="center",
                  fontsize=7.5, transform=ax_e.transAxes)

    # Group colour legend
    leg_y2 = leg_y - row_h * 1.2
    for i, (g, c) in enumerate(GROUP_COLOURS.items()):
        col_i = i % 5
        row_i = i // 5
        x0 = 0.04 + col_i * 0.19
        y0 = leg_y2 - row_i * row_h * 1.4
        ax_e.add_patch(mpatches.FancyBboxPatch(
            (x0, y0 - row_h * 0.4), 0.018, row_h * 0.7,
            boxstyle="square,pad=0", fc=c, ec="none",
            transform=ax_e.transAxes, clip_on=False))
        ax_e.text(x0 + 0.024, y0 - row_h * 0.05, g,
                  ha="left", va="center", fontsize=6.5,
                  transform=ax_e.transAxes)

    return fig


# =============================================================================
# 8.  RUN
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Loading and preprocessing data …")
    df_norm, df_raw, df_step1 = load_and_preprocess(DATA_PATH)

    print(f"\nProtein counts per fraction (normalised, intensity > 0):")
    for c in QUANT_COLS:
        print(f"  {c:25s}: {(df_norm[c] > 0).sum()}")

    print("\nGO Slim category breakdown (all detected):")
    detected_any = df_norm[(df_norm[QUANT_COLS] > 0).any(axis=1)]
    for cat, n in detected_any["go_slim"].value_counts().items():
        print(f"  {cat:35s}: {n}")

    print("\nGenerating QC figure …")
    fig_qc = make_qc_figure(df_norm, df_raw, df_step1)
    fig_qc.savefig(OUT_QC_PNG, dpi=200, bbox_inches="tight",
                   facecolor="white")
    fig_qc.savefig(OUT_QC_PDF, dpi=200, bbox_inches="tight",
                   facecolor="white")
    plt.close(fig_qc)
    print(f"  Saved → {OUT_QC_PNG} / {OUT_QC_PDF}")

    print("Generating main figure …")
    fig_main = make_main_figure(df_norm, df_raw)
    fig_main.savefig(OUT_MAIN_PNG, dpi=200, bbox_inches="tight",
                     facecolor="white")
    fig_main.savefig(OUT_MAIN_PDF, dpi=200, bbox_inches="tight",
                     facecolor="white")
    plt.close(fig_main)
    print(f"  Saved → {OUT_MAIN_PNG} / {OUT_MAIN_PDF}")

    print("\nDone.")
    print("=" * 60)
