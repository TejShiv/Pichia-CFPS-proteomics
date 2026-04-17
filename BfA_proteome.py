# figure_BfA_proteome.py
#
# Shivakumar et al. — Opening the black box: proteomic and single-vesicle dissection of microsomes in a Pichia pastoris cell-free system
#
# Generates comparative proteomics figures for BfA-treated vs control
# K. phaffii SuperMan5 cell-free lysate fractions (Figure 5).
#
# Outputs
# -------
# Figure_Proteomics_BD_Panel.pdf / .png
#   Panel B  Hierarchical clustering heatmap of row-z-scored log2
#            intensities across all six fractions. Rows are Ward-clustered;
#            column order is fixed as control then BfA.
#   Panel D  Pairwise scatter plot of log10 normalised intensities,
#            BfA microsomes vs control microsomes. Proteins exclusive to
#            one condition are plotted at an imputed minimum intensity.
#            Highlighted categories: ER Resident, Glycosylation,
#            Chaperone/Refolding, Stress Response.
#
# Inputs
# ------
# Proteomics_raw_data.xlsx        — raw LC-MS intensity table
# pichia_to_cerevisiae_orthologs.csv — K. phaffii to S. cerevisiae mapping
# pichia_uniprot_entries.json        — UniProt GO annotations
#
# Usage:  python figure_BfA_proteome.py
#
# Font note: if P052 .otf files are present in the working directory they
# are loaded automatically; otherwise matplotlib falls back to serif.
#
# Dependencies: pandas, numpy, matplotlib, scipy, matplotlib-venn

import pandas as pd
import numpy as np
import re, os, json, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_venn import venn2
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.stats import zscore

# =====================================================================
# FONT AND STYLE SETUP
# =====================================================================
for f in os.listdir('.'):
    if f.endswith('.otf') and 'P052' in f:
        fm.fontManager.addfont(f)
FONT = 'P052'

# Color palette
DEEP_BLUE  = '#1B3A5C'; MED_BLUE = '#2B6CB0'; LIGHT_BLUE = '#4A90D9'
SOFT_BLUE  = '#7FB3DE'; PALE_BLUE = '#B8D4E8'
DARK_GREY  = '#2D3436'; MED_GREY = '#636E72'; LIGHT_GREY = '#B2BEC3'; PALE_GREY = '#DFE6E9'
RUST_RED   = '#C0392B'; TEAL = '#0E8C7F'; GOLD = '#D4A843'; PLUM = '#6C3483'

# Scatter plot styles
HIGHLIGHT_STYLES = {
    'ER Resident':         {'color': RUST_RED,  'marker': 'o', 'size': 180},
    'Glycosylation':       {'color': MED_BLUE,  'marker': 'o', 'size': 180},
    'Chaperone/Refolding': {'color': GOLD,      'marker': 'o', 'size': 180},
    'Stress Response':     {'color': PLUM,      'marker': 'o', 'size': 180},
}

BG_CAT_COLORS = {
    'Ribosome/Translation': '#E8D5B7',
    'Metabolism':           '#C8E6C9',
    'Golgi/Vesicular':      '#F3E5F5',
    'Cytoskeleton':         '#FFECB3',
    'Signaling':            '#B3E5FC',
    'Other':                PALE_GREY,
}

CMAP_RWB = LinearSegmentedColormap.from_list(
    'RdWtBu',
    [(0.0, '#1565C0'), (0.15, '#42A5F5'), (0.35, '#BBDEFB'),
     (0.5, '#FFFFFF'), (0.65, '#FFCDD2'), (0.85, '#EF5350'), (1.0, '#C62828')],
    N=256
)

plt.rcParams.update({
    'font.family': FONT, 'font.size': 11, 'axes.linewidth': 1.0,
    'figure.dpi': 300, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# =====================================================================
# DATA LOADING
# =====================================================================
print("Loading data...")
df = pd.read_excel('Proteomics_raw_data.xlsx', sheet_name='Sheet1')
df.columns = df.columns.str.strip()
df['gene_name'] = df['protein.Description'].apply(
    lambda d: (re.search(r'GN=(\S+)', d) or type('',(),{'group':lambda s,x:None})()).group(1))

bfa_cols = ['BfA_lysate', 'BfA_supernatant', 'BfA_microsomes']
ctrl_cols = ['lysate', 'supernatant', 'microsomes']
df['detected_BfA'] = (df[bfa_cols] > 0).any(axis=1)
df['detected_ctrl'] = (df[ctrl_cols] > 0).any(axis=1)
df['BfA_exclusive'] = df['detected_BfA'] & ~df['detected_ctrl']

# Load ortholog data
ortho_df = pd.read_csv('pichia_to_cerevisiae_orthologs.csv')

# Load and parse UniProt annotations
def parse_uniprot(data):
    info = {'go_terms': [], 'interpro': [], 'subcellular': '', 'gene_name': '', 'keywords': []}
    genes = data.get('genes', [])
    if genes:
        gn = genes[0].get('geneName', {})
        info['gene_name'] = gn.get('value', '') if isinstance(gn, dict) else str(gn)
    for x in data.get('uniProtKBCrossReferences', []):
        db = x.get('database', '')
        props = {p['key']: p['value'] for p in x.get('properties', [])}
        if db == 'GO': info['go_terms'].append(f"{x.get('id','')}: {props.get('GoTerm','')}")
        elif db == 'InterPro': info['interpro'].append(f"{x.get('id','')}: {props.get('EntryName','')}")
    for c in data.get('comments', []):
        if c.get('commentType') == 'SUBCELLULAR LOCATION':
            locs = c.get('subcellularLocations', [])
            info['subcellular'] = '; '.join([l.get('location',{}).get('value','') 
                                              for l in locs if isinstance(l.get('location',{}), dict)])
    for kw in data.get('keywords', []):
        info['keywords'].append(kw.get('name', '') if isinstance(kw, dict) else str(kw))
    return info

with open('pichia_uniprot_entries.json', 'r') as f:
    uniprot_raw = json.load(f)
parsed_uniprot = {acc: parse_uniprot(d) for acc, d in uniprot_raw.items()}

print(f"Loaded: {len(df)} proteins, {len(ortho_df)} orthologs, {len(parsed_uniprot)} UniProt entries")

# =====================================================================
# PREPARE VENN DATA
# =====================================================================
ctrl_lys = set(df[df['lysate'] > 0]['protein.Accession'])
ctrl_mic = set(df[df['microsomes'] > 0]['protein.Accession'])
bfa_lys = set(df[df['BfA_lysate'] > 0]['protein.Accession'])
bfa_mic = set(df[df['BfA_microsomes'] > 0]['protein.Accession'])

# =====================================================================
# PREPARE HEATMAP DATA (with user-specified column order)
# =====================================================================
heatmap_data = df[['protein.Accession'] + bfa_cols + ctrl_cols].copy()
heatmap_data = heatmap_data[(heatmap_data.iloc[:, 1:] > 0).any(axis=1)]
heatmap_vals = np.log2(heatmap_data.iloc[:, 1:].values + 1)
heatmap_zscored = zscore(heatmap_vals, axis=1)
heatmap_zscored[np.isnan(heatmap_zscored).any(axis=1)] = 0

# Column order: Control Lys, Ctrl Sup, Ctrl Mic, BfA Lys, BfA Sup, BfA Mic
# Original: BfA_lys(0), BfA_sup(1), BfA_mic(2), Ctrl_lys(3), Ctrl_sup(4), Ctrl_mic(5)
desired_order = [3, 4, 5, 0, 1, 2]
col_labels = ['Control\nLysate', 'Control\nSupernatant', 'Control\nMicrosomes',
              'BfA\nLysate', 'BfA\nSupernatant', 'BfA\nMicrosomes']

heatmap_df = pd.DataFrame(heatmap_zscored, index=heatmap_data['protein.Accession'].values)
heatmap_reordered = heatmap_df.iloc[:, desired_order].copy()
heatmap_reordered.columns = col_labels

# Row clustering only
Z_rows = linkage(heatmap_reordered.values, method='ward', metric='euclidean')
row_order = leaves_list(Z_rows)
heatmap_ordered = heatmap_reordered.iloc[row_order, :]

print(f"Heatmap: {heatmap_ordered.shape[0]} proteins × {len(col_labels)} conditions")

# =====================================================================
# PREPARE SCATTER DATA (with ortholog + UniProt GO classification)
# =====================================================================
scatter = df[['protein.Accession', 'protein.Description', 'gene_name',
              'BfA_microsomes', 'microsomes']].copy()
scatter = scatter[(scatter['BfA_microsomes'] > 0) | (scatter['microsomes'] > 0)]
scatter = scatter.merge(ortho_df[['pichia_accession', 'cerevisiae_gene_name', 'cerevisiae_function']],
                        left_on='protein.Accession', right_on='pichia_accession', how='left')

min_ctrl = scatter[scatter['microsomes'] > 0]['microsomes'].min()
impute_val = 1000
scatter['log10_BfA'] = np.log10(scatter['BfA_microsomes'].replace(0, impute_val))
scatter['log10_ctrl'] = np.log10(scatter['microsomes'].replace(0, impute_val))

np.random.seed(42)
scatter['log10_BfA_plot'] = scatter['log10_BfA'].copy()
scatter['log10_ctrl_plot'] = scatter['log10_ctrl'].copy()
bfa_only_mask = scatter['microsomes'] == 0
ctrl_only_mask = scatter['BfA_microsomes'] == 0
scatter.loc[bfa_only_mask, 'log10_ctrl_plot'] += np.random.uniform(-0.2, 0.2, bfa_only_mask.sum())
scatter.loc[ctrl_only_mask, 'log10_BfA_plot'] += np.random.uniform(-0.2, 0.2, ctrl_only_mask.sum())

def classify_enhanced(row):
    desc = str(row['protein.Description']).lower()
    gene = str(row['gene_name']).upper() if pd.notna(row['gene_name']) else ''
    acc = row['protein.Accession']
    sc_gene = str(row.get('cerevisiae_gene_name', '')).upper() if pd.notna(row.get('cerevisiae_gene_name')) else ''
    sc_func = str(row.get('cerevisiae_function', '')).lower() if pd.notna(row.get('cerevisiae_function')) else ''
    u = parsed_uniprot.get(acc, {})
    go_text = ' '.join([str(g).lower() for g in u.get('go_terms', [])])
    kw_text = ' '.join([str(k).lower() for k in u.get('keywords', [])])
    subcellular = str(u.get('subcellular', '')).lower()
    combined = f"{desc} {sc_func} {go_text}"

    er_genes = {'PDI1','PDI','SEC61','SEC62','SEC63','SAR1','SSH1','CNE1','KAR2','ERO1','DSL1','SEC39','ERV29','ERV25','GET3','EGD1','EGD2'}
    if gene in er_genes or sc_gene in er_genes: return 'ER Resident'
    if any(k in combined for k in ['protein disulfide-isomerase','dsl1','endoplasmic reticulum','translocon']): return 'ER Resident'
    if 'endoplasmic reticulum' in subcellular: return 'ER Resident'

    glyco_genes = {'ALG1','ALG3','DPM1','OST1','STT3','WBP1','SEC53','MNS1','CWH41','ROT2'}
    if gene in glyco_genes or sc_gene in glyco_genes: return 'Glycosylation'
    if any(k in combined for k in ['glycosyl','mannosyl','glucosidase','dolichol','oligosaccharyl','phosphomannomutase','mannosidase','glycosylated']): return 'Glycosylation'

    chap_genes = {'SSA1','SSA2','SSA3','SSA4','SSB1','SSB2','SSC1','SSE1','SSE2','SSZ1','HSP82','HSC82','HSP104','HSP60','HSP10','CPR1','ZUO1'}
    if gene in chap_genes or sc_gene in chap_genes: return 'Chaperone/Refolding'
    if any(k in combined for k in ['chaperone','chaperonin','heat shock protein','hsp90','hsp70','hsp60','peptidyl-prolyl','protein folding','bag domain','foldase']): return 'Chaperone/Refolding'
    if any(k in go_text for k in ['protein folding','unfolded protein binding','chaperone']): return 'Chaperone/Refolding'

    stress_genes = {'TSA1','TSA2','AHP1','TRX1','TRX2','TRR1','GLR1','SOD1','SOD2','CTT1','STM1'}
    if gene in stress_genes or sc_gene in stress_genes: return 'Stress Response'
    if any(k in combined for k in ['thioredoxin','peroxiredoxin','glutathione','peroxidase','superoxide','catalase']): return 'Stress Response'

    if any(k in combined for k in ['ribosom','translat','eif','elongation factor','initiation factor']): return 'Ribosome/Translation'
    if any(k in go_text for k in ['ribosom','structural constituent of ribosome','translation']): return 'Ribosome/Translation'
    if any(k in kw_text for k in ['ribosomal protein','ribonucleoprotein']): return 'Ribosome/Translation'

    if any(k in combined for k in ['dehydrogenase','synthase','oxidoreductase','reductase','oxidase','enolase','pyruvate','atp synth','nadp','nadh','aldolase','transketolase']): return 'Metabolism'
    if any(k in kw_text for k in ['oxidoreductase','hydrolase','transferase']): return 'Metabolism'

    if any(k in combined for k in ['golgi','vesicl','coat','copii','copi','snare','secretory','clathrin','vacuol']): return 'Golgi/Vesicular'
    if 'golgi' in subcellular: return 'Golgi/Vesicular'
    if any(k in combined for k in ['actin','tubulin','cytoskelet']): return 'Cytoskeleton'
    if any(k in combined for k in ['kinase','phosphatase','calmodulin','gtp']): return 'Signaling'
    return 'Other'

scatter['highlight'] = scatter.apply(classify_enhanced, axis=1)

def resolve_gene(row):
    gene = row['gene_name'] if pd.notna(row['gene_name']) else ''
    sc = row.get('cerevisiae_gene_name', '')
    if isinstance(gene, str) and (gene.startswith('PAS_chr') or gene.startswith('PP7435') or gene == ''):
        if pd.notna(sc) and sc != '': return sc
        desc = str(row['protein.Description']).lower()
        mapping = {'disulfide':'PDI1', 'dsl1':'DSL1', 'tethering':'DSL1', 'glycosylated':'STE13',
                   'hsp60':'HSP60', 'chaperonin':'HSP60', '10 kda heat shock':'HSP10',
                   'thioredoxin reductase':'TRR1', 'peroxiredoxin':'TSA1',
                   'non-chaperonin molecular chaperone':'SSB2',
                   'atpase involved in protein folding and the response':'SSA3',
                   'atpase involved in protein folding and nuclear':'SSA1',
                   'hsp90':'HSP82', 'hsp70 protein that interacts with zuo1':'SSZ1'}
        for key, name in mapping.items():
            if key in desc: return name
        # Hardcoded accession overrides for uncharacterised K. phaffii proteins
        # that only resolve via S. cerevisiae ortholog mapping
        accession_overrides = {
            'C4R0P3': 'SEC39',  # DSL complex subunit (Uncharacterized → S.cer SEC39)
        }
        if row.get('protein.Accession', '') in accession_overrides:
            return accession_overrides[row['protein.Accession']]
        # Hardcoded accession overrides for uncharacterised K. phaffii proteins
        # that only resolve via S. cerevisiae ortholog mapping
        accession_overrides = {
            'C4R0P3': 'SEC39',  # DSL complex subunit (Uncharacterized -> S.cer SEC39)
        }
        if row.get('protein.Accession', '') in accession_overrides:
            return accession_overrides[row['protein.Accession']]
        return None
    return gene

scatter['resolved_gene'] = scatter.apply(resolve_gene, axis=1)
scatter['clean_gene'] = scatter['resolved_gene'].apply(
    lambda x: 'PDI1' if x is not None and x.lower() in ['pdi','pdi1'] else 
              (None if x is not None and (x.startswith('ATY40_') or x.startswith('BA75_')) else x))

# Best per gene for labeling
hl = scatter[scatter['highlight'].isin(HIGHLIGHT_STYLES.keys())].copy()
hl['total'] = hl['BfA_microsomes'] + hl['microsomes']
best_per_gene = hl.dropna(subset=['clean_gene']).sort_values('total', ascending=False).drop_duplicates('clean_gene', keep='first')

# Compute stats
n_shared = ((scatter['BfA_microsomes'] > 0) & (scatter['microsomes'] > 0)).sum()
shared_fc = np.log2(scatter.loc[(scatter['BfA_microsomes'] > 0) & (scatter['microsomes'] > 0), 'BfA_microsomes'] /
                    scatter.loc[(scatter['BfA_microsomes'] > 0) & (scatter['microsomes'] > 0), 'microsomes'])

det_limit = np.log10(min_ctrl) - 0.1
imp_val_log = np.log10(impute_val)
fold2 = np.log10(2)

print(f"Scatter: {len(scatter)} proteins (BfA-only: {bfa_only_mask.sum()}, Ctrl-only: {ctrl_only_mask.sum()}, Shared: {n_shared})")

# =====================================================================
# GENERATE B+D PANEL
# =====================================================================
print("\nGenerating B+D panel...")
fig_bd = plt.figure(figsize=(26, 13), facecolor='white')
gs_bd = gridspec.GridSpec(1, 2, wspace=0.22, left=0.05, right=0.97, top=0.92, bottom=0.10, width_ratios=[1, 1.35])

# Panel B: Heatmap
gs_b = gs_bd[0, 0].subgridspec(1, 2, width_ratios=[0.06, 0.94], wspace=0.003)
ax_dendro = fig_bd.add_subplot(gs_b[0, 0])
dendrogram(Z_rows, orientation='left', no_labels=True, ax=ax_dendro, color_threshold=0, above_threshold_color='#666666')
ax_dendro.set_xticks([]); ax_dendro.set_yticks([])
for spine in ax_dendro.spines.values(): spine.set_visible(False)

ax_heat = fig_bd.add_subplot(gs_b[0, 1])
im = ax_heat.imshow(heatmap_ordered.values, aspect='auto', cmap=CMAP_RWB, interpolation='nearest', vmin=-2, vmax=2)
short_cols = [c.replace('\n', ' ') for c in heatmap_ordered.columns]
ax_heat.set_xticks(range(len(short_cols)))
ax_heat.set_xticklabels(short_cols, rotation=45, ha='right', fontsize=15, fontweight='bold', fontfamily=FONT, color=DARK_GREY)
ax_heat.set_yticks([])
cbar = plt.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, shrink=0.8)
cbar.set_label('Z-score', fontsize=17, fontfamily=FONT, color=DARK_GREY, fontweight='bold')
cbar.set_ticks([-2, -1, 0, 1, 2])
cbar.ax.tick_params(labelsize=14)
ax_heat.set_title('Hierarchical Clustering of Protein Abundances', fontsize=19, fontweight='bold', fontfamily=FONT, color=DARK_GREY, pad=12)
fig_bd.text(0.02, 0.95, 'B', fontsize=34, fontweight='bold', va='top', fontfamily=FONT, color=DARK_GREY)

# Panel D: Scatter
ax_d = fig_bd.add_subplot(gs_bd[0, 1])
ax_d.axvspan(imp_val_log - 0.35, det_limit, color='#FFF8E1', alpha=0.35, zorder=0)
ax_d.axhspan(imp_val_log - 0.35, det_limit, color='#E8EAF6', alpha=0.35, zorder=0)

bg_cats = ['Ribosome/Translation', 'Metabolism', 'Golgi/Vesicular', 'Cytoskeleton', 'Signaling', 'Other']
for cat in bg_cats:
    mask = scatter['highlight'] == cat
    if mask.sum() > 0:
        ax_d.scatter(scatter.loc[mask, 'log10_ctrl_plot'], scatter.loc[mask, 'log10_BfA_plot'],
                    c=BG_CAT_COLORS.get(cat, PALE_GREY), s=28, alpha=0.6, edgecolors='none', zorder=1)

for cat, style in HIGHLIGHT_STYLES.items():
    mask = scatter['highlight'] == cat
    if mask.sum() > 0:
        ax_d.scatter(scatter.loc[mask, 'log10_ctrl_plot'], scatter.loc[mask, 'log10_BfA_plot'],
                    c=style['color'], s=style['size']*1.2, alpha=0.9, edgecolors=DARK_GREY,
                    linewidths=0.7, marker=style['marker'], zorder=3,
                    label=f"{cat} (n = {mask.sum()})")

ax_d.plot([2.5, 7.6], [2.5, 7.6], '--', color=DARK_GREY, linewidth=0.8, alpha=0.3)
ax_d.plot([2.5, 7.6], [2.5+fold2, 7.6+fold2], ':', color=RUST_RED, linewidth=1.2, alpha=0.35)
ax_d.plot([2.5, 7.6], [2.5-fold2, 7.6-fold2], ':', color=MED_BLUE, linewidth=1.2, alpha=0.35)
ax_d.set_xlim(2.5, 7.6); ax_d.set_ylim(2.2, 7.45)

# Manually positioned labels for all 16 highlighted proteins
LABEL_POSITIONS = {
    # BfA-only left cluster — fanned to the right of the cluster
    'SSZ1':  (3.55, 6.40),
    'PDI1':  (3.55, 5.65),
    'TRR1':  (3.55, 5.28),
    'STE13': (3.55, 4.92),
    'SEC39': (3.55, 4.55),
    'DSL1':  (3.55, 4.18),
    # Shared right cluster
    'EGD1':  (4.80, 7.28),
    'STM1':  (5.40, 7.28),
    'HSP60': (5.95, 7.28),
    'SSB2':  (6.50, 7.28),
    'SSC1':  (7.05, 5.35),
    'SSA3':  (7.05, 6.38),
    'SSA1':  (7.05, 6.12),
    'HSP82': (7.05, 5.75),
    # Bottom ctrl-only
    'HSP10': (4.50, 3.55),
    'TSA1':  (6.55, 3.50),
}

for _, r in best_per_gene.iterrows():
    cat = r['highlight']
    if cat not in HIGHLIGHT_STYLES:
        continue
    gene = r['clean_gene']
    if not gene or gene not in LABEL_POSITIONS:
        continue
    tc = HIGHLIGHT_STYLES[cat]['color']
    px, py = r['log10_ctrl_plot'], r['log10_BfA_plot']
    tx, ty = LABEL_POSITIONS[gene]
    ax_d.annotate(gene, xy=(px, py), xytext=(tx, ty),
                  fontsize=28, fontweight='bold', fontstyle='italic',
                  fontfamily=FONT, color=tc,
                  arrowprops=dict(arrowstyle='-', color=MED_GREY, lw=0.8, alpha=0.65),
                  zorder=5,
                  bbox=dict(boxstyle='round,pad=0.12', facecolor='none',
                            alpha=0.0, edgecolor='none'),
                  annotation_clip=False)


ax_d.set_xlabel('Control Microsomes (log$_{10}$ intensity)', fontsize=38, fontweight='bold', fontfamily=FONT, color=DARK_GREY, labelpad=22)
ax_d.set_ylabel('BfA Microsomes (log$_{10}$ intensity)', fontsize=38, fontweight='bold', fontfamily=FONT, color=DARK_GREY, labelpad=22)
# Title removed

leg = ax_d.legend(loc='upper right', fontsize=24, framealpha=0.5, edgecolor=LIGHT_GREY, fancybox=False,
                  title='Protein Category', title_fontsize=24, bbox_to_anchor=(1.0, 0.59),
                  markerscale=1.6, handletextpad=0.8)
for t in leg.get_texts(): t.set_fontfamily(FONT); t.set_color(DARK_GREY)
leg.get_title().set_fontfamily(FONT); leg.get_title().set_color(DARK_GREY)
leg.get_title().set_fontweight('bold'); leg.get_frame().set_alpha(0.5)

ax_d.spines['top'].set_visible(False); ax_d.spines['right'].set_visible(False)
ax_d.spines['left'].set_color(LIGHT_GREY); ax_d.spines['bottom'].set_color(LIGHT_GREY)
ax_d.tick_params(axis='both', colors=DARK_GREY, labelsize=22)
ax_d.grid(alpha=0.06, color=LIGHT_GREY, linewidth=0.3); ax_d.set_axisbelow(True)


fig_bd.text(0.50, 0.95, 'D', fontsize=34, fontweight='bold', va='top', fontfamily=FONT, color=DARK_GREY)

plt.savefig('Figure_Proteomics_BD_Panel.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('Figure_Proteomics_BD_Panel.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: Figure_Proteomics_BD_Panel.png/pdf")

print("\nDone!")
