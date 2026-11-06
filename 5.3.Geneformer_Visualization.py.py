# conda activate Geneformer
# git lfs install
# git clone https://huggingface.co/ctheodoris/Geneformer
# cd Geneformer
# pip install .

import os
import glob
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon
from scipy import stats
from statsmodels.stats.multitest import multipletests
from itertools import combinations

BASE_DIR = "/WData/WorkSpace/xiangjun2024/LiJunHu/Data/Geneformer"
DICT_DIR = "/SoftWare/Geneformer"

PERTURB_DIR = os.path.join(BASE_DIR, "geneformer_perturb_output_single_UP")
GENE_DIR = os.path.join(BASE_DIR, "gene_lists_from_UP")
SIG_DIR = os.path.join(BASE_DIR, "geneformer_significant_affected_genes")

# 图片输出目录改到这里
PLOT_DIR = "/WData/WorkSpace/xiangjun2024/LiJunHu/Figure/8.Gene.Knockdown"

for d in [SIG_DIR, PLOT_DIR]:
    os.makedirs(d, exist_ok=True)

groups = ["cluster1", "cluster2", "cluster3", "cluster4", "cluster5"]

colors = {
    "cluster1": "#C3D0E3",
    "cluster2": "#A4E1E8",
    "cluster3": "#79D5C3",
    "cluster4": "#F9CBAA",
    "cluster5": "#F2B5B5",
}

# 原始字典：通常是 ensembl -> token
gene_id_token_dict = pickle.load(
    open(os.path.join(DICT_DIR, "geneformer/token_dictionary_gc104M.pkl"), "rb")
)

# 反转：token -> ensembl
token_to_gene_id = {v: k for k, v in gene_id_token_dict.items()}

# symbol -> ensembl
gene_name_id_dict = pickle.load(
    open(os.path.join(DICT_DIR, "geneformer/gene_name_id_dict_gc104M.pkl"), "rb")
)

ensembl_to_symbol = {v: k for k, v in gene_name_id_dict.items() if v is not None}

print("Step 1: extracting significant affected genes ...")

summary_effects = {}

for group in groups:
    usable_file = os.path.join(GENE_DIR, f"{group}_usable_UP_genes.csv")
    if not os.path.exists(usable_file):
        print(f"skip {group}: usable gene file not found")
        continue
    usable_df = pd.read_csv(usable_file)
    group_targets = usable_df["ensembl_id"].dropna().unique().tolist()
    subtype_frames = []
    subtype_effect_values = []
    for gid in group_targets:
        in_dir = os.path.join(PERTURB_DIR, group, gid)
        pkl_files = glob.glob(os.path.join(in_dir, "*.pickle"))
        if len(pkl_files) == 0:
            continue
        pair_to_values = {}
        for pf in pkl_files:
            with open(pf, "rb") as fh:
                data = pickle.load(fh)
            for k, v in data.items():
                # gene-related输出：tuple key，且排除 cell_emb
                if isinstance(k, tuple) and len(k) >= 2 and "cell_emb" not in str(k[1]):
                    pair_to_values.setdefault(k, []).extend(v)
                # 保留整体 effect 用于绘图
                if "cell_emb" in str(k):
                    vals = np.array(v, dtype=float)
                    vals = vals[~np.isnan(vals)]
                    subtype_effect_values.extend(1 - vals)
        rows = []
        for k, vals in pair_to_values.items():
            vals = np.array(vals, dtype=float)
            vals = vals[~np.isnan(vals)]
            if vals.size < 10:
                continue
            mean_cos = vals.mean()
            effect = 1 - mean_cos
            try:
                stat, p = wilcoxon(vals - 1.0, alternative="two-sided")
            except Exception:
                p = np.nan
            pert_token = k[0]
            aff_token = k[1]
            pert_ensembl = token_to_gene_id.get(pert_token, np.nan)
            aff_ensembl = token_to_gene_id.get(aff_token, np.nan)
            rows.append({
                "CellSubtype": group,
                "Perturbed_Ensembl_ID": pert_ensembl,
                "Perturbed_gene": ensembl_to_symbol.get(pert_ensembl, pert_ensembl),
                "Affected_Ensembl_ID": aff_ensembl,
                "Affected_gene": ensembl_to_symbol.get(aff_ensembl, aff_ensembl),
                "N": vals.size,
                "mean_cosine_sim": mean_cos,
                "effect": effect,
                "pvalue": p
            })
        if len(rows) == 0:
            continue
        df = pd.DataFrame(rows)
        df["FDR"] = multipletests(df["pvalue"].fillna(1), method="fdr_bh")[1]
        df = df.sort_values(["FDR", "effect"], ascending=[True, False])
        df.to_csv(
            os.path.join(SIG_DIR, f"{group}_{gid}_affected_genes_full.csv"),
            index=False
        )
        sig = df.loc[
            (df["FDR"] < 0.05) &
            (df["effect"] > 0.01) &
            (df["N"] >= 10),
            :
        ].copy()
        sig.to_csv(
            os.path.join(SIG_DIR, f"{group}_{gid}_affected_genes_significant.csv"),
            index=False
        )
        subtype_frames.append(sig)
        print(f"{group} | {gid} | significant affected genes: {sig.shape[0]}")
    if len(subtype_frames) > 0:
        all_sig = pd.concat(subtype_frames, axis=0, ignore_index=True)
        all_sig.to_csv(
            os.path.join(SIG_DIR, f"{group}_all_significant_affected_genes.csv"),
            index=False
        )
    else:
        pd.DataFrame().to_csv(
            os.path.join(SIG_DIR, f"{group}_all_significant_affected_genes.csv"),
            index=False
        )
    summary_effects[group] = np.array(subtype_effect_values, dtype=float)

print("Significant affected genes exported.")

print("Step 2: generating per-subtype top affected gene tables ...")

for group in groups:
    all_sig_file = os.path.join(SIG_DIR, f"{group}_all_significant_affected_genes.csv")
    if not os.path.exists(all_sig_file):
        continue
    df = pd.read_csv(all_sig_file)
    if df.shape[0] == 0:
        print(f"{group}: no significant affected genes")
        continue
    agg = (
        df.groupby(["Affected_Ensembl_ID", "Affected_gene"], dropna=False)
        .agg(
            mean_effect=("effect", "mean"),
            min_FDR=("FDR", "min"),
            n_perturbations=("Perturbed_gene", "nunique")
        )
        .reset_index()
        .sort_values(["mean_effect", "min_FDR"], ascending=[False, True])
    )
    agg.to_csv(os.path.join(SIG_DIR, f"{group}_top_affected_genes_summary.csv"), index=False)
    print(f"{group}: exported top affected gene summary")

print("Step 3: plotting ...")

# 直方图保存为 PDF
fig, axes = plt.subplots(2, 3, figsize=(12, 10))
axes = axes.flatten()

for ax, group in zip(axes, groups):
    effect = summary_effects.get(group, np.array([]))
    if effect.size == 0:
        ax.set_title(f"{group} - no data")
        ax.axis("off")
        continue
    ax.hist(
        effect,
        bins=40,
        color=colors[group],
        edgecolor="black",
        alpha=0.8
    )
    ax.axvline(
        effect.mean(),
        color="black",
        linestyle="--",
        label=f"mean = {effect.mean():.5f}"
    )
    ax.set_title(f"{group} - perturbation effect")
    ax.set_xlabel("Perturbation effect (1 - cosine similarity)")
    ax.set_ylabel("Count")
    ax.legend()

plt.tight_layout()
plt.savefig(
    os.path.join(PLOT_DIR, "1.1.UP_gene_perturbation_hist_MP1_MP4.pdf"),
    format="pdf",
    bbox_inches="tight"
)
plt.close()

# 箱线图保存为 PDF
plot_data = []
plot_labels = []

for group in groups:
    effect = summary_effects.get(group, np.array([]))
    if effect.size > 0:
        plot_data.append(effect)
        plot_labels.append(group)

fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True)

for patch, label in zip(bp["boxes"], plot_labels):
    patch.set_facecolor(colors[label])
    patch.set_alpha(0.8)

ax.set_ylabel("Perturbation effect (1 - cosine similarity)")
ax.set_xlabel("CellSubtype")
ax.set_title("Perturbation effect across MP1-MP4")

plt.tight_layout()
plt.savefig(
    os.path.join(PLOT_DIR, "1.2.UP_gene_perturbation_boxplot_MP1_MP4.pdf"),
    format="pdf",
    bbox_inches="tight"
)
plt.close()

print("PDF figures saved to:", PLOT_DIR)

print("Step 4: statistical tests ...")

valid_effects = []
valid_groups = []

for group in groups:
    effect = summary_effects.get(group, np.array([]))
    if effect.size > 0:
        valid_effects.append(effect)
        valid_groups.append(group)

if len(valid_effects) >= 2:
    f_stat, p_val = stats.f_oneway(*valid_effects)
    print(f"\nANOVA: F={f_stat:.4f}, p={p_val:.4e}")
    with open(os.path.join(SIG_DIR, "anova_result.txt"), "w") as f:
        f.write(f"ANOVA: F={f_stat:.4f}, p={p_val:.4e}\n")
else:
    print("Not enough groups for ANOVA.")

pairwise_rows = []
for g1, g2 in combinations(groups, 2):
    e1 = summary_effects.get(g1, np.array([]))
    e2 = summary_effects.get(g2, np.array([]))
    if len(e1) == 0 or len(e2) == 0:
        print(f"{g1} vs {g2}: skip")
        continue
    t_stat, p_val = stats.ttest_ind(e1, e2, equal_var=False)
    pairwise_rows.append({
        "group1": g1,
        "group2": g2,
        "t_stat": t_stat,
        "pvalue": p_val
    })
    print(f"{g1} vs {g2}: t={t_stat:.4f}, p={p_val:.4e}")

pd.DataFrame(pairwise_rows).to_csv(
    os.path.join(SIG_DIR, "pairwise_ttest_results.csv"),
    index=False
)

print("All analysis finished.")
