
import os
import shutil
import glob
import pickle
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

from geneformer import TranscriptomeTokenizer, InSilicoPerturber

# =========================
# 0. 
# =========================
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
os.environ["PYTORCH_ALLOC_CONF"] = "expandable_segments:True"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
warnings.filterwarnings("ignore")

BASE_DIR = "/WData/WorkSpace/xiangjun2024/LiJunHu/Data/Geneformer"
DICT_DIR = "/SoftWare/Geneformer"
MODEL_DIR = "/SoftWare/Geneformer/Geneformer-V2-104M"

MTX_FILE = os.path.join(BASE_DIR, "malignant_counts.mtx")
GENE_FILE = os.path.join(BASE_DIR, "malignant_genes.csv")
META_FILE = os.path.join(BASE_DIR, "malignant_meta.csv")
H5AD_FILE = os.path.join(BASE_DIR, "malignant.h5ad")

UP_DIR = os.path.join(BASE_DIR, "up_genes")
CLEAN_DIR = os.path.join(BASE_DIR, "clean_h5ad_input")
PROC_DIR = os.path.join(BASE_DIR, "geneformer_processed")
PERTURB_DIR = os.path.join(BASE_DIR, "geneformer_perturb_output_single_UP")
GENE_DIR = os.path.join(BASE_DIR, "gene_lists_from_UP")

for d in [CLEAN_DIR, PROC_DIR, PERTURB_DIR, GENE_DIR]:
    os.makedirs(d, exist_ok=True)

groups = ["cluster1", "cluster2", "cluster3", "cluster4", "cluster5"]

# =========================
# 1. mtx + genes + meta -> h5ad
# =========================
print("Step 1: building h5ad ...")

X = mmread(MTX_FILE).tocsr()
genes = pd.read_csv(GENE_FILE)
meta = pd.read_csv(META_FILE, index_col=0)

gene_names = genes.iloc[:, 0].astype(str).values

print("Matrix shape:", X.shape)
print("Gene number:", len(gene_names))
print("Meta rows:", meta.shape[0])

assert X.shape[0] == len(gene_names), "gene number mismatch"
assert X.shape[1] == meta.shape[0], "cell number mismatch"

adata = sc.AnnData(X=X.T,obs=meta.copy(),var=pd.DataFrame(index=gene_names))

adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.write(H5AD_FILE)
print("Saved:", H5AD_FILE)

# =========================
# 2. gene symbol -> ensembl id
# =========================
print("Step 2: mapping gene symbols to Ensembl IDs ...")

adata = sc.read_h5ad(H5AD_FILE)

if "n_counts" not in adata.obs.columns:
    adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()

gene_name_id = pickle.load(
    open(os.path.join(DICT_DIR, "geneformer/gene_name_id_dict_gc104M.pkl"), "rb")
)

token_dict = pickle.load(
    open(os.path.join(DICT_DIR, "geneformer/token_dictionary_gc104M.pkl"), "rb")
)

adata.var["ensembl_id"] = [gene_name_id.get(x, "") for x in adata.var_names]
adata = adata[:, adata.var["ensembl_id"] != ""].copy()
adata.var_names = list(adata.var["ensembl_id"])
adata.var_names_make_unique()

print("Converted gene number:", adata.n_vars)
print("CellCluster counts:")
print(adata.obs["CellCluster"].value_counts())

# =========================
# 3. 按 cluster1-cluster5 分组保存 h5ad
# 并导出每个 MP 的全部基因
# =========================
print("Step 3: splitting by CellCluster ...")

for group in groups:
    sub = adata[adata.obs["CellCluster"] == group].copy()
    out = os.path.join(CLEAN_DIR, f"Malignant_{group}.h5ad")
    sub.write_h5ad(out)
    pd.DataFrame({"ensembl_id": sub.var_names}).to_csv(
        os.path.join(GENE_DIR, f"{group}_all_genes.csv"),
        index=False
    )
    print(f"{group}: cells={sub.n_obs}, genes={sub.n_vars}, saved={out}")

# =========================
# 4. tokenize 每个 MP
# =========================
print("Step 4: tokenizing each MP ...")

for group in groups:
    src_h5ad = os.path.join(CLEAN_DIR, f"Malignant_{group}.h5ad")
    group_dir = os.path.join(CLEAN_DIR, group)
    os.makedirs(group_dir, exist_ok=True)
    dst_h5ad = os.path.join(group_dir, f"Malignant_{group}.h5ad")
    if not os.path.exists(dst_h5ad):
        shutil.copy(src_h5ad, dst_h5ad)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"CellCluster": "CellCluster"},
        nproc=1,
        model_input_size=2048,
        special_token=True,
        gene_median_file=os.path.join(
            DICT_DIR, "geneformer/gene_median_dictionary_gc104M.pkl"
        ),
        token_dictionary_file=os.path.join(
            DICT_DIR, "geneformer/token_dictionary_gc104M.pkl"
        )
    )
    tk.tokenize_data(
        data_directory=group_dir,
        output_directory=os.path.join(PROC_DIR, group),
        output_prefix=f"Malignant_{group}",
        file_format="h5ad"
    )
    print(f"Tokenized: {group}")

# =========================
# 5. 读取每个 MP 的 UP genes
# 筛选真正可用于 Geneformer 的候选扰动基因
# 条件：
#   1) 能映射 Ensembl
#   2) 在 token dict
#   3) 在该 MP 数据中存在
# =========================
print("Step 5: preparing usable perturbation genes from UP genes ...")

usable_targets = {}

for group in groups:
    up_file = os.path.join(UP_DIR, f"{group}_UP_genes.csv")
    adata_file = os.path.join(CLEAN_DIR, f"Malignant_{group}.h5ad")
    up_df = pd.read_csv(up_file)
    adata_sub = sc.read_h5ad(adata_file)
    rows = []
    for gene in up_df["Symbol"].astype(str).tolist():
        ensembl_id = gene_name_id.get(gene, None)
        in_token = ensembl_id in token_dict if ensembl_id else False
        in_data = ensembl_id in adata_sub.var_names if ensembl_id else False
        rows.append({
            "gene_symbol": gene,
            "ensembl_id": ensembl_id,
            "in_token_dict": in_token,
            "in_subtype_data": in_data,
            "usable_for_geneformer": in_token and in_data
        })
    rows = pd.DataFrame(rows)
    rows.to_csv(os.path.join(GENE_DIR, f"{group}_UP_genes_mapped.csv"), index=False)
    usable = rows.loc[rows["usable_for_geneformer"], :].copy()
    usable.to_csv(os.path.join(GENE_DIR, f"{group}_usable_UP_genes.csv"), index=False)
    usable_targets[group] = usable["ensembl_id"].dropna().unique().tolist()
    usable[["gene_symbol"]].drop_duplicates().to_csv(
        os.path.join(GENE_DIR, f"{group}_usable_UP_gene_symbols_only.csv"),
        index=False
    )
    print(f"{group}: total UP genes={rows.shape[0]}, usable={len(usable_targets[group])}")

# =========================
# 6. 
# =========================
print("Step 6: running single-gene perturbation ...")
for group in groups:
    input_file = os.path.join(PROC_DIR, group, f"Malignant_{group}.dataset")
    if not os.path.exists(input_file):
        print(f"skip {group}: dataset not found")
        continue
    group_targets = usable_targets[group]
    if len(group_targets) == 0:
        print(f"skip {group}: no usable UP genes")
        continue
    for gid in group_targets:
        out_dir = os.path.join(PERTURB_DIR, group, gid)
        os.makedirs(out_dir, exist_ok=True)
        old_pkls = glob.glob(os.path.join(out_dir, "*.pickle"))
        if len(old_pkls) > 0:
            print(f"skip existing: {group} | {gid}")
            continue
        print(f"Running {group} | {gid}")
        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[gid],
            combos=0,
            model_type="Pretrained",
            emb_mode="cls_and_gene",
            max_ncells=50,
            token_dictionary_file=os.path.join(
                DICT_DIR, "geneformer/token_dictionary_gc104M.pkl"
            )
        )
        isp.perturb_data(
            model_directory=MODEL_DIR,
            input_data_file=input_file,
            output_directory=out_dir,
            output_prefix=f"{group}_{gid}_KO_genelevel"
        )    

print("All perturbation jobs finished.")
