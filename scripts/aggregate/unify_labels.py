#!/usr/bin/env python3
# Unify the aggregate typing tiers into one per-cell label table for all downstream analysis.
# Left-joins tier-1 coarse (compartment) + tier-2 immune-rescued + tier-2 stroma-rescued onto
# the full 3.27M-cell cohort, then resolves a single `cell_subtype` keyed by compartment:
#   immune -> immune_subtype | stroma -> stroma_subtype | tumor -> "Tumor" | else "unassigned".
# Tumor carries no finer subtype by design (single 4T1 clone; de-novo subtyping rejected as
# panel noise -- see plan-aggregate.md). Asserts every immune/stroma compartment cell has a
# tier-2 label (the tier-2 tables are exact subsets of their compartment by construction).
# Args: <coarse.parquet> <immune_rescued.parquet> <stroma_rescued.parquet> <out.parquet> <summary_tsv>
import sys
import pandas as pd

coarse_pq, immune_pq, stroma_pq, out_pq, summ_tsv = sys.argv[1:6]

base = pd.read_parquet(coarse_pq)                                  # cell, leiden, coarse_label, compartment, dataset, slide_id
imm  = pd.read_parquet(immune_pq, columns=["cell", "immune_subtype", "rescued"]) \
         .rename(columns={"immune_subtype": "_imm", "rescued": "_imm_resc"})
strm = pd.read_parquet(stroma_pq, columns=["cell", "stroma_subtype", "rescued"]) \
         .rename(columns={"stroma_subtype": "_strm", "rescued": "_strm_resc"})

n0 = len(base)
df = base.merge(imm, on="cell", how="left").merge(strm, on="cell", how="left")
assert len(df) == n0, (len(df), n0)                               # left join must not duplicate

# coverage: tier-2 tables are exact subsets of their compartment
imm_cells  = df["compartment"].eq("immune")
strm_cells = df["compartment"].eq("stroma")
assert df.loc[imm_cells,  "_imm"].notna().all(),  "immune cell(s) without a tier-2 immune label"
assert df.loc[strm_cells, "_strm"].notna().all(), "stroma cell(s) without a tier-2 stroma label"

# resolve single fine label keyed by compartment
df["cell_subtype"] = "unassigned"
df.loc[df["compartment"].eq("tumor"), "cell_subtype"] = "Tumor"
df.loc[imm_cells,  "cell_subtype"] = df.loc[imm_cells,  "_imm"].to_numpy()
df.loc[strm_cells, "cell_subtype"] = df.loc[strm_cells, "_strm"].to_numpy()

src = pd.Series("coarse_unassigned", index=df.index)
src[df["compartment"].eq("tumor")] = "tumor_identity"
src[imm_cells]  = "immune_tier2"
src[strm_cells] = "stroma_tier2"
df["subtype_source"] = src.to_numpy()

df["rescued"] = (df["_imm_resc"] == True) | (df["_strm_resc"] == True)   # NaN (tumor/unassigned) -> False

assert df["cell_subtype"].notna().all() and df["compartment"].notna().all()

out = df[["cell", "dataset", "slide_id", "leiden", "coarse_label",
          "compartment", "cell_subtype", "subtype_source", "rescued"]]
out.to_parquet(out_pq, index=False)

summ = (out.groupby(["compartment", "cell_subtype"]).size()
        .rename("n_cells").reset_index().sort_values("n_cells", ascending=False))
summ["frac"] = (summ["n_cells"] / len(out)).round(4)
summ.to_csv(summ_tsv, sep="\t", index=False)

print(f"unified labels: {len(out)} cells")
print("\n--- compartment ---")
print(out["compartment"].value_counts(dropna=False).to_string())
print("\n--- cell_subtype ---")
print(out["cell_subtype"].value_counts(dropna=False).to_string())
print(f"\nrescued cells: {int(out['rescued'].sum())}")
print(f"\nwrote {out_pq}\nwrote {summ_tsv}")
