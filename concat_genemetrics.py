#!/usr/bin/env python3

import sys
import glob
import os
import pandas as pd
import numpy as np

def infer_cn(log2):
    if pd.isna(log2):
        return np.nan, np.nan
    if log2 <= -1.1:
        return 0, "HOMDEL"
    elif log2 <= -0.4:
        return 1, "DEL"
    elif log2 < 0.3:
        return 2, "NEUTRAL"
    elif log2 < 0.7:
        return 3, "GAIN"
    else:
        return 4, "AMP"

def main(genemetrics_dir, prefix):
    files = sorted(glob.glob(os.path.join(genemetrics_dir, "*.genemetrics.tsv")))
    if not files:
        sys.exit("ERROR: No *.genemetrics.tsv files found")

    log2_df = None
    cn_df = None
    infer_df = None

    for f in files:
        sample = os.path.basename(f).replace(".genemetrics.tsv", "")

        df = pd.read_csv(f, sep="\t")
        df = df[["gene", "chromosome", "start", "end", "log2"]]

        cn_vals = df["log2"].apply(infer_cn)
        df["cn"] = [x[0] for x in cn_vals]
        df["infer"] = [x[1] for x in cn_vals]

        base = df[["gene", "chromosome", "start", "end"]]

        if log2_df is None:
            log2_df = base.copy()
            cn_df = base.copy()
            infer_df = base.copy()

        log2_df[sample] = df["log2"]
        cn_df[sample] = df["cn"]
        infer_df[sample] = df["infer"]

    sort_cols = ["chromosome", "start", "end"]
    log2_df = log2_df.sort_values(sort_cols)
    cn_df = cn_df.sort_values(sort_cols)
    infer_df = infer_df.sort_values(sort_cols)

    log2_df.to_csv(f"{prefix}_log2.tsv", sep="\t", index=False)
    cn_df.to_csv(f"{prefix}_cn.tsv", sep="\t", index=False)
    infer_df.to_csv(f"{prefix}_infer.tsv", sep="\t", index=False)

    print("Written:")
    print(f"  {prefix}_log2.tsv")
    print(f"  {prefix}_cn.tsv")
    print(f"  {prefix}_infer.tsv")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <genemetrics_dir> <output_prefix>")
    main(sys.argv[1], sys.argv[2])
