# -*- coding: utf-8 -*-
"""建立 TPM 路径的样本×科 RRA 矩阵并缓存 (复现论文口径)
   过滤条件: length>=150, evalue<1e-5(or 1e-8), pident>=thr, TPM>0
"""
import os
import io, os, re, sys, glob, pickle
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd, numpy as np

DATA = os.environ.get("MUNTJAC_DATA", ".")
# 中间产物（per-sample blast_annotation.csv + tpm_abundance.csv，约 1.6 GB）
# 不随数据包分发；如需要可向作者索取，或由 02_ASM_assembly_blast_tpm.sbatch 重新生成。
BL = os.environ.get("MUNTJAC_INTERMEDIATE", os.path.join(DATA, "intermediate"))
DICT = os.path.join(DATA, "06_species_family_dictionary", "species_family_dictionary.csv")
CACHE = os.path.join(DATA, "_cache", "tpm_cache.pkl")

HOST = {}
for l in io.open(os.path.join(DATA, "05_sample_metadata", "sample_species_map.txt"), encoding="utf-8"):
    p = l.split()
    if len(p) >= 2:
        HOST[p[0]] = "black" if p[1] == "黑麂" else "reevesi"

d = pd.read_csv(DICT, encoding="gbk")
s2f = {str(r.Plant_lower).strip().lower(): str(r.Family).strip() for r in d.itertuples()}
g2f = {}
for k, v in s2f.items():
    g2f.setdefault(k.split()[0], v)
FIX = {"ltea": "Itea", "lndocalamus": "Indocalamus", "llex": "Ilex", "ctinidia": "Actinidia"}
# R 代码硬编码：这三个属字典里没有，直接指定科
FAM_FIX = {"itea": "Iteaceae", "indocalamus": "Poaceae", "ilex": "Aquifoliaceae"}
SEQ_FIX = {"seq_695": ("Kadsura longipedunculata", "Schisandraceae")}

def parse(sseqid):
    s = re.sub(r"_{2,}", "_", str(sseqid))
    s = re.sub(r"^_*seq_\d+_*", "", s)
    parts = [p for p in s.split("_") if p != ""]
    if not parts:
        return None, None
    g = FIX.get(parts[0].lower(), parts[0].capitalize())
    sp = "_".join(parts[1:]) if len(parts) > 1 else ""
    return (g + " " + sp).strip().lower(), g.lower()

def build(pident_min, evalue_max=1e-5):
    frames = []
    for f in sorted(glob.glob(os.path.join(BL, "*_blast_annotation.csv"))):
        s = os.path.basename(f).split("_")[0]
        if s not in HOST:
            continue
        b = pd.read_csv(f, usecols=["qseqid", "sseqid", "pident", "length", "evalue"])
        b = b[(b.length >= 150) & (b.evalue < evalue_max) & (b.pident >= pident_min)]
        if b.empty:
            continue
        b["Sample"] = s
        frames.append(b)
    blast = pd.concat(frames, ignore_index=True)
    tf = []
    for s, grp in blast.groupby("Sample"):
        t = pd.read_csv(os.path.join(BL, f"{s}_tpm_abundance.csv"))
        t = t[t.Contig.isin(set(grp.qseqid)) & (t.TPM > 0)].copy()
        t["Sample"] = s
        tf.append(t[["Contig", "TPM", "Sample"]])
    tpm = pd.concat(tf, ignore_index=True)
    df = blast.merge(tpm, left_on=["qseqid", "Sample"], right_on=["Contig", "Sample"])
    df["Pl"] = [parse(x)[0] for x in df.sseqid]
    df["Gl"] = [parse(x)[1] for x in df.sseqid]
    df["Family"] = [s2f.get(p) if p else None for p in df.Pl]
    m = df.Family.isna()
    df.loc[m, "Family"] = [g2f.get(g) if g else None for g in df.loc[m, "Gl"]]
    # R 代码的硬编码修正
    df.loc[df.Gl.isin(FAM_FIX), "Family"] = df.loc[df.Gl.isin(FAM_FIX), "Gl"].map(FAM_FIX)
    for sid, (nm, fam) in SEQ_FIX.items():
        df.loc[df.sseqid == sid, "Family"] = fam
    df = df[df.Family.notna()].copy()
    df["Host"] = df.Sample.map(HOST)
    agg = df.groupby(["Sample", "Family"], as_index=False).TPM.sum()
    mat = agg.pivot(index="Sample", columns="Family", values="TPM").fillna(0.0)
    meta = df.groupby("Sample").Host.first()
    return mat, meta

if __name__ == "__main__":
    store = {}
    for thr in [90, 95, 99]:
        for ev in [1e-5, 1e-8]:
            store[(thr, ev)] = build(thr, ev)
            mat, meta = store[(thr, ev)]
            print(f"pident>={thr} evalue<{ev}: 样本 {mat.shape[0]} × 科 {mat.shape[1]}")
    pickle.dump(store, open(CACHE, "wb"))
    print("cached ->", CACHE)
