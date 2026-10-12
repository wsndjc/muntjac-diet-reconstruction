# -*- coding: utf-8 -*-
"""Downstream diet indicators for the 65 samples analysed in this study
(54 Reeves' muntjac + 11 black muntjac)."""
import os
import csv, os, math
import numpy as np

RAW = os.environ.get("MUNTJAC_RAW", "raw")  # directory of per-sample blast annotation + TPM files
OUT = os.environ.get("MUNTJAC_OUT", "output")
os.makedirs(OUT, exist_ok=True)

# Sample assignment used in this study: 54 Reeves + 11 black
REEVES = """D11 D14 D15 D18 D19 D23 D3 D34 D38 D46 D48 D50 D51 D5 D52 D53 D57 D61 D6
D62 D63 D66 D67 D68 D69 D70 D71 D7 D72 D86 E20 E21 E23 E24 E32 E33 E35 E36 E37
E42 E46 E56 E57 E65 E75 E76 E81 E82 E83 E84 E85 E89 E91 E92""".split()
BLACK = """D17 D21 D22 D24 D25 D27 D28 D29 D30 D31 D78""".split()

src = open(os.path.join(os.path.dirname(__file__), "build_diet_tables.py"), encoding="utf-8").read()
seg = src[src.index('G2F = {'):src.index('def load_sample')]
ns = {}; exec(seg, ns)
family_of = ns['family_of']

def load_family(sample, thr):
    blast = os.path.join(RAW, f"{sample}_blast_annotation.csv")
    tpmf  = os.path.join(RAW, f"{sample}_tpm_abundance.csv")
    hits = []
    with open(blast, encoding='utf-8-sig', errors='replace') as fh:
        rd = csv.reader(fh); next(rd)
        for row in rd:
            if len(row) < 4: continue
            try: pi = float(row[2]); al = float(row[3])
            except ValueError: continue
            fam = family_of(row[1])
            if fam is not None and pi >= thr and al >= 150:
                hits.append((row[0], fam))
    needed = {q for q, _ in hits}
    tpm = {}; reads = {}
    with open(tpmf, encoding='utf-8-sig', errors='replace') as fh:
        rd = csv.reader(fh); next(rd)
        for row in rd:
            if len(row) < 4: continue
            if row[0] in needed:
                try: tpm[row[0]] = float(row[3]); reads[row[0]] = float(row[2])
                except ValueError: pass
    ft = {}; fr = {}
    for q, fam in hits:
        if q in tpm:
            ft[fam] = ft.get(fam, 0.0) + tpm[q]
            fr[fam] = fr.get(fam, 0.0) + reads[q]
    return ft, fr

def rra_per_sample(sample, thr):
    ft, _ = load_family(sample, thr)
    tot = sum(ft.values())
    if tot <= 0: return {}
    return {f: v / tot * 100.0 for f, v in ft.items()}

def counts_per_sample(sample, thr):
    _, fr = load_family(sample, thr)
    return fr

def total_reads(sample):
    tpmf = os.path.join(RAW, f"{sample}_tpm_abundance.csv")
    tot = 0.0
    with open(tpmf, encoding='utf-8-sig', errors='replace') as fh:
        rd = csv.reader(fh); next(rd)
        for row in rd:
            if len(row) < 3: continue
            try: tot += float(row[2])
            except ValueError: pass
    return tot

def rarefied_family_counts(sample, thr, R, rng):
    _, fr = load_family(sample, thr)
    T = total_reads(sample)
    if T <= 0 or R <= 0: return {f: 0.0 for f in fr}
    fams = sorted(fr)
    other = T - sum(fr.values())
    probs = np.array([other] + [fr[f] for f in fams], dtype=float)
    probs = probs / probs.sum()
    draw = rng.multinomial(int(round(R)), probs)
    return {f: float(draw[i + 1]) for i, f in enumerate(fams)}

def shannon_pielou(rra):
    fams = [f for f, v in rra.items() if v > 0]
    S = len(fams)
    if S == 0: return 0.0, float('nan')
    p = [rra[f] / 100.0 for f in fams]
    H = -sum(x * math.log(x) for x in p)
    if S == 1: return H, float('nan')
    return H, H / math.log(S)

def build(thr):
    reeves = list(REEVES)
    black  = list(BLACK)
    rra = {}; cnt = {}
    for host, samples in (("M. reevesi", reeves), ("M. crinifrons", black)):
        for s in samples:
            rra[(host, s)] = rra_per_sample(s, thr)
            cnt[(host, s)] = counts_per_sample(s, thr)
    ret_r = [(s, rra[("M. reevesi", s)]) for s in reeves if rra[("M. reevesi", s)]]
    ret_b = [(s, rra[("M. crinifrons", s)]) for s in black if rra[("M. crinifrons", s)]]
    fams = sorted({f for _, rr in ret_r + ret_b for f in rr})
    with open(os.path.join(OUT, f"RRA_ASM{thr}_new.csv"), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh); w.writerow(['SampleID', 'Host'] + fams)
        for host, rows in (("M. reevesi", ret_r), ("M. crinifrons", ret_b)):
            for s, rr in rows:
                w.writerow([s, host] + [round(rr.get(f, 0.0), 6) for f in fams])
    with open(os.path.join(OUT, f"Counts_ASM{thr}_new.csv"), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh); w.writerow(['SampleID', 'Host'] + fams)
        for host, rows in (("M. reevesi", ret_r), ("M. crinifrons", ret_b)):
            for s, _rr in rows:
                cc = cnt[(host, s)]; w.writerow([s, host] + [round(cc.get(f, 0.0), 3) for f in fams])
    rng = np.random.default_rng(1)
    retained = [s for s, _ in ret_r] + [s for s, _ in ret_b]
    R = int(round(min(total_reads(s) for s in retained)))
    with open(os.path.join(OUT, f"Counts_rarefied_ASM{thr}_new.csv"), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh); w.writerow(['SampleID', 'Host'] + fams)
        for host, rows in (("M. reevesi", ret_r), ("M. crinifrons", ret_b)):
            for s, _rr in rows:
                rc = rarefied_family_counts(s, thr, R, rng)
                w.writerow([s, host] + [round(rc.get(f, 0.0), 3) for f in fams])
    def diet_table(rows, tag):
        N = len(rows); agg = {}
        for s, rr in rows:
            for f, v in rr.items():
                d = agg.setdefault(f, [0, 0.0]); d[0] += 1; d[1] += v
        out = []
        for f, (occ, ssum) in agg.items():
            foo = occ / N * 100.0; mrra = ssum / N
            out.append([f, occ, round(foo, 6), round(mrra, 6), N])
        out.sort(key=lambda x: -x[3])
        with open(os.path.join(OUT, f"Diet_ASM{thr}_{tag}.csv"), 'w', newline='', encoding='utf-8-sig') as fh:
            w = csv.writer(fh); w.writerow(['Family', 'Occurrence', 'FOO_pct', 'mRRA_pct', 'S_valid'])
            for f, occ, foo, mrra, n in out:
                w.writerow([f, occ, foo, mrra, n])
        return out
    dr = diet_table(ret_r, 'Reeves')
    db = diet_table(ret_b, 'black')
    print(f"ASM{thr} NEW: Reeves {len(ret_r)}/{len(reeves)} retained, black {len(ret_b)}/{len(black)} retained, {len(fams)} families")
    print(f"  Reeves top: " + ", ".join(f"{f}={mrra:.2f}" for f,_,_,mrra,_ in dr[:8]))
    print(f"  black  top: " + ", ".join(f"{f}={mrra:.2f}" for f,_,_,mrra,_ in db[:8]))

if __name__ == '__main__':
    for thr in (95, 90):
        build(thr)
    print("DONE ->", OUT)
