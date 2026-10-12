# -*- coding: utf-8 -*-
import os
import io, sys, csv, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
RAW = os.environ.get('MUNTJAC_RAW', 'raw')  # directory of per-sample blast annotation + TPM files

src = open(os.path.join(os.path.dirname(__file__), "build_diet_tables.py"), encoding="utf-8").read()
seg = src[src.index('G2F = {'):src.index('def load_sample')]
ns = {}; exec(seg, ns)
family_of = ns['family_of']

REEVES = 'D11 D14 D15 D18 D19 D23 D3 D34 D38 D46 D48 D50 D51 D5 D52 D53 D57 D61 D6 D62 D63 D66 D67 D68 D69 D70 D71 D7 D72 D86 E20 E21 E23 E24 E32 E33 E35 E36 E37 E42 E46 E56 E57 E65 E75 E76 E81 E82 E83 E84 E85 E89 E91 E92'.split()
BLACK = 'D17 D21 D22 D24 D25 D27 D28 D29 D30 D31 D78'.split()

def diet(samples, thr):
    ret = []
    for s in samples:
        bp = os.path.join(RAW, f'{s}_blast_annotation.csv')
        tp = os.path.join(RAW, f'{s}_tpm_abundance.csv')
        tpm = {}
        with open(tp, encoding='utf-8-sig', errors='replace') as fh:
            rd = csv.reader(fh); next(rd)
            for r in rd:
                if len(r) < 4: continue
                try: tpm[r[0]] = float(r[3])
                except ValueError: pass
        ft = {}
        with open(bp, encoding='utf-8-sig', errors='replace') as fh:
            rd = csv.reader(fh); next(rd)
            for r in rd:
                if len(r) < 4: continue
                try: pi = float(r[2]); al = float(r[3])
                except ValueError: continue
                f = family_of(r[1])
                if f and pi >= thr and al >= 150 and r[0] in tpm and tpm[r[0]] > 0:
                    ft[f] = ft.get(f, 0.0) + tpm[r[0]]
        tot = sum(ft.values())
        if tot > 0:
            ret.append({f: v / tot * 100.0 for f, v in ft.items()})
    N = len(ret); agg = {}
    for rr in ret:
        for f, v in rr.items():
            agg[f] = agg.get(f, 0.0) + v
    out = [(f, round(agg[f] / N, 2)) for f in agg]
    out.sort(key=lambda x: -x[1])
    return out, N

for sp, samples in [('Reeves', REEVES), ('black', BLACK)]:
    print(f'===== {sp} =====')
    tables = {}
    for thr in (90, 95, 99):
        out, N = diet(samples, thr)
        tables[thr] = {f: v for f, v in out}
        print(f'  ASM{thr}: {len(out)} families, {N} retained samples')
        print('     top8: ' + ', '.join(f'{f}={v:.1f}' for f, v in out[:8]))
    # families lost going 90->95 and 95->99, with their ASM90 mRRA (to show they are minor)
    for lo, hi in [(90, 95), (95, 99)]:
        lost = set(tables[lo]) - set(tables[hi])
        lost_desc = [(f, tables[lo][f]) for f in sorted(lost, key=lambda f: -tables[lo][f])]
        print(f'  {lo}->{hi}: lost {len(lost)} families; their mRRA at ASM{lo}: ' + ', '.join(f'{f}={v:.2f}' for f, v in lost_desc[:12]))
