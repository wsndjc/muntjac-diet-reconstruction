# -*- coding: utf-8 -*-
"""Reconstruct family-level ASM diet tables (mRRA, FOO). Single-pass, fast, writes Diet_*.csv."""
import os
import csv, os, sys, io, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

RAW = os.environ.get("MUNTJAC_RAW", "raw")  # directory of per-sample blast annotation + TPM files
OUT = "results"
os.makedirs(OUT, exist_ok=True)

REEVES = """D11 D14 D15 D18 D19 D23 D3 D34 D38 D46 D48 D50 D51 D5 D52 D53 D57 D61 D6
D62 D63 D66 D67 D68 D69 D70 D71 D7 D72 D86 E20 E21 E23 E32 E33 E35 E36 E37
E42 E46 E56 E57 E65 E75 E76 E82 E83 E84 E85 E89 E91 E92""".split()
BLACK = """D17 D21 D22 D24 D25 D27 D28 D29 D30 D31 D78 E24 E81""".split()

G2F = {
 'acer':'Sapindaceae','adinandra':'Pentaphylacaceae','akebia':'Lardizabalaceae',
 'alangium':'Cornaceae','alniphyllum':'Styracaceae','alyxia':'Apocynaceae',
 'amelanchier':'Rosaceae','ampelopsis':'Vitaceae','anoectochilus':'Orchidaceae',
 'aphananthe':'Cannabaceae','arachniodes':'Dryopteridaceae','ardisia':'Primulaceae',
 'arisaema':'Araceae','aster':'Asteraceae','betula':'Betulaceae',
 'broussonetia':'Moraceae','callerya':'Fabaceae','callicarpa':'Lamiaceae',
 'camellia':'Theaceae','carpinus':'Betulaceae','castanopsis':'Fagaceae',
 'celtis':'Cannabaceae','chengiodendron':'Oleaceae','choerospondias':'Anacardiaceae',
 'cladrastis':'Fabaceae','clematis':'Ranunculaceae','cleyera':'Pentaphylacaceae',
 'commelina':'Commelinaceae','cornus':'Cornaceae','corylopsis':'Hamamelidaceae',
 'daphne':'Thymelaeaceae','daphniphyllum':'Daphniphyllaceae','dendropanax':'Araliaceae',
 'dioscorea':'Dioscoreaceae','diospyros':'Ebenaceae','diplopterygium':'Gleicheniaceae',
 'dryopteris':'Dryopteridaceae','elaeocarpus':'Elaeocarpaceae','elatostema':'Urticaceae',
 'eriobotrya':'Rosaceae','erythroxylum':'Erythroxylaceae','euonymus':'Celastraceae',
 'eurya':'Pentaphylacaceae','fagus':'Fagaceae','ficus':'Moraceae',
 'frangula':'Rhamnaceae','gardenia':'Rubiaceae','holboellia':'Lardizabalaceae',
 'houttuynia':'Saururaceae','huperzia':'Lycopodiaceae','hydrangea':'Hydrangeaceae',
 'idesia':'Salicaceae','illicium':'Schisandraceae','kadsura':'Schisandraceae',
 'lindera':'Lauraceae','liparis':'Orchidaceae','liquidambar':'Altingiaceae',
 'liriodendron':'Magnoliaceae','lithocarpus':'Fagaceae','llex':'Aquifoliaceae',
 'lndocalamus':'Poaceae','lonicera':'Caprifoliaceae','lophatherum':'Poaceae',
 'loropetalum':'Hamamelidaceae','ltea':'Iteaceae','lysimachia':'Primulaceae',
 'maclura':'Moraceae','magnolia':'Magnoliaceae','mahonia':'Berberidaceae',
 'mallotus':'Euphorbiaceae','malus':'Rosaceae','meliosma':'Sabiaceae',
 'metasequoia':'Cupressaceae','microtropis':'Celastraceae','miscanthus':'Poaceae',
 'morus':'Moraceae','ohwia':'Fabaceae','ophiopogon':'Asparagaceae',
 'oplismenus':'Poaceae','osmanthus':'Oleaceae','paraprenanthes':'Asteraceae',
 'parathelypteris':'Thelypteridaceae','parthenocissus':'Vitaceae','pellionia':'Urticaceae',
 'perkinsiodendron':'Styracaceae','photinia':'Rosaceae','phyllostachys':'Poaceae',
 'pilea':'Urticaceae','pileostegia':'Hydrangeaceae','pinus':'Pinaceae',
 'plagiogyria':'Plagiogyriaceae','polystichum':'Dryopteridaceae','prunus':'Rosaceae',
 'pteridium':'Dennstaedtiaceae','pyrrosia':'Polypodiaceae','quercus':'Fagaceae',
 'reynoutria':'Polygonaceae','rhaphiolepis':'Rosaceae','rhododendron':'Ericaceae',
 'rosa':'Rosaceae','rubus':'Rosaceae','sabia':'Sabiaceae','sageretia':'Rhamnaceae',
 'sargentodoxa':'Lardizabalaceae','schima':'Theaceae','schoepfia':'Schoepfiaceae',
 'scutellaria':'Lamiaceae','sedum':'Crassulaceae','selaginella':'Selaginellaceae',
 'semiaquilegia':'Ranunculaceae','sinomenium':'Menispermaceae','smilax':'Smilacaceae',
 'solanum':'Solanaceae','spiraea':'Rosaceae','stachyurus':'Stachyuraceae',
 'stauntonia':'Lardizabalaceae','styphnolobium':'Fabaceae','styrax':'Styracaceae',
 'symplocos':'Symplocaceae','syzygium':'Myrtaceae','taxus':'Taxaceae',
 'ternstroemia':'Pentaphylacaceae','tetradium':'Rutaceae','trachelospermum':'Apocynaceae',
 'trachycarpus':'Arecaceae','tsuga':'Pinaceae','ulmus':'Ulmaceae',
 'uncaria':'Rubiaceae','vaccinium':'Ericaceae','veronicastrum':'Plantaginaceae',
 'viburnum':'Adoxaceae','viola':'Violaceae','vitis':'Vitaceae',
 'wisteria':'Fabaceae','woodwardia':'Blechnaceae','seq_695':'Schisandraceae',
}

def family_of(sseqid):
    p = sseqid.split('_')
    genus = p[2].lower() if len(p) >= 4 else sseqid.lower()
    return G2F.get(genus)

def load_sample(sample):
    """Read blast (filtered) + tpm once -> {fam: tpm} for pident 90 and 95.
    Returns dict thr -> {fam: tpm}. Also returns all_contigs_tpm_sum (unused)."""
    blast_path = os.path.join(RAW, f"{sample}_blast_annotation.csv")
    tpm_path = os.path.join(RAW, f"{sample}_tpm_abundance.csv")
    # read blast: (qseqid, pident, length, family)
    hits = []
    with open(blast_path, encoding='utf-8-sig', errors='replace') as fh:
        rd = csv.reader(fh)
        header = next(rd)
        for row in rd:
            if len(row) < 4:
                continue
            q, sseqid = row[0], row[1]
            try:
                pident = float(row[2]); alen = float(row[3])
            except ValueError:
                continue
            fam = family_of(sseqid)
            hits.append((q, pident, alen, fam))
    # read tpm only for contigs that appear in blast (small subset)
    needed = {q for q, _, _, fam in hits if fam is not None}
    tpm = {}
    with open(tpm_path, encoding='utf-8-sig', errors='replace') as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            if len(row) < 4:
                continue
            c = row[0]
            if c in needed:
                try:
                    tpm[c] = float(row[3])
                except ValueError:
                    pass
    out = {95: {}, 90: {}}
    for q, pident, alen, fam in hits:
        if fam is None or q not in tpm:
            continue
        for thr in (95, 90):
            if pident >= thr and alen >= 150:
                out[thr][fam] = out[thr].get(fam, 0.0) + tpm[q]
    return out

def build_table(samples, pident_thr, per_sample_rra):
    """per_sample_rra: dict sample -> {fam: rra}"""
    retained = [s for s in samples if per_sample_rra.get(s)]
    n = len(retained)
    fams = {}
    for s in retained:
        for f, v in per_sample_rra[s].items():
            fams.setdefault(f, {'present': 0, 'rra_sum': 0.0})
            fams[f]['present'] += 1
            fams[f]['rra_sum'] += v
    out = []
    for f, d in fams.items():
        out.append((f, d['present'], d['present'] / n * 100.0, d['rra_sum'] / n))
    out.sort(key=lambda x: -x[3])
    return out, retained, n

if __name__ == '__main__':
    # single pass: load every sample once
    fam_tpm = {}  # sample -> {thr -> {fam: tpm}}
    for s in REEVES + BLACK:
        fam_tpm[s] = load_sample(s)

    results = {}
    for thr in (95, 90):
        for d73 in ('',):
                                    def rra(samples):
                out = {}
                for s in samples:
                    ft = fam_tpm[s][thr]
                    tot = sum(ft.values())
                    out[s] = {f: v / tot * 100.0 for f, v in ft.items()} if tot > 0 else {}
                return out
            rra_r = rra(sr); rra_b = rra(sb)
            tr, rr, nr = build_table(sr, thr, rra_r)
            tb, rb, nb = build_table(sb, thr, rra_b)
            results[(thr, d73)] = {'Reeves': (tr, nr, rr), 'black': (tb, nb, rb)}
            print(f"\n===== ASM{thr} {d73} : Reeves {nr}/{len(sr)} retained, black {nb}/{len(sb)} retained =====")
            print("--- Reeves top 12 ---")
            for f, pres, foo, mrra in tr[:12]:
                print(f"  {f:20s} FOO {foo:6.2f}% ({pres}/{nr}) mRRA {mrra:6.2f}%")
            print("--- black top 12 ---")
            for f, pres, foo, mrra in tb[:12]:
                print(f"  {f:20s} FOO {foo:6.2f}% ({pres}/{nb}) mRRA {mrra:6.2f}%")

    # write Diet_*.csv in the format expected by the R script
    def write_csv(thr, d73):
        tok = ''
        for sp, tag in (('Reeves', 'Reeves'), ('black', 'black')):
            tr, nr, rr = results[(thr, d73)][sp]
            path = os.path.join(OUT, f"Diet_ASM{thr}_{tok}_{tag}.csv")
            with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
                w = csv.writer(fh)
                w.writerow(['Family', 'mRRA_pct', 'FOO_pct', 'Occurrence', 'S_valid'])
                for f, pres, foo, mrra in tr:
                    w.writerow([f, round(mrra, 6), round(foo, 6), pres, nr])
            print("wrote", path)
    for thr in (95, 90):
        for d73 in ('',):
            write_csv(thr, d73)
    print("DONE")
