# -*- coding: utf-8 -*-
"""Top-20 / Top-10 Jaccard via a scientifically-defensible ranking:
rank families by MEAN mRRA across the 6 molecular datasets (absent=0), restrict all
7 datasets (Lit + DIR-L/M/S + ASM90/95/99) to top-N, compute pairwise Jaccard.
NEW mitochondrial assignment: 54 Reeves + 11 black."""
import os
import io, sys, csv, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

RAW = os.environ.get('MUNTJAC_RAW', 'raw')  # directory of per-sample blast annotation + TPM files
DIET = os.environ.get('MUNTJAC_DIR', 'dir')  # DIR coverage files
SUPP = os.path.join(os.environ.get("MUNTJAC_DATA", "."), "04_statistics", "family_presence_by_method.csv")

G2F={'acer':'Sapindaceae','actinidia':'Actinidiaceae','adina':'Rubiaceae','adinandra':'Pentaphylacaceae','ainsliaea':'Asteraceae','akebia':'Lardizabalaceae','alangium':'Cornaceae','albizia':'Fabaceae','alniphyllum':'Styracaceae','alyxia':'Apocynaceae','amelanchier':'Rosaceae','ampelopsis':'Vitaceae','anoectochilus':'Orchidaceae','aphananthe':'Cannabaceae','arachniodes':'Dryopteridaceae','aralia':'Araliaceae','ardisia':'Primulaceae','arisaema':'Araceae','artemisia':'Asteraceae','aster':'Asteraceae','berchemia':'Rhamnaceae','betula':'Betulaceae','broussonetia':'Moraceae','callerya':'Fabaceae','callicarpa':'Lamiaceae','camellia':'Theaceae','carex':'Cyperaceae','carpinus':'Betulaceae','castanea':'Fagaceae','castanopsis':'Fagaceae','celastrus':'Celastraceae','celtis':'Cannabaceae','cephalotaxus':'Taxaceae','chamaecyparis':'Cupressaceae','chengiodendron':'Oleaceae','chloranthus':'Chloranthaceae','choerospondias':'Anacardiaceae','cinnamomum':'Lauraceae','cladrastis':'Fabaceae','clematis':'Ranunculaceae','clerodendrum':'Lamiaceae','clethra':'Clethraceae','cleyera':'Pentaphylacaceae','clinopodium':'Lamiaceae','commelina':'Commelinaceae','coptosapelta':'Rubiaceae','cornus':'Cornaceae','corylopsis':'Hamamelidaceae','crataegus':'Rosaceae','cunninghamia':'Cupressaceae','cyclocarya':'Juglandaceae','dalbergia':'Fabaceae','damnacanthus':'Rubiaceae','daphne':'Thymelaeaceae','daphniphyllum':'Daphniphyllaceae','dendropanax':'Araliaceae','deutzia':'Hydrangeaceae','dicranopteris':'Gleicheniaceae','dioscorea':'Dioscoreaceae','diospyros':'Ebenaceae','diplopterygium':'Gleicheniaceae','diplospora':'Rubiaceae','distylium':'Hamamelidaceae','dryopteris':'Dryopteridaceae','elaeagnus':'Elaeagnaceae','elaeocarpus':'Elaeocarpaceae','elatostema':'Urticaceae','eriobotrya':'Rosaceae','erythroxylum':'Erythroxylaceae','euonymus':'Celastraceae','eurya':'Pentaphylacaceae','euscaphis':'Staphyleaceae','fagus':'Fagaceae','ficus':'Moraceae','firmiana':'Malvaceae','frangula':'Rhamnaceae','fraxinus':'Oleaceae','gamblea':'Araliaceae','gardenia':'Rubiaceae','glochidion':'Phyllanthaceae','hedera':'Araliaceae','holboellia':'Lardizabalaceae','houttuynia':'Saururaceae','hovenia':'Rhamnaceae','huperzia':'Lycopodiaceae','hydrangea':'Hydrangeaceae','hydrocotyle':'Araliaceae','idesia':'Salicaceae','ilex':'Aquifoliaceae','illicium':'Schisandraceae','kadsura':'Schisandraceae','lasianthus':'Rubiaceae','lespedeza':'Fabaceae','ligustrum':'Oleaceae','lindera':'Lauraceae','liparis':'Orchidaceae','liquidambar':'Altingiaceae','liriodendron':'Magnoliaceae','lithocarpus':'Fagaceae','litsea':'Lauraceae','lonicera':'Caprifoliaceae','lophatherum':'Poaceae','loropetalum':'Hamamelidaceae','ltea':'Iteaceae','itea':'Iteaceae','lycopodium':'Lycopodiaceae','lyonia':'Ericaceae','lysimachia':'Primulaceae','machilus':'Lauraceae','maclura':'Moraceae','maesa':'Primulaceae','magnolia':'Magnoliaceae','mahonia':'Berberidaceae','mallotus':'Euphorbiaceae','malus':'Rosaceae','meliosma':'Sabiaceae','metasequoia':'Cupressaceae','michelia':'Magnoliaceae','microtropis':'Celastraceae','miscanthus':'Poaceae','morella':'Myricaceae','morus':'Moraceae','myrsine':'Primulaceae','neolitsea':'Lauraceae','nyssa':'Nyssaceae','ohwia':'Fabaceae','ophiopogon':'Asparagaceae','oplismenus':'Poaceae','oreocharis':'Gesneriaceae','osmanthus':'Oleaceae','osmunda':'Osmundaceae','paederia':'Rubiaceae','paraprenanthes':'Asteraceae','parathelypteris':'Thelypteridaceae','parthenocissus':'Vitaceae','paulownia':'Paulowniaceae','pellionia':'Urticaceae','perilla':'Lamiaceae','perkinsiodendron':'Styracaceae','phoebe':'Lauraceae','photinia':'Rosaceae','phyllanthus':'Phyllanthaceae','phyllostachys':'Poaceae','picrasma':'Simaroubaceae','pieris':'Ericaceae','pilea':'Urticaceae','pileostegia':'Hydrangeaceae','pinus':'Pinaceae','pistacia':'Anacardiaceae','pittosporum':'Pittosporaceae','plagiogyria':'Plagiogyriaceae','platycarya':'Juglandaceae','polystichum':'Dryopteridaceae','premna':'Lamiaceae','prunus':'Rosaceae','pteridium':'Dennstaedtiaceae','pterostyrax':'Styracaceae','pyrrosia':'Polypodiaceae','quercus':'Fagaceae','reynoutria':'Polygonaceae','rhaphiolepis':'Rosaceae','rhododendron':'Ericaceae','rhus':'Anacardiaceae','rosa':'Rosaceae','rubus':'Rosaceae','sabia':'Sabiaceae','sageretia':'Rhamnaceae','sargentodoxa':'Lardizabalaceae','schima':'Theaceae','schoepfia':'Schoepfiaceae','scutellaria':'Lamiaceae','sedum':'Crassulaceae','selaginella':'Selaginellaceae','semiaquilegia':'Ranunculaceae','sinoadina':'Rubiaceae','sinomenium':'Menispermaceae','skimmia':'Rutaceae','smilax':'Smilacaceae','solanum':'Solanaceae','sorbus':'Rosaceae','spiraea':'Rosaceae','stachyurus':'Stachyuraceae','stauntonia':'Lardizabalaceae','stranvaesia':'Rosaceae','styphnolobium':'Fabaceae','styrax':'Styracaceae','symplocos':'Symplocaceae','syzygium':'Myrtaceae','tarenna':'Rubiaceae','tashiroea':'Melastomataceae','taxus':'Taxaceae','ternstroemia':'Pentaphylacaceae','tetradium':'Rutaceae','tilia':'Malvaceae','toona':'Meliaceae','toxicodendron':'Anacardiaceae','trachelospermum':'Apocynaceae','trachycarpus':'Arecaceae','tsuga':'Pinaceae','ulmus':'Ulmaceae','uncaria':'Rubiaceae','vaccinium':'Ericaceae','vernicia':'Euphorbiaceae','veronicastrum':'Plantaginaceae','viburnum':'Adoxaceae','viola':'Violaceae','vitis':'Vitaceae','weigela':'Caprifoliaceae','wikstroemia':'Thymelaeaceae','wisteria':'Fabaceae','woodwardia':'Blechnaceae','yushania':'Poaceae','zanthoxylum':'Rutaceae','zelkova':'Ulmaceae','zoysia':'Poaceae','lndocalamus':'Poaceae','indocalamus':'Poaceae','gynostemma':'Cucurbitaceae'}

REEVES=['D11','D14','D15','D18','D19','D23','D3','D34','D38','D46','D48','D50','D51','D5','D52','D53','D57','D61','D6','D62','D63','D66','D67','D68','D69','D70','D71','D7','D72','D86','E20','E21','E23','E24','E32','E33','E35','E36','E37','E42','E46','E56','E57','E65','E75','E76','E81','E82','E83','E84','E85','E89','E91','E92']
BLACK=['D17','D21','D22','D24','D25','D27','D28','D29','D30','D31','D78']

def genus_of(name):
    p=name.lstrip('_').split('_')
    return p[2].lower() if len(p)>=3 and p[0]=='seq' and p[1].isdigit() else p[0].lower()

# ---- DIR: per-sample family RRA (numreads currency) + family set ----
def dir_sample(s, R, C):
    agg={}
    with open(os.path.join(DIET, f'{s}_coverage.txt'), encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#'): continue
            c=line.rstrip('\n').split('\t')
            if len(c)<7: continue
            try: nr=float(c[3]); cv=float(c[5])
            except ValueError: continue
            g=genus_of(c[0])
            if g not in agg: agg[g]=[0.0,0.0]
            agg[g][0]+=nr; agg[g][1]=max(agg[g][1],cv)
    # family numreads (sum over genera passing threshold)
    famreads={}
    for g,(sr,mc) in agg.items():
        f=G2F.get(g)
        if f and sr>=R and mc>=C:
            famreads[f]=famreads.get(f,0.0)+sr
    tot=sum(famreads.values())
    rra={f:v/tot*100.0 for f,v in famreads.items()} if tot>0 else {}
    return rra, set(rra)

# ---- ASM: per-sample family TPM -> RRA + family set ----
def asm_sample(s, thr):
    bp=os.path.join(RAW,f'{s}_blast_annotation.csv'); tp=os.path.join(RAW,f'{s}_tpm_abundance.csv')
    if not os.path.exists(bp) or not os.path.exists(tp): return {}, set()
    tpm={}
    with open(tp,encoding='utf-8-sig',errors='replace') as fh:
        rd=csv.reader(fh); next(rd)
        for row in rd:
            if len(row)<4: continue
            try: tpm[row[0]]=float(row[3])
            except ValueError: pass
    ft={}
    with open(bp,encoding='utf-8-sig',errors='replace') as fh:
        rd=csv.reader(fh); next(rd)
        for row in rd:
            if len(row)<4: continue
            try: pi=float(row[2]); al=float(row[3])
            except ValueError: continue
            f=G2F.get(genus_of(row[1]))
            if f and pi>=thr and al>=150 and row[0] in tpm and tpm[row[0]]>0:
                ft[f]=ft.get(f,0.0)+tpm[row[0]]
    tot=sum(ft.values())
    rra={f:v/tot*100.0 for f,v in ft.items()} if tot>0 else {}
    return rra, set(rra)

# ---- load Lit lists ----
LIT={}
with open(SUPP,encoding='utf-8-sig',errors='replace') as fh:
    rd=csv.reader(fh); next(rd)
    for r in rd:
        if r[1].strip()=='Traditional Survey':
            if 'reevesi' in r[0]: LIT['Reeves']=set(f.strip() for f in r[4].split(','))
            elif 'crinifrons' in r[0]: LIT['black']=set(f.strip() for f in r[4].split(','))

def J(a,b):
    if not a and not b: return 0.0
    return len(a&b)/len(a|b)*100.0

DS=['Lit','DIR-L','DIR-M','DIR-S','ASM90','ASM95','ASM99']
THR_DIR={'DIR-L':(5,0.5),'DIR-M':(20,1.0),'DIR-S':(50,5.0)}
THR_ASM={'ASM90':90,'ASM95':95,'ASM99':99}

def analyze(sp):
    samples = REEVES if sp=='Reeves' else BLACK
    # per-dataset family mRRA (mean over samples of within-sample RRA)
    mrra = {k: {} for k in DS}
    sets = {k: set() for k in DS}
    sets['Lit'] = LIT[sp]
    for k in ['DIR-L','DIR-M','DIR-S']:
        R,C = THR_DIR[k]
        acc={}
        for s in samples:
            rra, fam = dir_sample(s, R, C)
            for f in rra:
                acc.setdefault(f, []).append(rra[f])
        for f, v in acc.items():
            mrra[k][f] = sum(v)/len(samples)   # mean over ALL samples (absent=0)
            sets[k].add(f)
    for k in ['ASM90','ASM95','ASM99']:
        thr = THR_ASM[k]
        acc={}
        for s in samples:
            rra, fam = asm_sample(s, thr)
            for f in rra:
                acc.setdefault(f, []).append(rra[f])
        for f, v in acc.items():
            mrra[k][f] = sum(v)/len(samples)
            sets[k].add(f)
    # rank by mean mRRA across 6 molecular datasets
    allf = set().union(*[sets[k] for k in DS[1:]])
    score = {}
    for f in allf:
        score[f] = sum(mrra[k].get(f, 0.0) for k in DS[1:]) / 6.0
    order = sorted(allf, key=lambda f: -score[f])
    print(f'\n===== {sp} — top-20 families by mean mRRA (across 6 molecular datasets) =====')
    for i, f in enumerate(order[:20], 1):
        print(f'  {i:2d}. {f:20s} mean mRRA={score[f]:6.2f}%   det6={sum(1 for k in DS[1:] if f in sets[k])}')
    return sets, score, order

for sp in ('Reeves','black'):
    sets, score, order = analyze(sp)
    print(f'\n--- {sp}: top-20/top-10 Jaccard (all 7 datasets restricted to top-N) ---')
    for N in (20, 10):
        top = set(order[:N])
        restr = {k: (sets[k] & top) for k in DS}
        vals = [J(restr[DS[i]], restr[DS[j]]) for i in range(7) for j in range(i+1,7)]
        print(f'  top-{N:2d}: min={min(vals):5.1f}%  max={max(vals):5.1f}%   sizes={ {k:len(restr[k]) for k in DS} }')
