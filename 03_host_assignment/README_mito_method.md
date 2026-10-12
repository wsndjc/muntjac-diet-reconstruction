# 麂属粪样线粒体宿主判定 —— 服务器运行包

## 一、为什么必须换方法

稿件 §2.4 现用**全基因组映射率**判定宿主（取映射率高者）。该方法在两个参考组装之间不可靠。

### 参考基因组的实际参数（NCBI 官方）

| 参数 | 小麂 GCA_020226045.1 | 黑麂 GCA_020276665.1 |
|---|---|---|
| **染色体数** | **23** | **4** |
| **Scaffold N50** | **113,316,579** | **646,574,120**（5.71×）|
| Contig N50 | 18,600,000 | 25,225,000 |
| Scaffold 数 | 460 | 836 |
| 基因组总长 | 2,494,014,298 bp | 2,489,495,709 bp |
| 覆盖度 | 76× | 122× |

麂属是著名的极端染色体重排类群（*M. crinifrons* 2n=4）。在这套参数下，reads 的
比对率主要由**哪套组装更容易比对上**决定，而非样本真正来自哪个物种。

### 实测证据

对 65 个样本的核基因组映射率做统计：

| 指标 | 实测值 |
|---|---|
| 判定裕度中位数 | **0.02 个百分点** |
| 裕度 < 0.05 pp 的样本 | **52/65（80%）** |
| 两率完全相等 | **7 个** |
| 小麂率更高 | **46/65（71%）** |

最后一行是关键：若两参考质量相当，正负应各约 50%。71% 倒向一侧，说明存在
系统性偏倚。被判为黑麂的 13 个样本中，11 个裕度 ≤0.77 pp；D31 两率完全相同
（12.07% vs 12.07%）却仍被判为黑麂。

### 文献依据

> **Atağ G, Vural KB, Kaptan D, et al. (2023)** MTaxi: A comparative tool for taxon
> identification of ultra low coverage ancient genomes. *Open Research Europe* 2:100.
> doi:10.12688/openreseurope.14936.3

该文用绵羊/山羊作为模型，**实测了同一个决策规则**（比对到两个候选物种的核参考基因组、
数比对上的碱基数），结论是：

> "Out of 9 sheep datasets with known species identity, **only 4 showed a higher number
> of bases aligning to the sheep reference** relative to the goat reference."
>
> "**taxon identification using only the number of aligned bases in comparative
> alignment is not reliable**"
>
> 原因："could be related to variability in **reference genome qualities**（the N50 of
> the goat genome ARS1 is 26,244,591 while that for the sheep genome Oar_v3.1 is 40,376）"

绵羊与山羊的分化程度**远大于**同属的两种麂；绵羊/山羊尚且 9 例判错 5 例。

> **Borry M, Cordova B, Perri A, et al. (2020)** CoproID predicts the source of
> coprolites and paleofeces using microbiome composition and host DNA content.
> *PeerJ* 8:e9001. doi:10.7717/peerj.9001

面对完全相同的场景（粪便、多个候选宿主基因组、需判定来源物种），CoproID **刻意不做
原始映射数比较**：把能比上多个基因组的 reads 标记为 "common reads" 并丢弃，再按
基因组大小与物种特异的宿主 DNA 内源比例归一化。原因是内源宿主 DNA 比例是物种属性，
跨物种相差一个数量级（人类 1.24% vs 狗 0.11%，Kruskal–Wallis p < 0.0001）。

---

## 二、本方法：read 碱基投票

**核心思想：不数 reads 总数，读每条 read 在种间差异位点上的碱基。**

种间差异位点由两条线粒体参考做全局比对得出：

| 项目 | 数值 |
|---|---|
| 线粒体基因组 | 小麂 16,354 bp / 黑麂 16,357 bp |
| 序列相似度 | **94.21%** |
| **种间差异位点** | **947 个** |
| ├ 转换 transition | 864 |
| └ 颠换 transversion | 83（MTaxi 严格模式只用这一类）|

一条 150 bp 的 read，平均覆盖 150 × 5.79% ≈ **8.7 个**差异位点：

| | 错配数 |
|---|---|
| read 比对上**正确**参考 | ~0.15（仅测序错误，Q30）|
| 同一 read 比对上**错误**参考 | **~8.7** |
| **分离度** | **约 58 倍** |

因此即便发生**交叉比对**（read 同时能比上两个参考），read 自身携带的碱基仍然可靠。
原方法把 84% 落在保守区的 reads 计入总数——这些 reads 不携带任何判别信息。

**注意**：MTaxi 只用颠换，是因为古 DNA 有 C→T 脱氨基损伤。本研究样本是现代粪便
DNA，无此损伤，**转换可安全使用**，判别位点因此从 83 增至 947（11 倍）。脚本同时
输出两种口径，可互为校验。

### 交叉验证

用麂属其他物种验证比对可靠性：

| 比较 | 分化度 |
|---|---|
| 小麂 vs 小麂（不同个体 MZ895085）| **0.59%** ← 种内水平 |
| 小麂 vs 黑麂（本研究）| **5.79%** |
| 小麂 vs 赤麂（外群）| 6.02% |

种内 0.59% 与种间 5.79% 相差约十倍，符合预期，证明参考序列与标签正确。

---

## 三、文件清单

```
05_mito_hostID.sbatch                     主流程（SLURM）
refs/muntjac_mito_combined.fna            合并参考（两条线粒体）
sites_all.tsv                             947 个位点，两染色体共 1894 行
sites_all.bed                             同上，BED 格式（mpileup -l 用）
sites_transversion_only.tsv               仅 83 个颠换位点
sites_transversion_only.bed
meta.json                                 位点统计
README.md                                 本文件
```

位点表列名：`chrom  pos  base_reevesi  base_crinifrons  kind`

合并参考含两条染色体，同一分化事件在两个坐标系各出现一次。位点表统一使用
物种命名（而非 self/alt），因此无论 read 落在哪条染色体上，其碱基都能被正确
归类，不会互相抵消。

---

## 四、运行

```bash
# 1) 上传本目录到服务器家目录
#    ~/mito_id_bundle/

# 2) 确认 clean reads 位置（脚本默认 ~/260412）
#    支持 *_R1.fq.gz / *_1.fq.gz 两种命名

# 3) 建立输出目录并提交
cd ~/mito_id_bundle
mkdir -p slurm_logs
sbatch 05_mito_hostID.sbatch

# 4) 实时查看
tail -f slurm_logs/mito_*.out
```

若 clean reads 不在 `~/260412`，修改脚本中的 `CLEAN_DIR`。

**依赖**：`conda activate genome_index` 环境（含 bwa、samtools、python3）。
python3 仅用标准库，无需额外安装。

### 资源与耗时

| 项目 | 设置 |
|---|---|
| CPU | `-c 24` |
| 内存 | `--mem 32G`（参考仅 16.4 kb，无需 128G）|
| 墙钟 | `-t 06:00:00` |

按每样本 10 GB（压缩后，约 1.14 亿 reads）、共 65 个样本估算：

| 并行度 | 每样本线程 | 预计耗时 |
|---|---|---|
| 串行（默认 1）| 24 | 3.1 小时 |
| **6（脚本默认）** | **4** | **0.5 小时** |
| 8 | 3 | 0.4 小时 |

**并行度可用环境变量调整**：

```bash
PARALLEL=8 sbatch 05_mito_hostID.sbatch
```

### 为什么并行度是关键

线粒体参考仅 16.4 kb，**比对本身几乎不耗时**；真正的瓶颈是 `bwa mem` 内部
**单线程 gzip 解压**，而 `-t 24` 只加速比对。若一次只跑一个样本并给满 24 线程，
则 23 个核全程闲置等待解压。改为 6 个样本并行、每样本 4 线程，可近线性加速。

脚本另含一项关键优化：**在 `samtools sort` 之前用 `-F 4` 丢弃未比对 reads**。
`bwa mem` 默认把未比对 reads 也写进 SAM，1.14 亿 reads 会产生约 40 GB SAM，
而其中能落在 947 个判别位点上的只有几百到几千条；65 个样本累计约 2.6 TB
废数据流。过滤后每样本 BAM 仅数十 MB。

---

## 五、输出

```
~/mito_hostID_results/
├── host_call_allsites.tsv        全部 947 位点的判定
├── host_call_transversion.tsv    仅 83 颠换位点的判定（MTaxi 严格口径）
├── log/<sample>.mpileup          判别位点的 pileup（可用于复核）
└── vote_sites.py                 投票脚本
```

结果表列：

| 列 | 含义 |
|---|---|
| `SampleID` | 样本名 |
| `votes_reevesi` | 支持小麂的碱基投票数 |
| `votes_crinifrons` | 支持黑麂的碱基投票数 |
| `ambiguous` | 两参考之外的碱基（测序错误）|
| `n_votes` | 有效投票总数 |
| `sites_hit` | 命中的判别位点数 |
| `p_value` | 精确二项检验双侧 p（H₀: p=0.5）|
| `call` | `reevesi` / `crinifrons` / `ambiguous` / `NA` |

BAM 在 mpileup 后自动清理以节省空间（mpileup 已含所需全部信息）。

---

## 六、预期结果与功效

基于现有 67 个样本的实际覆盖度（平均深度中位数约 4×，最高 81×）做的功效分析：

| 假设的投票噪声 | 功效 ≥80% 的样本 |
|---|---|
| p=0.95（保守）| **67/67 (100%)** |
| p=0.90（很保守）| **67/67 (100%)** |
| p=0.80（极端保守）| 65/67 (97%) |

即便最弱的样本（E32，仅 20 条线粒体 reads）功效仍有 99.5%。

**对照**：现有全基因组映射率法可判定 **6/67（9%）**。

### 本地模拟验证

已用模拟数据端到端测试投票逻辑：

| 真实物种 | 错读率 | 判定结果 |
|---|---|---|
| reevesi | 0% | reevesi ✓ |
| reevesi | 5% | reevesi ✓ |
| reevesi | 15% | reevesi ✓ |
| crinifrons | 0% | crinifrons ✓ |
| crinifrons | 5% | crinifrons ✓ |
| crinifrons | 15% | crinifrons ✓ |

含 `^]` 起始标记、`$` 结束标记、`+2AG` indel 指示的解析均已验证。

---

## 七、跑完后必须核对的三件事

**1. 样本数对齐。** Table S2 是 **65** 个样本，而线粒体汇总表是 **67** 个
。分析集为 65 个样本。

**2. 与核基因组判定的一致性。** 现有核判定为 小麂 46 / 黑麂 12 / 平局 7，
其中 **21 个样本与线粒体（旧法）冲突**。以下两个尤其需要解释：

| 样本 | 核判定（裕度）| 线粒体旧法（裕度）|
|---|---|---|
| **E81** | 黑麂（5.97 pp）| 小麂（9.5%）|
| **E21** | 小麂（4.54 pp）| 黑麂（3.3%）|

**3. 黑麂样本 n=6 的确认。** D17 / D21 / D24 / D28 / D29 / D78 是种间食性比较的
实际数据来源。线粒体旧法已全部支持黑麂，新法应确认这一点。

补充：归属最可疑的三个黑麂样本（D31、D27、E81）目前**全部零检出**，未进入食性
计算——若新法结果有变，需重新评估。

---

## 八、结果回传

跑完后请回传：

1. `~/mito_hostID_results/host_call_allsites.tsv`
2. `~/mito_hostID_results/host_call_transversion.tsv`

用于与核基因组判定、食性数据做三方交叉核对。

---

## 附：建议写入 Methods 的表述

> Host identity was determined by competitive mapping of clean reads to the complete
> mitochondrial genomes of *Muntiacus reevesi* (NC_004069.1) and *M. crinifrons*
> (NC_004577.1). Because the two reference assemblies differ substantially in
> chromosome number (2n = 46 vs 2n = 4) and scaffold N50 (113 Mb vs 647 Mb), whole-genome
> mapping rates are governed by reference quality rather than by the true source species,
> and are not a reliable discriminator for closely related taxa (Atağ et al., 2023;
> Borry et al., 2020). Each sample was therefore assigned by the bases carried by
> individual reads at the 947 interspecific diagnostic sites (5.79% sequence divergence;
> 947 substitutions, of which 83 are transversions), tested with an exact two-sided
> binomial test (H₀: p = 0.5). A read of 150 bp spans on average 8.7 diagnostic sites,
> giving an approximately 58-fold separation between correct and incorrect reference
> assignment, so that read-level assignment remains reliable even where cross-mapping
> occurs. Samples with p ≥ 0.05 after Bonferroni correction were reported as ambiguous.
