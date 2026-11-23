# Project log: microglial gene networks in Alzheimer's disease (SEA-AD, MTG)

*Working notes, not the official report. English copy of `project_log.md`. Last updated: 21 Nov 2026.*

## 0. Question

Which microglial genes are **directly** linked to the severity of Alzheimer's pathology (CPS), and which only through
other genes — and how far do the three network reconstruction methods from the course (**ARACNE, PC, MIIC**) agree on this?

## 1. Data

- **SEA-AD** (Allen Institute), immune nuclei, June 2026 release:
  `SEAAD_Immune_..._final-nuclei.2026-06-22.h5ad` (3.27 GB): **240,651 nuclei × 36,601 genes**, 84 donors, 10 brain regions.
- Raw counts are in `layers/UMIs` (integers); `X` holds normalized data, which we **do not use** (mentor's advice).
- **CPS** comes from `sea-ad_global_and_local_cps_050126.xlsx`, sheet "Travaglini 2026"; we use `CPS_Local` for MTG
  (0 = almost no pathology, 1 = most severe in the cohort; **one value per donor**).
- Donor metadata descriptions: `sea-ad_donorindex_dataelementbriefdescriptions_063026.pdf`.

Scripts: `R/01`–`R/08`; `.h5ad` files are read without Python via `R/read_h5ad.R`.

## 2. Data exploration (`R/01_explore_data.R` → `results/01_explore_output.txt`)

- Little is captured per nucleus: median **~2,600 genes and ~5,300 molecules**. The most frequent genes are the nuclear RNAs MALAT1 and NEAT1.
- Many key genes are mostly zero: TREM2 is non-zero in ~23% of nuclei, GPNMB ~15%, LPL/CD9 ~6%.
- 127 metadata columns: donor, region, supertype (`Supertype`), Braak, Thal, ADNC, APOE, age, sex, etc.
  **There is no CPS column in the file** — it comes from a separate table.

## 3. Cell selection (`R/02_mtg_microglia.R` → `results/02_mtg_microglia_output.txt`)

Funnel: 240,651 → MTG 36,686 → without monocytes/lymphocytes 35,907 → passed QC **34,679 nuclei, 83 donors**.

Our interpretation of the supertypes from marker genes (fraction of nuclei in which the gene is detected):

| Supertype | Nuclei | Interpretation | Why |
|---|---|---|---|
| Micro-PVM_1 | 1,335 | perivascular macrophages → **excluded** | MRC1 0.73, LYVE1 0.64, CD163 0.91; P2RY12 0.21 |
| Micro-PVM_2 | 20,386 | homeostatic microglia | P2RY12 0.77, CX3CR1 0.58, GPNMB 0.10 |
| Micro-PVM_2_1, _2_3 | 8,294 | close to homeostatic | profile similar to PVM_2 |
| Micro-PVM_3-SEAAD | 4,259 | reactive (DAM) | GPNMB 0.43, higher TREM2, lower P2RY12 |
| Micro-PVM_4-SEAAD | 405 | strongly activated | SPP1 0.95, GPNMB 0.49, P2RY12 0.14 |

This matches Gabitto et al. 2024 (Fig. 6a: Micro-PVM_2 "homeostatic", Micro-PVM_3 "DAM").

CPS sanity check: median `CPS_Local` increases with ADNC (Not AD 0.06 → Low 0.14 → Intermediate 0.43 → High 0.68)
and with Braak stage (0.08 → 0.70). "Severely affected donors" (10 people) have severe pathology (median 0.69), not poor data quality.

## 4. Metacells (`R/03_build_dataset.R`)

**Why.** Only a small fraction of molecules is captured in a nucleus, so zeros are measurement chance. Network methods look for
coordinated changes of genes from row to row, and random zeros hide that coordination. Summing 30 nuclei averages out the noise
(TREM2 non-zero: 23% of nuclei → 99% of metacells).

**How.** Within each "donor × supertype" pair, nuclei are randomly split into groups of 30 and their raw counts are **summed**
(a sum of integers is still raw counts — no normalization). Without Micro-PVM_1: 33,344 nuclei → **983 metacells, 81 donors**
(1 to 32 metacells per donor).

**Technical node `total_UMIs`.** Metacells contain 78k to 374k molecules; sequencing depth raises all genes at once
(gene–depth correlation up to 0.78). Instead of normalizing, depth is included as a network variable.

## 5. Gene selection (`R/03` → `R/04_select_genes.R` → `results/04_selected_genes.csv`)

Filter: gene non-zero in ≥80% of metacells. **48 genes + CPS_Local + total_UMIs = 50 variables.**

| Source | Group | Genes |
|---|---|---|
| Biology (34) | homeostatic microglia (13) | P2RY12, P2RY13, CX3CR1, TMEM119, CSF1R, SALL1, GPR34, OLFML3, HEXB, SELPLG, TGFBR1, SLC2A5, CST3 |
| | disease response (10) | TREM2, TYROBP, APOE, SPP1, GPNMB, CD74, HLA-DRA, C1QA, C1QB, CD163 |
| | AD risk genes (8) | PLCG2, INPP5D, MS4A6A, PICALM, SORL1, BIN1, CD33, ABI3 |
| | transcription factors (3) | SPI1, IRF8, MEF2C |
| Data (14) | higher in reactive (PVM_3/4) vs homeostatic (PVM_2) | TPRG1, DIRC3, FAM20A, TGFBI, PLIN2, PTPRG, FGR, CD83, TBC1D2, IQGAP2, GAS2L3, ZNF804A, PLXNC1, CDK14 |

History: the first version, "most variable genes" (`network_dataset_v1_hvg.csv`), gave mostly unnamed non-coding RNAs,
stress genes (HSPA1A/B, DUSP1), genes of other cell types (FLT1) and genotype-driven genes (ARL17B) — rejected.
Data-driven genes "higher in homeostatic" turned out to be contamination from other cells (MMRN1, TNFRSF13C, CLEC9A) — replaced
by known homeostatic microglia genes. Normalization was used only for ranking. **Genes were not selected by their association with CPS.**

## 6. Networks on metacells (`R/05_networks.R` → `results/05_networks_output.txt`, `results/networks/`)

Plots: [ARACNE](../results/networks/net_aracne.png) · [PC](../results/networks/net_pc.png) ·
[MIIC](../results/networks/net_miic.png) · [consensus](../results/networks/net_consensus.png).
All edges and which methods found them: `results/networks/edges_comparison.csv`.

- Edges: ARACNE 71, PC 28 (see section 10 — the PC count is being recomputed), MIIC 111.
- **21 consensus edges** (all three methods), mostly known biology:
  C1QA–C1QB (two parts of the same complement complex), CD74–HLA-DRA (MHC-II and its chaperone),
  P2RY12–CX3CR1, P2RY12–SALL1, GPR34–OLFML3 (homeostatic microglia), APOE–TYROBP, APOE–C1QA (disease response),
  MS4A6A–PLCG2 (AD risk genes), GPNMB–IQGAP2, CD163–TGFBI/FGR, TPRG1–ZNF804A (reactive microglia),
  plus MEF2C/PICALM/PLCG2 — total_UMIs (depth).
- `total_UMIs` absorbed the genes most dependent on depth (MEF2C, PICALM, PLCG2 ~0.75), as intended.
- **CPS_Local:** ARACNE 1 neighbour (PTPRG), PC 0, **MIIC 19** — suspicious (many of these neighbours have ~0 correlation
  with CPS; MIIC marked 6 of these edges as "latent common cause").

## 7. Experiment B: single nuclei (`R/06a_build_nuclei_dataset.R` → `results/networks_nuclei/`)

Same 50 variables, but each row is one nucleus with raw counts; 12 random nuclei per donor (996 rows).
Report: `results/networks_nuclei/networks_report.txt`.
Plots: [ARACNE](../results/networks_nuclei/net_aracne.png) · [PC](../results/networks_nuclei/net_pc.png) ·
[MIIC](../results/networks_nuclei/net_miic.png) · [consensus](../results/networks_nuclei/net_consensus.png).

| | Metacells | Single nuclei |
|---|---|---|
| Edges ARACNE / PC / MIIC | 71 / 28* / 111 | 57 / 29* / 84 |
| Consensus edges | 21 | 18 |
| CPS neighbours (ARACNE / PC / MIIC) | 1 / 0* / **19** | **2 / 2* / 2** (PTPRG, TGFBI) |
| MIIC "latent cause" edges | 9 | 0 |
| total_UMIs neighbours (ARACNE / PC / MIIC) | 8 / 3* / 12 | 26 / 8* / 29 |

\* PC numbers are being recomputed (section 10).

- Gene–gene correlations weakened because of zeros (C1QA–C1QB 0.81 → 0.47). Of the 18 gene–gene consensus edges on metacells,
  **6 are found by all methods on nuclei as well** (APOE–C1QA, APOE–TYROBP, CD74–HLA-DRA, CX3CR1–P2RY12, GPNMB–IQGAP2, DIRC3–IQGAP2) —
  the most reliable links in the project; 5 are found by 1–2 methods; 7 disappeared (mostly sparse genes: OLFML3 in 7% of nuclei, GAS2L3 8%).
- Depth matters more in nuclei (1,342 to 36,686 molecules per nucleus, a 27-fold range).

## 8. CPS stability and control (`R/06b_cps_stability.R` → `results/06b_cps_stability_output.txt`, `_runs.csv`)

First, what we "did to the data": donors were **not mixed** — every metacell comes from the nuclei of a single person.
The problem is different: CPS is the same for all rows of a donor (983 rows but 81 distinct values), so for the algorithm CPS acts as
a **donor ID**; metacells remove noise and make differences between people (genetics, circumstances of death, tissue handling)
very sharp. MIIC captures any dependence, not only monotone trends, and therefore links CPS to everything that differs between people.

Test: 5 different samples of 12 nuclei per donor; on each, the real CPS and CPS randomly shuffled **between donors**
(each donor keeps one value, but someone else's: the link with disease is destroyed, the "donor ID" is kept). The same shuffle on metacells.

| Variant | ARACNE | PC | MIIC |
|---|---|---|---|
| Nuclei, real CPS (mean neighbours, 5 samples) | 2.2 | 3.0 | 2.6 |
| Nuclei, shuffled CPS | 0.2 | 0.8 | 0.8 |
| Metacells, shuffled CPS (2 repeats) | 2.5 | 1.5 | **18 and 18** |
| Metacells, real CPS (from section 6) | 1 | 0* | 19 |

**Conclusions:**
1. **The 19 MIIC CPS neighbours on metacells are a donor artefact.** With shuffled CPS, MIIC finds the same ~18 neighbours
   (APOE, BIN1, C1QB, CD74, CSF1R, P2RY12, SPP1 …), although there is no link with disease by construction.
2. **PTPRG is a robust CPS neighbour:** ARACNE 4/5, PC 5/5, MIIC 5/5 samples; never appears with shuffled CPS.
3. **IQGAP2** is the second most stable (ARACNE 3/5, PC 4/5, MIIC 2/5). **BIN1** — 2/5 (negative association).
4. **TGFBI is not robust:** in the first sample (06a) all three methods found it, but in 06b only 1 of 15 runs did.
   It was a feature of that particular nucleus sample.
5. Background noise: with random CPS on nuclei, 0–2 "neighbours" appear per run (different genes each time) — a reference
   for how many false neighbours to expect.

## 9. Comparison with the literature (`results/07_literature_gene_comparison.csv`)

Papers:
- **Gabitto et al. 2024, Nat Neurosci** (main SEA-AD paper), Fig. 6 — microglia: supertype Micro-PVM_3 = DAM, its proportion
  increases with CPS; supertype markers (Fig. 6c); genes increasing early and late along CPS (Fig. 6d); 4 early-response TFs (Fig. 6e).
- **Adeoye et al. 2025, Aging Dis** — the same SEA-AD MTG data, hdWGCNA on metacells; microglial modules MIC-M1…M5 (Fig. 8A).
  Disease as groups (no / early / late pathology), not continuous CPS; no explicit gene–gene edges in the text, only modules and hubs (Suppl. Table 13).
- **Hawrylycz et al. 2024, Nat Aging** — description of the SEA-AD resource, no genes; used only as a data reference.

**Agreements:**
- **PTPRG, TPRG1, ZNF804A** — our genes selected *from the data* — are among the DAM supertype markers in Gabitto (Fig. 6c).
  In our network they are connected to each other (TPRG1–ZNF804A: all three methods on metacells; PTPRG–ZNF804A: ARACNE+MIIC).
  **PTPRG is our robust CPS neighbour**, consistent with the increase of the DAM supertype with disease in Gabitto.
- Gabitto, increasing early along CPS: **CSF1R, CD74, C1QA, C1QB** (in our set); late: **APOE** (in our set).
- Adeoye MIC-M2 (up in AD): **APOE, CST3** — we have an **APOE–CST3** edge (ARACNE+MIIC, on metacells and nuclei): the module is reproduced.
- Adeoye MIC-M5: GPNMB; MIC-M4: PICALM; shared hubs: MEF2C, TGFBR1.

**Not in our set although present in the data:** almost all genes from the papers pass our filter (only sparse ones fail: IL1B, FCGR3B, LYZ, LINGO1, etc.).
Candidates that increase with CPS both in the papers and in our data: **CTSD** (ρ≈+0.33), NHSL1, TANC2, CADM1, WIPF3, MITF, STARD13, CD81.

**Discrepancies and explanations:**
- Gabitto's early genes (CSF1R, C1QA/B, IRF1/7, Fc receptors, TFs RUNX1/IKZF1/NFATC2/MAF) barely correlate with CPS in our data
  (ρ from −0.2 to +0.17). In Gabitto they rise **at the early stage** and then plateau; a single CPS node (monotone association
  over the whole range) blurs such non-linear changes.
- Correlations in table 07 were computed on normalized metacells and are overconfident (one CPS value per donor) — for orientation only.

## 10. Known issues and caveats

- **Bug in PC edge counting (found 21 Nov).** In `05_networks.R` a PC edge was counted only when the adjacency matrix held the value 1;
  edges marked while resolving orientation conflicts were lost. 06b counted all edges — there PC has ~90–105 edges,
  not ~28. Fixed; **the 05 results for metacells and nuclei must be rerun** and the PC numbers in sections 6–7 updated.
- CPS is a donor-level variable: there are ~81–83 independent observations for CPS links, not ~1,000.
- An edge is not causation; PC/MIIC arrows are hypotheses.
- Limitation of a single CPS node: non-linear (early) changes are poorly captured.
- The choice of 48 genes out of 36,601 affects the result → experiments with other gene sets (section 11).

## 11. Next (`R/08_gene_sets.R`)

Two alternative gene sets, each on metacells and on nuclei (12 per donor, the same sample as 06b rep 1):
1. **"Literature"** — only genes from Gabitto Fig. 6 and Adeoye Fig. 8 (~70 after filtering).
2. **"Data only"** — genes with the highest variability **between metacells of the same donor** (60 by default):
   this selects cell-state programmes rather than differences between people.

## Files

| What | Where |
|---|---|
| Network tables (metacells / nuclei) | `data/processed/network_dataset.csv`, `network_dataset_nuclei.csv` (+ `.tsv` for the MIIC server) |
| Metacell / nucleus metadata | `data/processed/metacell_info.csv`, `nucleus_info.csv` |
| Selected genes | `results/04_selected_genes.csv` |
| Networks | `results/networks/`, `results/networks_nuclei/` |
| CPS stability | `results/06b_cps_stability_output.txt`, `results/06b_cps_stability_runs.csv` |
| Literature comparison | `results/07_literature_gene_comparison.csv` |
