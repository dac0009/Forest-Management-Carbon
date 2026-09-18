# Phase 0 QC log — tree dataset rebuild (Sept 18, 2026)

Source: project file `Above_ground_biomass_data_` (pasted per-tree inventory).

| Gate | Result | Status |
|---|---|---|
| Raw rows, all watersheds | 1,058 (WS5 321, WS2 231, F38-1 189, WS4 160, WS1 157) | matches prior reconstruction |
| Core rows WS2/WS4/WS5 | 712 | matches |
| WS4 plot count | 7 labels; `WS4-A-11` = 9 trees, max DBH 16.3 in, 226 kg total vs 1,170–2,174 kg for full plots | **excluded per Dom (partial plot); 6 full plots retained** |
| Rows after exclusion | 703 | — |
| Null DBH dropped | WS2 0, WS4 2, WS5 7 (9 total) | matches prior counts |
| Units | raw `diameter` in **inches** (1.0–36.4); converted ×2.54 → 2.5–92.5 cm | 92.5 cm max in 110-yr unmanaged stand plausible; 36.4 cm max is not |
| Legacy `agb` column | median(new/legacy) = 11.2× | consistent with legacy = allometry run on inches-as-cm; legacy column NOT used |
| Allometry | Jenkins et al. (2003): hard maple/oak/hickory/beech (β₀=−2.0127, β₁=2.4342); soft maple/birch (−1.9123, 2.3651); mixed hardwood (−2.4800, 2.4835). Stems: 404 hard / 159 soft / 131 mixed | species→group mapping in trees_clean.csv (`jgrp`) |
| Carbon fraction | 0.48 × biomass | consistent with prior pipeline |
| Plot area | 0.1 acre (404.686 m²) assumed; radius not in writing (data-manager confirmation verbal) | to be re-validated against 1958–2020 100% inventories in Phase 2 (prior check: WS5A 73.5 vs 72.3 Mg C/ha) |
| Reproduction check | REF 177.4 ± 20.7, DLC 91.0 ± 6.2, STS 92.2 ± 17.5 Mg C/ha | matches last chat's verified 177/89/91 within rounding/mapping tolerance |

Files: trees_clean.csv (694 stems w/ valid DBH), stand_by_plot.csv (18 plots × 14 vars), stand_stats.csv (11 variables × 4 tests + Dunnett).
Key results: AGB-C ANOVA P=0.002 (Welch 0.015, KW 0.011; Levene 0.029 so Welch/KW are the honest headline); Dunnett DLC P=0.003, STS P=0.004. Stems ≥27.9 cm P=0.42 (n.s.). Zero trees ≥50 cm in DLC (REF: 17 trees, 71% of pooled AGB; STS: 10 trees). Largest tree: REF 92.5, DLC 47.0, STS 71.4 cm.
