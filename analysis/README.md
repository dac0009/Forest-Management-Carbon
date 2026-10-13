# Fernow carbon notebooks — ArcGIS Pro / Jupyter

Copy these five files into the **root of the Forest-Management-Carbon repo** (next to `gdb_export/` and the `WV_*.csv` files).
In ArcGIS Pro: Catalog > add the folder > open each .ipynb. Default `arcgispro-py3` environment; no installs needed.

**FIA first (once, locally):** run `00_FIA_Filter_Local` against your full download (`C:\Users\DCifelli2\Downloads\WV_CSV`).
It writes `WV_CSV_filtered/` (compressed tables + `fia_cond_live` + `fia_cond_seedlings` + `MANIFEST.csv`). Copy that
folder into the repo as `FIA_WV_filtered/`, push it, and set `FIA_DIR = r'FIA_WV_filtered'` in notebook 04.

Then run in order (02 writes `fernow_census_periods.csv`, which 04 overlays):

| Notebook | What it does |
|---|---|
| 00_FIA_Filter_Local | Chunked filter of the full WV FIA download; live-tree C and seedling summaries per condition |
| 01_Fernow_Core_Stats | AGB-C, soil C/N/BD/roots, density fractions (LF-C primary), ecosystem C with POM vs MAOM |
| 02_Fernow_History_Terrain | Harvest-date gate, 1958-59 baseline, 1958-2022 storage vs sequestration, plot representativeness, terrain ANCOVA, lidar check |
| 03_Soil_Mechanisms | QC gates (WS5 soil labels, OH mass), literature-expectation tests, equivalent soil mass |
| 04_FIA_WV_Context | WV FIA: representativeness, regional harvest effect, remeasurement growth rates, measured dead wood/forest floor |

Every notebook has one CONFIG cell at the top — paths and the few analysis decisions live there.
`fernow_stats.py` holds all shared statistics so every notebook uses identical code.
