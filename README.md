# IndepMR-App
An ontology-driven, ancestry-aware interactive R Shiny framework for automated two-sample and multivariable Mendelian Randomization with mediation analysis using public GWAS summary data.

This release corresponds to the version of the software described in the accompanying manuscript submitted to BMC Bioinformatics.

**Archived version (Zenodo):** https://doi.org/10.5281/zenodo.21441281

## Key features

- **Ontology-driven trait search engine** that exploits EFO lineage paths to retrieve GWAS studies beyond simple keyword matching
- **Ancestry independence engine** that classifies each study's population description using a five-group continental hierarchy with per-main group subgroup resolution, correctly handling multi-ancestry populations via umbrella-term detection
- **Fully integrated five-tab workflow**: exposure/outcome study selection, two-sample MR (2SMR), multivariable MR (MVMR), and MR-based mediation analysis, with automated instrument extraction, LD clumping, allele harmonisation, sensitivity diagnostics, and exportable results

## Requirements

- R (≥ 4.2.0)
- Key R dependencies: `shiny`, `TwoSampleMR`, `ieugwasr`, `MVMR`, `DT`, `plotly`, `DiagrammeR`, `openxlsx`, `readxl` (see Additional File 1 in the manuscript for the full list)
- An OpenGWAS JWT token is required for live API access (instrument extraction, LD clumping, outcome retrieval); a local instrument mode is available as a fallback

## File structure

indepmr/
├── app.R                          # Entry point
├── global.R                       # Package loading, data paths, DATA object
├── data/
│   ├── Dataset01.xlsx             # GWAS Catalog metadata
│   ├── Dataset02.rds              # Pre-extracted SNP data
│   └── Dataset03.csv              # EFO ontology token matrix
└── R/
├── 00_utils.R                 # JWT management, API wrappers
├── 01_cache_helpers.R         # On-disk caching
├── 02_load_data.R             # Data loading and standardisation
├── 03_independence_engine.R   # Ancestry NLP pipeline
├── 04_search_helpers.R        # GWAS Catalog search logic
├── 05_search_module.R         # Reusable Shiny search module
├── 06_mr_helpers.R            # MR estimation wrappers
├── 08_local_instruments.R     # Local SNP retrieval
└── server_original.R          # Shiny server function

## License

MIT License

## Citation

If you use IndepMR in your research, please cite:

Yatawara TPB. IndepMR: An Ontology-Driven, Ancestry-Aware Interactive R 
Shiny Framework for Automated Two-Sample and Multivariable Mendelian 
Randomization with Mediation Analysis Using Public GWAS Summary Data. 
BMC Bioinformatics (submitted).

Source code and data archive: https://doi.org/10.5281/zenodo.21441281
