# IndepMR-App
Initial release of IndepMR, an ontology-driven, ancestry-aware interactive R Shiny framework for automated two-sample and multivariable Mendelian Randomization with mediation analysis using public GWAS summary data.

This release corresponds to the version of the software described in the accompanying manuscript submitted to BMC Bioinformatics.

Key features:

Ontology-driven trait search engine that exploits EFO lineage paths to retrieve GWAS studies beyond simple keyword matching
Ancestry independence engine that classifies each study's population description using a five-group continental hierarchy with per-main-group subgroup resolution, correctly handling multi-ancestry populations via umbrella-term detection
Fully integrated five-tab workflow: exposure/outcome study selection, two-sample MR (2SMR), multivariable MR (MVMR), and MR-based mediation analysis, with automated instrument extraction, LD clumping, allele harmonisation, sensitivity diagnostics, and exportable results
Requirements:

R (≥ 4.2.0)
Key R dependencies: shiny, TwoSampleMR, ieugwasr, MVMR, DT, plotly, DiagrammeR, openxlsx, readxl (see Additional File 1 for full list)
OpenGWAS JWT token required for live API access (instrument extraction, LD clumping, outcome retrieval); a local instrument mode is available as a fallback
Archived version: https://doi.org/10.5281/zenodo.21441281<img width="468" height="28" alt="image" src="https://github.com/user-attachments/assets/bb7a9155-c4c8-4e4b-941d-34a6a01ac82e" />

License: MIT License
