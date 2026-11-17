# Changelog

## Academic dashboard rebuild — 2026-11-17

- Replaced the single-purpose sample-discovery interface with a multi-workspace Shiny dashboard.
- Added EIGENSTRAT/PLINK inspection, exact population-label selection, and AADR metadata retrieval.
- Added reproducible f2 block caching with input hashes, exact group membership checks, MAF/autosome/transition filters, missingness, pseudohaploid adjustment, small-sample correction, allele-frequency-product mode and memory controls.
- Added f3, f4, direct-genotype D-statistics, f4-ratio and FST workflows.
- Added qpWave, qpAdm, qpAdm rotation, qpGraph and FindGraph with model, regularization, bootstrap, search and diagnostic controls.
- Added Relate/GhostBuster execution, per-chromosome pickle generation/aggregation and LD-curve dating controls.
- Added per-run RDS/CSV/JSON provenance, graph-input archiving, external GhostBuster output checksums and ZIP export.
- Added Docker, Compose, Render template, installation script, CI workflow, documentation and smoke tests.

### Validation boundary

The authoring environment did not contain R/Rscript. Static delimiter checks, Python compilation, notebook JSON validation, YAML parsing, shell syntax checks and a mocked GhostBuster launcher test were run. ADMIXTOOLS 2 numerical integration, Shiny runtime, Docker build, GhostBuster and deployment remain unvalidated; see `docs/VALIDATION_STATUS.md`.
