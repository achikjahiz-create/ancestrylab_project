# Validation status

## Implemented in the dashboard

- Project-oriented Shiny navigation and dataset inspection.
- EIGENSTRAT and PLINK triplet validation.
- Exact population-label validation; missing labels cause a hard stop rather than proxy substitution.
- f2 cache with input hashes, allele-frequency-product mode, pseudohaploid adjustment and analysis settings in its cache signature; downstream cache freshness and exact group membership checks.
- ADMIXTOOLS 2 wrappers for qp3pop, qpDstat f4 and direct-genotype D-statistics, qpf4ratio, FST, qpWave, qpAdm, qpAdm rotation, qpGraph and FindGraph.
- qpAdm weights, fitted f4, rankdrop and popdrop diagnostics exposed in the interface.
- f2 configuration controls for autosomes, transition/transversion inclusion, MAF bounds, missingness, small-sample correction, pseudohaploid adjustment, allele-frequency products and memory limits.
- FST singleton-population guard with explicit exploratory override; jackknife/bootstrap options.
- qpWave/qpAdm/qpGraph bootstrap and regularization controls; FindGraph stopping and residual-search controls.
- Direct-genotype qpAdm path for `allsnps=TRUE`; cached f2 and direct-genotype paths are kept distinct because genotype-only options cannot be applied to an already-computed f2 array.
- Relate/GhostBuster input validation, per-chromosome and pickle-aggregation modes, selectable Conda Python executable, and separate LD-curve dating action.
- Per-run RDS/CSV/JSON persistence and ZIP archive export.
- Synthetic exact-label and design-validation smoke tests.

## Not yet validated

- The R/Shiny runtime and ADMIXTOOLS 2 calculations have not been executed in the current authoring environment because R/Rscript is unavailable there.
- The CI workflow must be run after these dashboard changes and must pass before merging/deploying.
- No full genotype-based integration test has yet compared f-statistics, qpWave or qpAdm results against a published benchmark dataset.
- GhostBuster and LD-curve dating require a separately installed upstream environment and valid Relate trees; no end-to-end GhostBuster run has been performed here.
- qpGraph/FindGraph results have not been benchmarked against a known graph and published residuals.
- No deployment, authentication, multi-user isolation, load test or security review has been performed.

## Required pre-publication validation

1. Run the CI tests and parse all R source files.
2. Run a small public ADMIXTOOLS 2 example dataset through each statistical module and compare with upstream examples.
3. Verify sample-to-population mapping, relatedness, SNP intersection, coverage and ascertainment using a pinned dataset release.
4. Reproduce at least one published f3/f4, qpWave and qpAdm result with the exact published cohort and parameter design.
5. Test qpAdm alternative right sets, source rotations, rankdrop/popdrop and block-resampling uncertainty.
6. Validate qpGraph residuals and FindGraph topology stability across seeds and search limits.
7. Run GhostBuster only with a valid Relate-tree example and matching genetic map/build.
8. Review storage, permissions, consent, retention, encryption and deployment security before using identifiable human genomic data.
