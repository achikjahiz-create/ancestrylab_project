# AncestryLab

**An academic ancient-DNA population-genetics workspace** for explicit sample selection, reproducible f-statistics, qpWave/qpAdm modelling, admixture-graph exploration and Relate/GhostBuster workflows.

> **Status: research prototype, not a validated commercial or clinical product.** The dashboard and engine are orchestration software. Numerical results are produced by upstream methods and must be independently checked against their documentation and the relevant study supplements. The app does not certify an ancestry interpretation.

## What is in the dashboard

| Workspace | Capabilities | Key controls |
|---|---|---|
| Project & data | EIGENSTRAT / PLINK intake, AADR metadata download, exact group inventory | Genotype prefix, panel, metadata-only retrieval, population labels |
| f-statistics | f2 block cache, qp3pop f3, qpDstat f4 / D, qpf4ratio | Exact population tuples, block/MAF/transition filters, missingness, pseudohaploid and small-sample correction, cache invalidation |
| FST | Pairwise differentiation | Population sets and jackknife/bootstrap uncertainty |
| qpWave | Rank testing | Left/right sets, covariance fudge and bootstrap |
| qpAdm | Mixture estimation | Target, `left` sources, `right` outgroups, constrained weights, fudge-twice, fitted-f4 diagnostics, jackknife/bootstrap, cached-f2 or direct-genotype/allsnps path |
| qpAdm rotation | Alternative left/right model sensitivity | `leftright` candidate set, fixed right set, full results |
| qpGraph | Fit a specified ADMIXTOOLS 2 graph | Graph object in RDS format, random starts/seed, covariance regularization, objective, constraints, residuals and bootstrap |
| FindGraph | Candidate topology search | Root/outgroup, admixture bounds, stopping score, residual optimisation, z-threshold, search generations and time limit |
| GhostBuster | Local ancestry / genealogical analysis | Relate trees, population labels, matching recombination map, genome build, sample IDs, chromosomes and EM settings |
| LD curve dating | One-pulse, two-pulse and continuous models | GhostBuster output prefix and genome build |
| Run archive | Persisted result objects, CSV tables and JSON manifests | Per-run provenance and downloadable ZIP archive |

## Scientific guardrails

- **No proxy substitution.** Population labels must match the input metadata exactly after trimming surrounding whitespace. A missing group is an error, not an invitation to substitute a nearby or genetically similar group.
- **Population labels are not archaeological validation.** Before modelling ancient groups, verify sample IDs, site, date, archaeological context, contamination, coverage, relatedness, SNP overlap and the primary publication/supplement.
- **f2 is the reusable computational base.** The cache is keyed to input file hashes, selected population set, block size, missingness threshold, allele-frequency-product setting, pseudohaploid adjustment, small-sample correction, autosome/transition filters, MAF bounds, memory limit and ADMIXTOOLS 2 version. Changing one of these creates a different cache. Downstream actions are blocked if the current genotype prefix or shared f2 settings differ from those used to build the in-memory cache.
- **f4, D and f4-ratio are not interchangeable.** D-statistics are computed directly from genotype files (`f4mode=FALSE`) and do not inherit the shared f2 cache's pseudohaploid/small-sample corrections; f4 uses shared f2 blocks; f4-ratio uses the five-population formula shown in the dashboard. `allsnps` mode is optional and explicitly flagged as exploratory because it can inflate false positives under non-random missingness. `afprod` can improve f4 precision with substantial missingness and is permitted for outgroup f3 only when the same fixed outgroup is used in every tuple.
- **qpAdm is conditional.** Inspect weights and standard errors, fitted f4 diagnostics, rankdrop, popdrop, model feasibility, right-set robustness, SNP coverage and alternative source combinations. A high p-value alone does not establish a historically true model.
- **qpAdm rotation is not automatic model discovery.** The interface uses ADMIXTOOLS 2's `leftright` / `rightfix` design and requests full results so candidate models can be examined rather than ranked by one p-value alone.
- **FST and pseudohaploid data.** The Hudson FST estimator uses block-jackknife standard errors by default; block bootstrap is an optional sensitivity procedure. The dashboard exposes ADMIXTOOLS 2 pseudohaploid genotype adjustment, which should normally remain enabled for pseudo-haploid calls. Single-individual populations remain unsuitable for unbiased FST.
- **GhostBuster is not an f-statistics tool.** Upstream currently requires Relate-inferred trees. An AADR `.geno/.snp/.ind` panel by itself is not a valid GhostBuster input. Use a recombination map that matches the genome build and preserve the Relate/sample mapping.
- **Graph searches are not unique histories.** qpGraph/FindGraph solutions depend on the populations, constraints, search limits and graph assumptions supplied.
- **Missingness, SNP ascertainment and block design matter.** Record and justify the choices used for each analysis. The UI records the main controls but does not replace a study-specific protocol.

## Install and launch locally

Requirements: R 4.4 or later, Python 3, compiler/system libraries required by ADMIXTOOLS 2, and a working internet connection for package installation and AADR downloads.

```bash
git clone https://github.com/achikjahiz-create/ancestrylab.git
cd ancestrylab
Rscript scripts/install.R
python -m pip install -r requirements.txt
R -e "shiny::runApp('app', host='0.0.0.0', port=3838)"
```

Open `http://localhost:3838`.

Docker template:

```bash
docker compose up --build
```

The Docker image installs the R/Shiny analysis layer. **GhostBuster has its own Python 3.13/Conda/GSL requirements** and should be installed separately; see [docs/GHOSTBUSTER.md](docs/GHOSTBUSTER.md).

## AADR sample discovery — exact labels only

The companion notebook `notebooks/AADR_Exact_Sample_Discovery_No_Proxies.ipynb` retrieves AADR `.ind` and `.anno` metadata, inventories group labels, applies declared exact-label/site/entity/date/period rules and exports a complete selection audit. Missing requested cohorts remain missing. It is a sample-discovery step, not an ancestry model and not a reproduction of a specific paper's qpAdm analysis.

## Reproducing a published study

For *The Genomic History of the Bronze Age Southern Levant* (Agranat-Tamir et al., 2020), use the actual target and source samples, SNP intersection, outgroup set, qpWave/qpAdm design, block-jackknife or bootstrap procedure, rank diagnostics and sensitivity analyses in the paper and supplement. Do not substitute a generic least-squares mixture fit or silently change the cohort definitions.

See [docs/ACADEMIC_WORKFLOW.md](docs/ACADEMIC_WORKFLOW.md), [docs/DATASETS.md](docs/DATASETS.md), [docs/VALIDATION_STATUS.md](docs/VALIDATION_STATUS.md) and [docs/GHOSTBUSTER.md](docs/GHOSTBUSTER.md).

## Reproducibility and outputs

Each successful analysis writes to `results/runs/<run_id>/`:

- `result.rds` — full native R result object;
- `result.csv` — tabular view where the result can be represented as a data frame;
- `manifest.json` — UTC timestamp, method, input dataset, parameters, R version and ADMIXTOOLS 2 version.

The f2 cache is stored separately under `results/f2_cache/`. The Run archive tab packages the run folders into a ZIP. Graph input RDS files are copied into their run folders. GhostBuster output files remain at the selected output prefix; the run manifest records their paths, sizes and checksums rather than duplicating potentially large files. Keep a copy of the original genotype files and metadata with your study records; the manifest records their checksums but does not duplicate multi-gigabyte input data.

## Security and deployment

The current app is intended for local or controlled research use. It has **no user authentication, multi-tenant isolation, role-based permissions, encrypted object storage, job queue, resource quotas, audit-log signing or automated deletion policy**. Do not expose it publicly or upload sensitive human genomic data to a shared instance without implementing institutional security, consent, governance and data-protection controls. A production service requires a separate security and operations review.

## References

- Agranat-Tamir, L. et al. (2020). *The Genomic History of the Bronze Age Southern Levant*. Cell. https://doi.org/10.1016/j.cell.2020.04.024
- Allen Ancient DNA Resource (AADR). Harvard Dataverse. https://doi.org/10.7910/DVN/FFIDCW
- ADMIXTOOLS 2 documentation: https://uqrmaie1.github.io/admixtools/
- ADMIXTOOLS 2 source: https://github.com/uqrmaie1/admixtools
- GhostBuster (Myers Group): https://github.com/MyersGroup/GhostBuster
- Relate: https://myersgroup.github.io/relate/

AncestryLab orchestration code is MIT-licensed. ADMIXTOOLS 2, GhostBuster, datasets and other upstream components retain their own licences and citation requirements. No genotype dataset is bundled.
