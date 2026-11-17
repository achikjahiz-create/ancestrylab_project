# Academic workflow and interpretation protocol

## 1. Freeze the dataset before modelling

1. Record the source archive, release/version, download date and checksum for every genotype and metadata file.
2. Confirm the genotype format and inspect the `.ind`/`.fam` sample-to-population mapping. For EIGENSTRAT, the app reads population labels from column 3 of `.ind`. For PLINK, FID is only displayed as a candidate label; verify that FID really encodes population membership.
3. Confirm sample IDs and group assignments against the original publication and supplement. Do not rely on a group name alone to establish site, date or archaeological context.
4. Apply the paper-specific SNP intersection and filtering protocol. Preserve the exact marker list and record the number of overlapping SNPs for every population/model.
5. Keep all missing cohorts missing. No nearest-label, nearest-date, geographic or genetic proxy is inserted automatically.

## 2. Build the f2 cache

Select only the populations needed for the planned tests. Record:

- exact genotype prefix and file hashes;
- population labels included in the cache;
- block size in Morgans;
- missingness threshold;
- autosomal/chromosome and SNP filters applied upstream;
- ADMIXTOOLS 2 version.

The cache signature is generated from the input hashes, population set, block size, missingness, allele-frequency-product mode, pseudohaploid adjustment, SNP filters and package version. If any input changes, rebuild the cache. `afprod=TRUE` returns negative average allele-frequency products instead of f2; use it for f4/qpAdm with substantial missingness, or for outgroup f3 only when the same fixed outgroup is used in every tuple.

## 3. f3, f4 and D-statistics

The UI accepts one tuple per line:

- f3: `target,source1,source2`;
- f4 or D: `A,B,C,D`;
- f4 ratio: `A,B,C,D,E`, computing f4(A,B;C,D) / f4(A,B;E,D).

f4 is calculated from the shared f2 blocks. D-statistics are calculated directly from the genotype prefix with `f4mode=FALSE`; they cannot be reconstructed from the cached f2 array through this interface. The direct `f4blockdat_from_geno` path has a narrower filter set than `f2_from_geno` and does not inherit the shared f2 pseudohaploid/small-sample corrections. The `allsnps` option can increase power with missing data but can also create false positives; keep it off unless justified. Preserve tuple order. State whether an f3 is an outgroup or admixture f3 and define the f4/D sign convention in the study report. Use the relevant ADMIXTOOLS 2 options for pseudohaploid samples and missing data; do not compare results from incompatible SNP sets without explanation.

## 4. FST

The app computes the Hudson FST estimator directly from the genotype prefix. Block-jackknife standard errors are the default; the dashboard can switch to block bootstrap with a specified number of replicates. The pseudohaploid genotype adjustment is exposed and should normally remain enabled when such calls are present. ADMIXTOOLS 2 documentation notes that at least two independent observations per population are needed for an unbiased estimate. Pseudohaploid singletons need special handling; do not treat a one-individual group as a reliable FST population. Report sample counts and the exact resampling procedure.

## 5. qpWave and qpAdm

For qpWave, specify the left populations and the right/reference populations. Interpret the rank table as a set of nested rank tests, not a direct count of named ancestral populations.

For qpAdm, define:

- one target population;
- a left/source set;
- a right/outgroup set;
- whether coefficients are constrained to be non-negative;
- whether the doubled diagonal `fudge_twice` option is used (ADMIXTOOLS 2 documents this as a way to better match original qpAdm p-values);
- whether fitted f4 diagnostics are returned;
- block-jackknife or block-bootstrap uncertainty and the number of bootstrap replicates;
- whether f2 blocks are reused or f4 blocks are computed directly from the genotype prefix with `allsnps=TRUE`. The direct-genotype path does not inherit the f2 cache's pseudohaploid correction and MAF/transition filters; document this trade-off when using pseudo-haploid ancient samples.

The dashboard displays the weights, `rankdrop` and `popdrop` results. Inspect standard errors, rank, model fit, nested tests, feasibility and sensitivity to alternative right sets and source combinations. Record negative weights and failed models rather than suppressing them. Do not call the largest p-value the “best historical model” without additional justification.

### qpAdm rotation

ADMIXTOOLS 2 `qpadm_rotate()` distributes the `leftright` candidates between the model's left and right sets while keeping `target` and `rightfix` fixed. This dashboard exposes those exact concepts. Use rotation to assess design sensitivity; it does not establish that the candidate set is complete or historically correct.

## 6. qpGraph and FindGraph

qpGraph accepts a graph object serialized in RDS format using the ADMIXTOOLS 2 graph representation. Save the graph source file alongside the run. Examine residuals, worst-fitting f-statistics, graph constraints and alternative topologies.

FindGraph explores candidate topologies under the supplied outgroup, admixture count, maximum admixture count, stopping score, worst-residual option, z-threshold and search limits. Search output is heuristic and can depend on initialisation and time limits. Re-run with alternative settings and compare graph fit and residual structure. Do not interpret one returned topology as uniquely identified history.

## 7. GhostBuster and LD-curve dating

GhostBuster is a separate genealogical workflow. It currently requires Relate-inferred trees, population labels, recombination maps and sample IDs. It is not run on AADR genotype panels directly. Confirm that the map, tree coordinates and genome build are consistent. Record Relate version and inference settings, sample/haplotype mapping, chromosome coverage, GhostBuster commit, EM restarts, iteration count, random seed and HMM setting.

LD-curve dating is downstream of GhostBuster outputs. The dashboard offers one-pulse, two-pulse and continuous-admixture model runs. Report model choice, fit diagnostics, generation-time assumption and uncertainty; do not convert a fitted time into a calendar date without stating the years-per-generation assumption.

## 8. Reproducibility package

Archive the full genotype metadata and checksums, sample-selection audit, SNP intersection, exact population lists, configuration, software versions, raw outputs and run manifests. The dashboard stores result objects and manifests, but the researcher remains responsible for archiving large source data and the exact study supplement.
