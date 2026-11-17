# Reproducible analysis workflow

The preferred interface is the Shiny dashboard (`app/`). The engine wrappers in `R/engine.R` can also be called directly.

## 1. Inspect samples and groups

```bash
python scripts/download_aadr.py --panel 1240K --out data/aadr --metadata-only
```

The downloader records the Dataverse version, file IDs, sizes and checksums. Use the dashboard to inspect a local EIGENSTRAT or PLINK triplet and select exact population labels. For PLINK, confirm that FID truly encodes populations.

## 2. Build f2 blocks

```r
source("R/engine.R")
prefix <- "data/aadr/<verified-genotype-prefix>"
blocks <- compute_f2(prefix, "results/f2_cache", pops=c("TARGET", "SOURCE_A", "SOURCE_B", "OUTGROUP"),
                     blgsize=0.05, maxmiss=0, afprod=FALSE,
                     adjust_pseudohaploid=TRUE)
```

`f2_from_geno()` returns a 3D array of blockwise f2 statistics. The cache signature includes input hashes, selected populations, block size, missingness, `afprod`, pseudohaploid adjustment and ADMIXTOOLS 2 version. Do not reuse a cache after changing any of those settings.

## 3. Standard statistics

```r
admixtools::qp3pop(blocks, pop1="TARGET", pop2="SOURCE_A", pop3="SOURCE_B")
admixtools::qpdstat(blocks, pop1="A", pop2="B", pop3="C", pop4="D", f4mode=TRUE)
admixtools::qpdstat(prefix, pop1="A", pop2="B", pop3="C", pop4="D", f4mode=FALSE)
admixtools::qpf4ratio(blocks, pops=matrix(c("A","B","C","D","E"), nrow=1))
admixtools::fst(prefix, pop1="POP_A", pop2="POP_B", boot=FALSE)
```

D-statistics require direct genotype input in this workflow. F4-ratio uses the ordered five-population formula shown in the dashboard. FST uses block-jackknife uncertainty by default and can use block bootstrap.

## 4. qpWave, qpAdm and rotation

```r
admixtools::qpwave(blocks, left=c("A","B","C"), right=c("R1","R2","R3"))
admixtools::qpadm(blocks, left=c("SOURCE_A","SOURCE_B"), right=c("R1","R2","R3"), target="TARGET")
admixtools::qpadm_rotate(blocks, leftright=c("L1","L2","L3","L4"), target="TARGET", rightfix=c("R1","R2"), full_results=TRUE)
```

Inspect qpAdm weights, standard errors, rankdrop and popdrop. Preserve infeasible and failed models. Rotation is a sensitivity analysis, not an automatic discovery of historical ancestry.

## 5. Graphs and GhostBuster

Supply an ADMIXTOOLS 2 graph object in RDS format to qpGraph. FindGraph is a heuristic topology search; retain all candidate outputs and search settings.

GhostBuster is a separate Relate-tree workflow. AADR genotype files alone are not valid GhostBuster input. Install GhostBuster separately and supply matching Relate trees, population labels, sample IDs, recombination map and genome build. The dashboard also exposes downstream LD-curve dating models.

## 6. Archive

The dashboard saves native RDS results, CSV tables and JSON manifests under `results/runs/<run_id>/`. Retain original data, file checksums, SNP filters, sample selection audit, software versions and the full paper/supplement model specification.
