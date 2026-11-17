# GhostBuster and Relate workflow

GhostBuster is a separate genealogy-based local ancestry workflow, not an f-statistics method. Upstream currently supports **Relate-inferred trees**. An AADR `.geno/.snp/.ind` triplet alone is not valid GhostBuster input.

Upstream repository and installation guide: https://github.com/MyersGroup/GhostBuster

## Installation

```bash
conda create -n ghostbuster python=3.13 -y
conda activate ghostbuster
conda install -c conda-forge gsl
python -m pip install --upgrade pip setuptools wheel
# from a cloned GhostBuster checkout
python -m pip install -r requirements.txt
```

Keep this environment separate from the R/Shiny application. The current Docker image does not install the full GhostBuster environment.

## Inputs required in the dashboard

- GhostBuster checkout directory containing `ghost_buster.py` and `ld_curve_dating.py`;
- Relate tree prefix (for chromosome `n`, the corresponding file must be `<prefix><n>.trees`);
- population-label file;
- recombination map compatible with the tree coordinates;
- genome build (`hg37` or `hg38`);
- haplotype/sample IDs and chromosome list.

The dashboard exposes the Python executable from the GhostBuster Conda environment, EM cluster count, random restarts, iteration count, seed, HMM setting, intermediate-pickle-only mode and a separate aggregation mode. In aggregation mode it omits `--trees` and expects the chromosome pickle files to exist already. It validates the selected tree files, population-label file, recombination-map path/prefix and upstream scripts before launching. It cannot validate that your tree/sample mapping is biologically correct; verify this independently.

Example upstream invocation:

```bash
python ghost_buster.py \
  --genome_build hg37 \
  --trees example/relate_chr \
  --poplabels example/poplabels.txt \
  --rec example/genetic_map_GRCh37_chr \
  --sample_id 0 1 2 3 \
  --chrs 22 \
  --output results/ghostbuster/relate \
  --num_clusters 2 --n_repeats 20 --num_iters 200 --seed 2
```

For chromosome-parallel processing, upstream supports `--only_make_pickle_files`; after generating per-chromosome pickle files, use the dashboard aggregation mode with the same output prefix and population/sample mapping. Aggregation omits `--trees`. Consult the current upstream README for exact cluster-specific execution.

## LD-curve dating

After a successful GhostBuster run, the dashboard can launch `ld_curve_dating.py` in one-pulse (`1date`), two-pulse (`2date`) or continuous-admixture mode. Selecting “Run all models” omits `--mode`, matching the upstream default. Report the chosen model, goodness-of-fit diagnostics, generation-time assumption and uncertainty. A model-derived generation count is not itself a calendar date.

## Output provenance

Archive the Relate version and tree-inference settings, input sample/haplotype mapping, recombination-map source/build, chromosome coverage, GhostBuster commit, EM settings, random seed and all output files. GhostBuster outputs include posterior local-ancestry membership by window, coalescence-rate tables and EM parameter files. Do not automatically equate inferred genealogical components with named archaeological populations.

GhostBuster is GPL-3.0; it is not bundled or copied into this repository. Cite the upstream project and follow its licence.
