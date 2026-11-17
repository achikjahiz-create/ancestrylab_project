# Dataset intake and provenance

AncestryLab does not bundle genotype data. It accepts a local EIGENSTRAT or PLINK triplet and can retrieve public AADR metadata from the official Harvard Dataverse dataset.

## AADR

- Dataset DOI: https://doi.org/10.7910/DVN/FFIDCW
- Use the dashboard's metadata-only option first to inventory sample/group annotations.
- The notebook `AADR_Exact_Sample_Discovery_No_Proxies.ipynb` downloads `.ind` and `.anno` metadata and writes exact sample/group audits.
- A genotype panel can be multi-gigabyte. Confirm release, file ID, size and checksum before analysis.
- Keep AADR release versions fixed for reproducible comparisons. A new release can change group labels, sample counts and genotype content.

## Accepted local inputs

### EIGENSTRAT

A matching prefix must resolve to `.geno`, `.snp` and `.ind`. The dashboard reads sample population labels from the third `.ind` column and reports the number of individuals assigned to each label.

### PLINK

A matching prefix must resolve to `.bed`, `.bim` and `.fam`. The dashboard displays FID values as candidate population labels. In PLINK, FID is not inherently a population label; confirm the encoding in your project before running a population-level analysis.

## Exact cohorts and archaeological metadata

A population label is a technical identifier, not proof of a shared archaeological context. For ancient samples, verify individual sample IDs, site, period, calibrated dates, laboratory processing, contamination, coverage, relatedness and the publication's sample inclusion criteria. If a requested group is absent, keep it absent. Never replace it with a geographically or genetically similar population without a separately justified and explicitly labelled model.

## Required provenance

Archive the source URL/DOI, release/version, download date, file checksums, sample annotation version, genotype format, SNP intersection and all sample exclusions. The app's per-run manifest stores the selected genotype prefix and model settings; it does not copy or checksum multi-gigabyte source data for you.
