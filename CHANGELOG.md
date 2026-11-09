# Changelog

Notable changes to this pipeline, newest first. Versions follow
[Semantic Versioning](https://semver.org/); each release is a git tag.
`git log` remains the complete record of what changed in each commit.

Each entry states whether outputs from the previous version must be regenerated.

Entry format:

```
## [x.y.z] - YYYY-MM-DD
### <module or file>
- Change: <what changed>
- Reason: <why; cite a source for methodological changes>
- Rerun required: <yes/no; which modules>
```

---

## [Unreleased]

## [0.1.0] - 2026-11-09

Initial public release.

- Data acquisition from the GDC (TCGA-HNSC), restricted to oral cavity SCC by
  tissue/organ of origin; each layer is downloaded only if its processed output
  is absent.
- Patient- and sample-level metadata. Primary tumor (`01`) and normal adjacent
  tissue (`11`) samples only. Mutation flags and burden from non-synonymous
  variants in primary tumors; tumor-intrinsic attributes are `NA` for NAT.
- limma-voom differential expression (mRNA: Tumor vs Normal, Late vs Early,
  Stage IV vs Stage I; miRNA: Tumor vs Normal), with within-patient correlation
  (`duplicateCorrelation`) for Tumor vs Normal. Genes keyed by Ensembl ID with
  GDC-supplied annotation; analysed gene types set in the config.
- Precursor-to-mature miRNA mapping and miRTarBase target pairs via multiMiR,
  restricted by support category.
- WGCNA split into a diagnostic stage (3A) and a network stage (3B). The
  soft-thresholding power is recorded by a person in `config/runs/<run_id>.json`
  with an input fingerprint; 3B refuses to run without a matching record.
- Signed network built in a single block; module-trait p-values use the
  per-pair number of non-missing samples; enrichment background is the network
  gene set.
- Configuration in `config/config.json`; timestamped logs with git commit hash.
