# Nuclear genetic control of mtDNA copy number and heteroplasmy in humans
[![DOI](https://zenodo.org/badge/448948131.svg)](https://zenodo.org/badge/latestdoi/448948131)

This repo contains pipeline files for the reference-aware mtSwirl pipeline as well as the code used to run, merge, and annotate the results.

## Citation and data

This pipeline was released as part of the manuscript: `Nuclear genetic control of mitochondrial DNA copy number and heteroplasmy in humans`, which can be found at [Nature](https://www.nature.com/articles/s41586-023-06426-5). If you use these resources in your work, please cite as `Gupta et al. 2023 Nature`:

```
Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA copy number and heteroplasmy in humans. Nature, in press. https://doi.org/10.1038/s41586-023-06426-5.
```

### Individual-level data

Individual level data corresponding to mtDNA copy number (before and after covariate correction) and the post-QC variant callset can be found:

- For UKB, via the [UKB data showcase](https://biobank.ndph.ox.ac.uk/ukb/). Note that final data return is currently in process.
- For AoU, as part of the `Nuclear genetic control of mtDNA copy number and heteroplasmy in humans` [workspace](https://workbench.researchallofus.org/workspaces/aou-rw-3273c7f0/nucleargeneticcontrolofmtdnacopynumberandheteroplasmyinhumans/data). Note that controlled tier access is required to clone this workspace.

### Summary statistics

Summary statistics from UKB are available:

- Via GWAS Catalog under ID [GCP000614](https://www.ebi.ac.uk/gwas/publications/37587338), where we have uploaded summary statistics corresponding to our largest analysis for each phenotype, corresponding to cross-ancestry meta-analyses when performed or EUR when no other populations had sufficient N for GWAS. These summary statistics are filtered to include only stringently "high_quality" variants; the full summary statistics including all otherwise QC-passing variants can be found on GCP (see below). **PLEASE NOTE: these data were corrected on 03/2024 as the `effect_allele` and `other_allele` columns were originally reversed. No other columns were changed. No data deposited in other locations (e.g., GCP, AllofUs; see below) required updating.**
- On Google Cloud Platform, in the `gs://mito-wgs-public-2023` bucket. Please note that this is a requester pays bucket. This bucket also contains `ukb_b37_b38_lifted_variants.tsv.bgz`, which maps GRCh37 coordinates in the UKB data to GRCh38. The summary statistics on GCP correspond to the same data, but are stored using the [Pan UKB schema](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files#per-phenotype-files). These files contain the cross-ancestry meta-analysis as well as per-ancestry association statistics as well (and thus are more comprehensive than those on GWAS Catalog). More information on the schema is described in the `README_ukb.md` file located in the `gs://mito-wgs-public-2023` bucket.

Summary statistics from AoU are available in the `Nuclear genetic control of mtDNA copy number and heteroplasmy in humans` [workspace](https://workbench.researchallofus.org/workspaces/aou-rw-3273c7f0/nucleargeneticcontrolofmtdnacopynumberandheteroplasmyinhumans/data) in the same format as UKB summary statistics found on GCP. Note that controlled tier access is required to clone this workspace.

See Supplementary table 1 for sample size information.

### AllofUs workspace access

Please note that at the time of writing, there is no mechanism by which custom workspaces in AoU can be made available to anyone with controlled tier access. Thus, we ask that in the interim, any users who desire to work with these data in AoU contact us to be added to the workspace. We are committed to making these data automatically available when this mechanism becomes available, and plan to beta-test this functionality when it is possible to do so.

## mtSwirl: Reference-aware quantification of mtDNA copy number and heteroplasmy using WGS

See the WDL folder for the self-contained WDL. The `v2.5_MongoSwirl_Single` folder contains the single-sample pipeline oriented for use with Cromwell. The `v2.6_MongoSwirl_Multi` folder contains a multi-sample pipeline for use on the UKB Research Analysis Platform using [dxCompiler](https://github.com/dnanexus/dxCompiler). This folder also contains supporting scripts and reference NUMTs used to generate nucDNA self-reference sequences. See manuscript Methods for more details.

## Generate multi-sample MatrixTables and perform QC

The `generate_mtdna_call_mt` folder contains code used to merge single-sample VCFs into Hail MatrixTables. This code was written originally as an extension of code previously released for mtDNA analysis (Laricchia et al. 2022 Genome Res). Scripts in the root of this folder work across any platform; scripts in each sub-folder are platform specific.

### dx

Run `dx_pipeline.sh` to run the merging pipeline.

### AoU

1. Currently, AoU does not have a central Cromwell implementation. Thus, we created `aou_mtdna_analysis_launcher.sh` to run the WDL. Tweak the parameters in the header for your configuration.
2. To combine per-base coverage into an MT use `aou_annotate_coverage.py`
3. To combine single-sample VCFs into an MT use `aou_combine_vcfs.py`

### Terra

1. To combine per-base coverage into an MT use `annotate_coverage.py`
2. To combine single-sample VCFs into an MT use `combine_vcfs.py`

### All platforms

1. To generate sample statistics after QC (e.g., mtCN), use `process_sample_stats.py`
2. To annotate the VCF MatrixTable, run QC, run VEP, and output a QC'd variant flat file, use `add_annotations.py`

## Genome-wide association study pipeline

### UKB

To run GWAS in UKB use the files in `gwas_ukb`. Using the outputs of QC, we run covariate correction with `generate_covariate_corrected_traits.Rmd` for mtCN (and for sensitivity analyses). To produce final heteroplasmy phenotypes, we use `produce_final_HL_traits.Rmd`. We use `saige_pan_ancestry_custom.py` to run SAIGE in UKB with `custom_load_custom_sumstats_into_mt.py` to combine results into an MT.

### AllofUs

We use the files in `gwas_aou` to run GWAS in AoU. To produce custom PCs by recomputing them per-ancestry, we use `run_per_ancestry_pca.py`. We run `aou_run_full_hl_gwas.py` to run the GWAS.
