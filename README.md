# mtSwirl: Reference-aware quantification of mtDNA copy number and heteroplasmy using WGS

This repo contains pipeline files for the reference-aware mtSwirl pipeline as well as the code used to run, merge, and annotate the results.

## WDL

See the WDL folder for the self-contained WDL. The `v2.5_MongoSwirl_Single` folder contains the single-sample pipeline oriented for use with Cromwell. The `v2.6_MongoSwirl_Multi` folder contains a multi-sample pipeline for use on the UKB Research Analysis Platform using dxCompiler (https://github.com/dnanexus/dxCompiler). This folder also contains supporting scripts and reference NUMTs used to generate nucDNA self-reference sequences. See manuscript Methods for more details.

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

## Run GWAS in UKB

To run GWAS in UKB use the files in `gwas_ukb`. Using the outputs of QC, we run covariate correction with `generate_covariate_corrected_traits.Rmd` for mtCN (and for sensitivity analyses). To produce final heteroplasmy phenotypes, we use `produce_final_HL_traits.Rmd`. We use `saige_pan_ancestry_custom.py` to run SAIGE in UKB with `custom_load_custom_sumstats_into_mt.py` to combine results into an MT.

## Run GWAS in AoU

We use the files in `gwas_aou` to run GWAS in AoU. To produce custom PCs by recomputing them per-ancestry, we use `run_per_ancestry_pca.py`. We run `aou_run_full_hl_gwas.py` to run the GWAS.