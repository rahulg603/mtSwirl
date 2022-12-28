#!/bin/bash

pip3 install gnomad
git clone https://github.com/rahulg603/gnomad-mitochondria.git
git clone https://github.com/broadinstitute/gnomad_qc.git
git clone https://github.com/broadinstitute/gnomad_methods.git
mv gnomad_qc gnomad_qc_hold
cd gnomad_qc_hold
mv gnomad_qc ../
cd ..
mv gnomad-mitochondria gnomad_mitochondria_hold
mv gnomad_mitochondria_hold/gnomad_mitochondria ./
cd gnomad_mitochondria/pipeline/

python dx_provision_sql.py --dx-init 220619_MitochondriaPipelineSwirl_v2_5_Multi_20k

python dx_collate_tables.py \
--pipeline-output-folder '220507_mitopipeline_v2_2_ukb_trial/' \
--vcf-merging-output 'tab_vcf_merging.tsv' \
--coverage-calling-output 'tab_coverage.tsv' \
--dx-init 220619_MitochondriaPipelineSwirl_v2_5_Multi_20k

dx upload tab_vcf_merging.tsv --path /220619_MitochondriaPipelineSwirl_v2_5_Multi_20k/merging/tab_vcf_merging.tsv
dx upload tab_coverage.tsv --path /220619_MitochondriaPipelineSwirl_v2_5_Multi_20k/merging/tab_coverage.tsv

python dx_annotate_coverage.py \
-i "tab_coverage.tsv" \
-o "coverage/testing_pipeline_trial_ukb_coverage.ht" \
-t "tmp/" \
--overwrite \
--dx-init 220507_mitopipeline_v2_2_ukb_trial

python dx_combine_vcfs.py \
-p "tab_vcf_merging.tsv" \
-c "coverage/testing_pipeline_trial_ukb_coverage.mt" \
-v vcf \
-a 'file:///mnt/project/reference/grch38_genome/blacklist_sites.hg38.chrM.bed' \
-a-ref GRCh38 \
-o "vcf/" \
-t "tmp/" \
-f 220506_testing_ukb_pipeline_variants \
--overwrite --include-extra-v2-fields \
--dx-init 220507_mitopipeline_v2_2_ukb_trial