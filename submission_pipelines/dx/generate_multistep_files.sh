#!/usr/bin/env bash
python3 munge_dx_arguments_msamp.py \
--files ./all_crams/ \
--out ./dx_jobs_2320_comparisons/ \
--n-each 40 \
--n-max 500000 \
--base-json ./inputs_batch_multi_new_liftover.json \
--make-name-col 'batch ID' \
--target-name-col 'sample_name' \
--job-base-name 'MitochondriaPipelineSwirl_Multi_'
