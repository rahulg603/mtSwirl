#!/usr/bin/env python3

import sys
import copy
import argparse
import logging
from datetime import date
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename='saige_pipeline.log')

from gnomad.utils import slack
from ukb_common import *
import time
import re
import hailtop.batch as hb

from ukbb_pan_ancestry import *
from ukb_common.utils.saige_pipeline import *

logger = logging.getLogger("saige_pan_ancestry_custom")
logger.addHandler(logging.StreamHandler(sys.stderr))
bucket_vcf = 'gs://ukb-diverse-pops'
root_vcf = f'{bucket_vcf}/results'
bucket = 'gs://rgupta-assoc'
root = f'{bucket}/saige_gwas'
pheno_folder = f'{bucket}/phenotype'
temp_bucket = 'gs://ukbb-diverse-temp-30day-multiregion'

DESCRIPTION_PATH = f'{pheno_folder}/field_descriptions.tsv'

HAIL_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/hail_utils:6.1'
#PHENO_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/rgupta_hail_utils'
#SAIGE_DOCKER_IMAGE = 'gcr.io/ukbb-diversepops-neale/saige:1.1.5'
PHENO_DOCKER_IMAGE = 'us-docker.pkg.dev/ukbb-diversepops-neale/pan-ukbb-docker-repo/rgupta-hail-utils'
SAIGE_DOCKER_IMAGE = 'us-docker.pkg.dev/ukbb-diversepops-neale/pan-ukbb-docker-repo/saige:1.1.5'
QQ_DOCKER_IMAGE = 'konradjk/saige_qq:0.2'

SCRIPT_DIR = '/ukb_common/saige'
# TODO: add binary_trait annotation to input table and remove this:
saige_pheno_types = {
    'continuous': 'quantitative',
    'biomarkers': 'quantitative',
    'categorical': 'binary',
    'icd': 'binary',
    'icd_first_occurrence': 'binary',
    'icd_all': 'binary',
    'phecode': 'binary',
    'prescriptions': 'binary'
}


def get_custom_ukb_pheno_mt_path(suffix):
    return f'{pheno_folder}/mt/phenotype_{suffix}.mt'


def get_custom_phenotype_summary_backup_path(suffix, curdate):
    return f'{pheno_folder}/summary/all_pheno_summary_{suffix}_before_{curdate}.txt.bgz'


def get_custom_phenotype_summary_path(suffix, extension = 'ht'):
    return f'{pheno_folder}/summary/phenotype_{suffix}.{extension}'


def get_custom_munged_pheno_path(suffix):
    return f'{pheno_folder}/mt/munged/munged_raw_phenotype_{suffix}.mt'


def get_custom_ukb_pheno_mt(suffix, pop: str = 'all', custom_covars: str = None):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(suffix))
    # mt = mt.annotate_rows(**get_ukb_meta(key_type=hl.tint32)[mt.row_key])
    mt = mt.annotate_rows(**get_covariates_with_custom(key_type=hl.int32, custom=custom_covars)[mt.row_key])
    if pop != 'all':
        mt = mt.filter_rows(mt.pop == pop)
    return mt


def get_covariates_with_custom(key_type, custom):
    new_covariates = get_covariates(key_type=key_type)
    if custom is None:
        return new_covariates
    else:
        custom_covariates = hl.import_table(custom, impute=True)
        custom_covariates = custom_covariates.key_by(s = key_type(custom_covariates.s))
        return new_covariates.annotate(**custom_covariates[new_covariates.key])#.checkpoint(f'{temp_bucket}/temp_covars_2.ht', overwrite=True)


def custom_get_phenos_to_run(suffix, pop, limit, specific_phenos, single_sex_only,
                             skip_case_count_filter, sex_stratified):
    ht = hl.read_table(get_custom_phenotype_summary_path(suffix))
    ht = ht.filter(ht.pop == pop)
    min_cases = MIN_CASES_EUR if pop == 'EUR' else MIN_CASES

    criteria = True
    if not skip_case_count_filter:
        criteria &= (ht.n_cases_by_pop >= min_cases)

    if single_sex_only:
        prop_female = ht.n_cases_females / (ht.n_cases_males + ht.n_cases_females)
        criteria &= ((prop_female <= 0.1) | (prop_female >= 0.9))

    ht = ht.filter(criteria).key_by()

    if sex_stratified:
        ht_sex_specific = ht.annotate(pheno_sex='males').union(ht.annotate(pheno_sex='females'))
        if sex_stratified == 'all':
            ht = ht.union(ht_sex_specific)
        else:
            ht = ht_sex_specific

    output = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select(*PHENO_KEY_FIELDS).collect()])
    if specific_phenos:
        specific_phenos = specific_phenos.split(',')
        output = [x for x in output if all(map(lambda y: y is not None, x)) and any([re.match(pcd, '-'.join(x)) for pcd in specific_phenos])]
    if limit:
        output = set(sorted(output)[:limit])

    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in output]
    return pheno_key_dict


def custom_summarize_data(suffix, overwrite):
    mt = hl.read_matrix_table(get_custom_ukb_pheno_mt_path(suffix))
    ht = mt.group_rows_by('pop').aggregate(
        stats=hl.agg.stats(mt.both_sexes),
        n_cases_by_pop=hl.if_else(hl.set({'continuous', 'biomarkers'}).contains(mt.trait_type),
                                  hl.agg.count_where(hl.is_defined(mt.both_sexes)),
                                  hl.int64(hl.agg.sum(mt.both_sexes)))
    ).entries()
    ht = ht.key_by('pop', *PHENO_KEY_FIELDS)
    ht = ht.checkpoint(get_custom_phenotype_summary_path(suffix), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.flatten().export(get_custom_phenotype_summary_path(suffix, 'tsv'))


def custom_add_description(mt):
    ht = hl.import_table(DESCRIPTION_PATH).key_by('phenotype_id')
    mt = mt.annotate_cols(description = ht[mt.phesant_pheno].description)
    mt = mt.annotate_cols(description = hl.if_else(hl.is_missing(mt.description), "", mt.description))
    return mt


def custom_load_custom_pheno(data_path, trait_type, modifier, source, sex: str = 'both_sexes', extension: str = 'txt', sample_col='s'):
    print(f'Loading {data_path}...')
    if extension == 'ht':
        ht = hl.read_table(data_path)
    else:
        if extension == 'tsv.gz':
            ht = hl.import_table(data_path, impute=True, force=True)
        else:
            ht = hl.import_table(data_path, impute=True)
        ht = ht.key_by(userId=ht[sample_col])
        if sample_col != 'userId':
            ht = ht.drop(sample_col)
        if trait_type == 'categorical':
            ht = ht.annotate(**{x: hl.bool(ht[x]) for x in list(ht.row_value)})

    mt = pheno_ht_to_mt(ht, trait_type, rekey=False).annotate_cols(data_type=trait_type)
    mt = custom_add_description(mt)
    mt = mt.key_cols_by(trait_type=trait_type, phenocode=mt.phesant_pheno, pheno_sex=sex, coding=NULL_STR_KEY,
                        modifier=modifier).drop('phesant_pheno')
    mt = mt.annotate_cols(category=source)
    return mt


def produce_custom_phenotype_mt(data_path, extn, suffix, trait_type, modifier, source, sample_col='s', append=True, overwrite=False, custom_covars=None):
    curdate = date.today().strftime("%y%m%d")
    mt = custom_load_custom_pheno(data_path, trait_type=trait_type, modifier=modifier, 
                                  source=source, sample_col=sample_col,
                                  extension=extn
                                  ).checkpoint(get_custom_munged_pheno_path(suffix), overwrite=overwrite)
    cov_ht = get_covariates_with_custom(hl.int32, custom_covars).persist()
    mt = combine_pheno_files_multi_sex_legacy({'custom': mt}, cov_ht)

    mt.group_rows_by('pop').aggregate(
        n_cases=hl.agg.count_where(mt.both_sexes == 1.0),
        n_controls=hl.agg.count_where(mt.both_sexes == 0.0),
        n_defined=hl.agg.count_where(hl.is_defined(mt.both_sexes))
    ).entries().drop(*[x for x in PHENO_COLUMN_FIELDS if x != 'description']).show(100, width=180)
    
    mt_path = get_custom_ukb_pheno_mt_path(suffix)
    if append and hl.hadoop_exists(f'{mt_path}/_SUCCESS'):
        original_mt = hl.read_matrix_table(mt_path)
        original_mt = original_mt.checkpoint(get_custom_ukb_pheno_mt_path(f'{suffix}_before_{curdate}'), overwrite=overwrite)
        original_mt.cols().export(get_custom_phenotype_summary_backup_path(suffix, curdate))
        original_mt.union_cols(mt, row_join_type='outer').write(mt_path, overwrite=overwrite)
    else:
        mt.write(mt_path, overwrite=overwrite)
    custom_summarize_data(suffix, overwrite=overwrite)


def export_pheno_serial_custom(output_path, pheno_keys, proportion_single_sex, suffix, pop, include_addl_covariates, pheno_sex='both_sexes'):
    binary_trait = saige_pheno_types.get(pheno_keys['trait_type']) != 'quantitative'
    mt = get_custom_ukb_pheno_mt(suffix, pop, include_addl_covariates)
    mt = mt.filter_cols(hl.all(lambda x: x, [mt[k] == pheno_keys[k] for k in PHENO_KEY_FIELDS if k != 'pheno_sex']))
    pheno_sex_mt = mt.filter_cols(mt.pheno_sex == pheno_sex)
    if pheno_sex_mt.count_cols() == 1:
        mt = pheno_sex_mt
    else:
        mt = mt.filter_cols(mt.pheno_sex == 'both_sexes')
    mt = mt.select_entries(value=mt[pheno_sex])
    if binary_trait:
        mt = mt.select_entries(value=hl.int(mt.value))
    if proportion_single_sex > 0:
        prop_female = mt.n_cases_females / (mt.n_cases_males + mt.n_cases_females)
        prop_female = prop_female.collect()[0]
        print(f'Female proportion: {prop_female}')
        if prop_female <= proportion_single_sex:
            print(f'{prop_female} less than {proportion_single_sex}. Filtering to males...')
            mt = mt.filter_rows(mt.sex == 1)
        elif prop_female >= 1 - proportion_single_sex:
            print(f'{prop_female} greater than {1 - proportion_single_sex}. Filtering to females...')
            mt = mt.filter_rows(mt.sex == 0)
    ht = mt.key_cols_by().select_cols().entries()
    ht.export(output_path)


def export_pheno_custom(p: Batch, output_path: str, pheno_keys, module: str, mt_loading_function: str,
                        docker_image: str, proportion_single_sex: float = 0.1, n_threads: int = 8, storage: str = '500Mi',
                        additional_args: str = ''):
    extract_task: Job = p.new_job(name='extract_pheno', attributes=copy.deepcopy(pheno_keys))
    extract_task.image(docker_image).cpu(n_threads).storage(storage)
    pheno_dict_opts = ' '.join([f"--{k} {shq(v)}" for k, v in pheno_keys.items()])
    python_command = f"""set -o pipefail; python3.8 {SCRIPT_DIR}/export_pheno.py
    --load_module {module} --load_mt_function {mt_loading_function}
    {pheno_dict_opts} {"--binary_trait" if saige_pheno_types.get(pheno_keys['trait_type']) != 'quantitative' else ""}
    --proportion_single_sex {proportion_single_sex}
    {"--additional_args " + additional_args if additional_args else ''}
    --output_file {extract_task.out}
    --n_threads {n_threads} | tee {extract_task.stdout}
    ; """.replace('\n', ' ')
    activate_service_account(extract_task)
    extract_task.command(python_command)

    p.write_output(extract_task.out, output_path)
    p.write_output(extract_task.stdout, f'{output_path}.log')
    return extract_task


def custom_fit_null_glmm(p: Batch, output_root: str, pheno_file: Resource, trait_type: str, covariates: str,
                         plink_file_root: str, docker_image: str, sparse_grm: Resource = None,
                         sparse_grm_extension: str = None, inv_normalize: bool = False, skip_model_fitting: bool = False,
                         min_covariate_count: int = 10,
                         n_threads: int = 16, storage: str = '10Gi', memory: str = '60G',
                         non_pre_emptible: bool = False, disable_loco: bool = False):
    analysis_type = "variant" if sparse_grm is None else "gene"
    pheno_col = 'value'
    user_id_col = 'userId'
    in_bfile = p.read_input_group(**{ext: f'{plink_file_root}.{ext}' for ext in ('bed', 'bim', 'fam')})
    fit_null_task = p.new_job(name=f'fit_null_model',
                              attributes={
                                  'analysis_type': analysis_type,
                                  'trait_type': trait_type
                              }).storage(storage).image(docker_image)
    if non_pre_emptible:
        fit_null_task._cpu = None
        fit_null_task._memory = None
        fit_null_task._machine_type = 'n1-highmem-8'
        fit_null_task._preemptible = False
    else:
        fit_null_task = fit_null_task.cpu(n_threads).memory(memory)
    output_files = {ext: f'{{root}}{ext if ext.startswith("_") else "." + ext}' for ext in
                   ('rda', f'{analysis_type}.varianceRatio.txt')}
    if analysis_type == 'gene':
        sparse_sigma_extension = sparse_grm_extension.replace("GRM", "Sigma")
        output_files[f'{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'] = \
            f'{{root}}.{analysis_type}.varianceRatio.txt{sparse_sigma_extension}'
    fit_null_task.declare_resource_group(null_glmm=output_files)
    bim_fix_command = f'perl -pi -e s/^chr// {in_bfile.bim}'
    # if trait_type == 'icd':
    #     bim_fix_command += (f"; zcat {pheno_file.gz} | perl -p -e 's/true/1/g' | perl -p -e 's/false/0/g' "
    #                         f"| gzip -c > {pheno_file.gz}.temp.gz; mv {pheno_file.gz}.temp.gz {pheno_file.gz}")
    loco_str = 'FALSE'
    command = (f'set -o pipefail; Rscript /usr/local/bin/step1_fitNULLGLMM.R '
               f'--plinkFile={in_bfile} '
               f'--phenoFile={pheno_file} '
               f'--covarColList={covariates} '
               f'--minCovariateCount={min_covariate_count} '
               f'--phenoCol={pheno_col} '
               f'--sampleIDColinphenoFile={user_id_col} '
               f'--traitType={saige_pheno_types[trait_type]} '
               f'--outputPrefix={fit_null_task.null_glmm} '
               f'--outputPrefix_varRatio={fit_null_task.null_glmm}.{analysis_type} '
               f'--skipModelFitting={str(skip_model_fitting).upper()} ')
    if inv_normalize:
        command += '--invNormalize=TRUE '
    if analysis_type == "gene":
        fit_null_task.declare_resource_group(sparse_sigma={sparse_sigma_extension: f'{{root}}.{sparse_sigma_extension}'})
        command += (f'--IsSparseKin=TRUE '
                    f'--sparseGRMFile={sparse_grm[sparse_grm_extension]} '
                    f'--sparseGRMSampleIDFile={sparse_grm[f"{sparse_grm_extension}.sampleIDs.txt"]} '
                    f'--isCateVarianceRatio=TRUE ')
    else:
        loco_str = 'FALSE' if disable_loco else 'TRUE'
    
    command += f'--nThreads={n_threads} --LOCO={loco_str} 2>&1 | tee {fit_null_task.stdout}'
    command = '; '.join([bim_fix_command, command])
    fit_null_task.command(command)
    p.write_output(fit_null_task.null_glmm, output_root)
    p.write_output(fit_null_task.stdout, f'{output_root}.{analysis_type}.log')
    # Runtimes: 8 threads: ~5 minutes of 100% CPU (~3G RAM), followed by ~9 minutes of 800% (~6G RAM)
    return fit_null_task


def custom_run_saige(p: Batch, output_root: str, model_file: str, variance_ratio_file: str,
                     vcf_file: ResourceGroup, samples_file: ResourceGroup,
                     docker_image: str,
                     group_file: str = None, sparse_sigma_file: str = None, use_bgen: bool = True,
                     trait_type: str = 'continuous',
                     chrom: str = 'chr1', min_mac: int = 1, min_maf: float = 0, max_maf: float = 0.5, cpu: int = 1,
                     memory: str = '', storage: str = '10Gi', add_suffix: str = '', log_pvalue: bool = False,
                     disable_loco: bool = False):

    analysis_type = "gene" if sparse_sigma_file is not None else "variant"
    run_saige_task: Job = p.new_job(name=f'run_saige',
                                    attributes={
                                        'analysis_type': analysis_type
                                    }).cpu(cpu).storage(storage).image(docker_image)  # Step 2 is single-threaded only

    if analysis_type == 'gene':
        run_saige_task.declare_resource_group(result={f'{add_suffix}gene.txt': '{root}',
                                                      f'{add_suffix}single.txt': '{root}_single'})
    else:
        run_saige_task.declare_resource_group(result={'single_variant.txt': '{root}'})

    loco_str = 'FALSE' if disable_loco else 'TRUE'
    command = (f'set -o pipefail; set +e; {MKL_OFF} Rscript /usr/local/bin/step2_SPAtests.R '
               f'--minMAF={min_maf} '
               f'--minMAC={min_mac} '
               f'--sampleFile={samples_file} '
               f'--GMMATmodelFile={model_file} '
               f'{"--IsOutputlogPforSingle=TRUE " if log_pvalue else ""}'
               f'--varianceRatioFile={variance_ratio_file} '
               f'--LOCO={loco_str} '
               f'--chrom={chrom} ')

    if use_bgen:
        command += (f'--bgenFile={vcf_file.bgen} '
                    f'--bgenFileIndex={vcf_file["bgen.bgi"]} ')
    else:
        command += (f'--vcfFile={vcf_file["vcf.gz"]} '
                    f'--vcfFileIndex={vcf_file["vcf.gz.tbi"]} '
                    f'--vcfField=GT ')
    if analysis_type == "gene":
        if saige_pheno_types[trait_type] == 'binary':
            command += f'--IsOutputPvalueNAinGroupTestforBinary=TRUE '
        command += (f'--groupFile={group_file} '
                    f'--sparseSigmaFile={sparse_sigma_file} '
                    f'--is_single_in_groupTest=TRUE '
                    f'--maxMAF_in_groxupTest={max_maf} '
                    f'--SAIGEOutputFile={run_saige_task.result} ')
                    #f'--IsOutputBETASEinBurdenTest=TRUE ')
    else:
        command += f'--SAIGEOutputFile=result_file '
    
    command += f' 2>&1 | tee {run_saige_task.stdout}; cat {run_saige_task.stdout}; ' # --IsOutputAFinCaseCtrl=TRUE
    if analysis_type == 'gene':
        command += f"input_length=$(wc -l {group_file} | awk '{{print $1}}'); " \
            f"output_length=$(wc -l {run_saige_task.result[f'{add_suffix}gene.txt']} | awk '{{print $1}}'); " \
            f"echo 'Got input:' $input_length 'output:' $output_length | tee -a {run_saige_task.stdout}; " \
            f"if [[ $input_length > 0 ]]; then echo 'got input' | tee -a {run_saige_task.stdout}; " \
            f"if [[ $output_length == 1 ]]; then echo 'but not enough output' | tee -a {run_saige_task.stdout}; " \
                   f"rm -f {run_saige_task.result[f'{add_suffix}gene.txt']}; exit 1; fi; fi"
    else:
        command += f"success_found=$(cat {run_saige_task.stdout} | tail -n1 | grep '^\[1\] \"Analysis done! The results have been saved to' | wc -l | awk '{{print $1}}'); " \
            f"dimnames_found=$(cat {run_saige_task.stdout} | tail -n40 | grep '^Error: length of .dimnames. \[2\] not equal to array extent' | wc -l | awk '{{print $1}}'); " \
            f"if [[ $success_found -eq 1 ]]; then echo 'Success found.' | tee -a {run_saige_task.stdout}; " \
                   f"cp result_file {run_saige_task.result}; " \
            f"elif [[ $dimnames_found -eq 1 ]]; then echo 'Failed due to dimnames.' | tee -a {run_saige_task.stdout}; exit 1; " \
            f"elif [[ $success_found -ne 1 ]]; then echo 'Failed due to other reason (aborted?).' | tee -a {run_saige_task.stdout}; exit 134; fi"
    run_saige_task.command(command)
    p.write_output(run_saige_task.result, output_root)
    p.write_output(run_saige_task.stdout, f'{output_root}.{analysis_type}.log')
    return run_saige_task


def custom_load_results_into_hail(p: Batch, output_root: str, pheno_keys, tasks_to_hold,
                                  vep_path: str, docker_image: str, gene_map_path: str = None, null_glmm_log: str = '',
                                  reference: str = 'GRCh38', saige_log: str = '', analysis_type: str = 'gene',
                                  n_threads: int = 8, storage: str = '10Gi', legacy_annotations: bool = False,
                                  log_pvalue: bool = False, overwrite: bool = True):
    load_data_task: Job = p.new_job(name=f'load_data', attributes=copy.deepcopy(pheno_keys)
                                    ).image(docker_image).cpu(n_threads).storage(storage).memory('standard')
    load_data_task.always_run().depends_on(*tasks_to_hold)
    pheno_dict_opts = ' '.join([f"--{k} {shq(v)}" for k, v in pheno_keys.items()])
    python_command = f"""python3.8 {SCRIPT_DIR}/load_results.py
    --input_dir {shq(output_root)}
    {"--null_glmm_log " + shq(null_glmm_log) if null_glmm_log else ''}
    --saige_run_log_format {saige_log}
    {pheno_dict_opts}
    {"--gene_map_ht_raw_path " + gene_map_path if gene_map_path else ''}
    {"--legacy_annotations" if legacy_annotations else ""}
    --ukb_vep_ht_path {vep_path}
    {"--overwrite" if overwrite else ""} --reference {reference}
    --analysis_type {analysis_type} {"--log_pvalue" if log_pvalue else ""}
    --n_threads {n_threads} | tee {load_data_task.stdout}
    ;""".replace('\n', ' ')

    python_command = python_command.replace('\n', '; ').strip()
    command = f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" ' + python_command
    load_data_task.command(command)
    activate_service_account(load_data_task)
    p.write_output(load_data_task.stdout, f'{output_root}/{pheno_keys["phenocode"]}_loading.log')
    return load_data_task


def custom_qq_plot_results(p: Batch, output_root: str, tasks_to_hold, export_docker_image: str, R_docker_image: str,
                    n_threads: int = 8, storage: str = '10Gi'):

    qq_export_task: Job = p.new_job(name='qq_export').image(export_docker_image).cpu(n_threads).storage(storage)
    qq_export_task.always_run().depends_on(*tasks_to_hold)

    python_command = f"""python3.8 {SCRIPT_DIR}/export_results_for_qq.py
    --input_dir {shq(output_root)}
    --output_file {qq_export_task.out}
    --n_threads {n_threads}
    ; """.replace('\n', ' ')

    command = f'set -o pipefail; PYTHONPATH=$PYTHONPATH:/ PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory={int(3 * n_threads)}g pyspark-shell" ' + python_command
    qq_export_task.command(command)
    activate_service_account(qq_export_task)

    qq_task: Job = p.new_job(name='qq_plot').image(R_docker_image).cpu(n_threads).storage(storage).always_run()
    qq_task.declare_resource_group(result={ext: f'{{root}}_Pvalue_{ext}'
                                           for ext in ('qqplot.png', 'manhattan.png', 'manhattan_loglog.png', 'qquantiles.txt')})
    R_command = f"/saige-pipelines/scripts/qqplot.R -f {qq_export_task.out} -o {qq_task.result} -p Pvalue; "
    qq_task.command(R_command)

    p.write_output(qq_task.result, output_root)
    return qq_export_task, qq_task


def custom_load_variant_data(directory: str, pheno_key_dict, ukb_vep_path: str, extension: str = 'single.txt',
                             n_cases: int = -1, n_controls: int = -1, heritability: float = -1.0,
                             saige_version: str = 'NA', inv_normalized: str = 'NA', log_pvalue: bool = False,
                             overwrite: bool = False, legacy_annotations: bool = False,
                             num_partitions: int = 1000):
    output_ht_path = f'{directory}/variant_results.ht'
    ht = hl.import_table(f'{directory}/*.{extension}', delimiter='\t', impute=True)
    print(f'Loading: {directory}/*.{extension} ...')
    marker_id_col = 'markerID' if extension == 'single.txt' else 'MarkerID'
    locus_alleles = ht[marker_id_col].split('_')
    if n_cases == -1: n_cases = hl.null(hl.tint)
    if n_controls == -1: n_controls = hl.null(hl.tint)
    if heritability == -1.0: heritability = hl.null(hl.tfloat)
    if saige_version == 'NA': saige_version = hl.null(hl.tstr)
    if inv_normalized == 'NA': inv_normalized = hl.null(hl.tstr)

    ht = ht.key_by(locus=hl.parse_locus(locus_alleles[0]), alleles=locus_alleles[1].split('/'),
                   **pheno_key_dict).distinct().naive_coalesce(num_partitions)
    if marker_id_col == 'MarkerID':
        ht = ht.drop('CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2')
    ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(
        n_cases=n_cases, n_controls=n_controls, heritability=heritability, saige_version=saige_version,
        inv_normalized=inv_normalized, log_pvalue=log_pvalue)
    ht = ht.drop('Tstat')
    ht = ht.annotate(**get_vep_formatted_data(ukb_vep_path, legacy_annotations=legacy_annotations)[
        hl.struct(locus=ht.locus, alleles=ht.alleles)])  # TODO: fix this for variants that overlap multiple genes
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite).drop('n_cases', 'n_controls', 'heritability')


def custom_get_cases_and_controls_from_log(log_format):
    """
    'gs://path/to/result_chr{chrom}_000000001.variant.log'
    """
    cases = controls = -1
    for chrom in range(10, 23):
        try:
            with hl.hadoop_open(log_format.format(chrom=chrom)) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('Analyzing'):
                        fields = line.split()
                        if len(fields) == 6:
                            try:
                                cases = int(fields[1])
                                controls = int(fields[4])
                                break
                            except ValueError:
                                logger.warn(f'Could not load number of cases or controls from {line}.')
                    elif line.endswith('samples will be used for analysis'):
                        fields = line.split()
                        try:
                            cases = int(fields[0])
                        except ValueError:
                            logger.warn(f'Could not load number of cases or controls from {line}.')
            return cases, controls
        except:
            pass
    return cases, controls


def main(args):
    hl.init(log='/tmp/saige_temp_hail.log')

    # num_pcs = 20
    num_pcs = 10
    start_time = time.time()
    if args.include_base_covariates:
        basic_covars = ['sex', 'age', 'age2', 'age_sex', 'age2_sex']
    else:
        basic_covars = []
    if args.include_addl_covariates is not None:
        addl_covars = [x for x in hl.import_table(args.include_addl_covariates).row if x not in ['s']]
    else:
        addl_covars = []
    covariates = ','.join(basic_covars + addl_covars + [f'PC{x}' for x in range(1, num_pcs + 1)])
    n_threads = args.n_threads
    analysis_type = "variant"
    chromosomes = list(map(str, range(1, 23))) + ['X']
    reference = 'GRCh37'
    chrom_lengths = hl.get_reference(reference).lengths
    iteration = 1
    pops = args.pops.split(',') if args.pops else POPS

    # if args.local_test:
    #     backend = hb.LocalBackend(gsa_key_file='/Users/konradk/.hail/ukb-diverse-pops.json')
    # else:

    # ensure the data exist
    if not hl.hadoop_exists(f'{get_custom_ukb_pheno_mt_path(args.suffix)}/_SUCCESS') or (args.overwrite_pheno_data or args.append):
        produce_custom_phenotype_mt(args.data_path, args.data_extn, args.suffix, 
                                    trait_type=args.trait_type, modifier=args.modifier,
                                    source=args.source, sample_col=args.sample_col, 
                                    append=args.append, overwrite=args.overwrite_pheno_data,
                                    custom_covars=args.include_addl_covariates)

    if not args.force_serial_phenotype_export:
        backend = hb.ServiceBackend(billing_project='ukb_diverse_pops',
                                    bucket=temp_bucket.split('gs://', 1)[-1])
    for pop in pops:
        if not args.force_serial_phenotype_export:
            p = hb.Batch(name=f'saige_pan_ancestry_{pop}', backend=backend, default_image=SAIGE_DOCKER_IMAGE,
                         default_storage='500Mi', default_cpu=n_threads)
        window = '1e7' if pop == 'EUR' else '1e6'
        logger.info(f'Setting up {pop}...')
        chunk_size = int(5e6) if pop != 'EUR' else int(1e6)
        phenos_to_run = custom_get_phenos_to_run(pop=pop, suffix=args.suffix, limit=int(args.local_test),
                                                 specific_phenos=args.phenos, single_sex_only=args.single_sex_only,
                                                 skip_case_count_filter=args.skip_case_count_filter,
                                                 sex_stratified=args.sex_stratified)
        logger.info(f'Got {len(phenos_to_run)} phenotypes...')
        if len(phenos_to_run) <= 20:
            logger.info(phenos_to_run)

        pheno_export_dir = f'{pheno_folder}/exported/{args.suffix}/{pop}'
        phenos_already_exported = {}
        if not args.overwrite_pheno_data and hl.hadoop_exists(pheno_export_dir):
            phenos_already_exported = {x['path'] for x in hl.hadoop_ls(pheno_export_dir)}
        pheno_exports = {}

        for pheno_key_dict in phenos_to_run:
            pheno_export_path = get_pheno_output_path(pheno_export_dir, pheno_key_dict, legacy=False)
            if args.force_serial_phenotype_export:
                export_pheno_serial_custom(pheno_export_path, pheno_key_dict, proportion_single_sex=0, 
                                           suffix=args.suffix, pop=pop, include_addl_covariates=args.include_addl_covariates)
                pheno_file = ''
            elif not args.overwrite_pheno_data and pheno_export_path in phenos_already_exported:
                pheno_file = p.read_input(pheno_export_path)
            else:
                # NOTE need to update module and function name
                addlargs = f'{args.suffix},{pop}' if args.include_addl_covariates is None else f'{args.suffix},{pop},{args.include_addl_covariates}'
                pheno_task = export_pheno_custom(p, pheno_export_path, pheno_key_dict, 'ukbb_pan_ancestry.saige_pan_ancestry_custom', 'get_custom_ukb_pheno_mt', 
                                                 PHENO_DOCKER_IMAGE, additional_args=addlargs, n_threads=args.n_cpu_pheno, proportion_single_sex=0)
                pheno_task.attributes.update({'pop': pop})
                pheno_file = pheno_task.out
            pheno_exports[stringify_pheno_key_dict(pheno_key_dict)] = pheno_file
        completed = Counter([isinstance(x, hb.resource.InputResourceFile) for x in pheno_exports.values()])
        logger.info(f'Exporting {completed[False]} phenos (already found {completed[True]})...')
        
        if args.force_serial_phenotype_export:
            break
        
        overwrite_null_models = args.create_null_models
        null_model_dir = f'{root}/null_glmm/{args.suffix}/{pop}'
        null_models_already_created = {}
        if not overwrite_null_models and hl.hadoop_exists(null_model_dir):
            null_models_already_created = {x['path'] for x in hl.hadoop_ls(null_model_dir)}
        null_models = {}

        for pheno_key_dict in phenos_to_run:
            null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, '', legacy=False)
            model_file_path = f'{null_glmm_root}.rda'
            variance_ratio_file_path = f'{null_glmm_root}.{analysis_type}.varianceRatio.txt'

            if not overwrite_null_models and model_file_path in null_models_already_created and \
                    variance_ratio_file_path in null_models_already_created:
                model_file = p.read_input(model_file_path)
                variance_ratio_file = p.read_input(variance_ratio_file_path)
            else:
                if args.skip_any_null_models: continue
                fit_null_task = custom_fit_null_glmm(p, null_glmm_root, pheno_exports[stringify_pheno_key_dict(pheno_key_dict)],
                                                     pheno_key_dict['trait_type'], covariates,
                                                     get_ukb_grm_plink_path(pop, iteration, window), SAIGE_DOCKER_IMAGE,
                                                     inv_normalize=args.force_inv_normalize, n_threads=n_threads, min_covariate_count=1,
                                                     non_pre_emptible=args.non_pre_emptible, storage='100Gi', disable_loco=args.disable_loco)
                fit_null_task.attributes.update({'pop': pop})
                fit_null_task.attributes.update(copy.deepcopy(pheno_key_dict))
                model_file = fit_null_task.null_glmm.rda
                variance_ratio_file = fit_null_task.null_glmm[f'{analysis_type}.varianceRatio.txt']
            null_models[stringify_pheno_key_dict(pheno_key_dict)] = (model_file, variance_ratio_file)

        completed = Counter([isinstance(x[0], hb.resource.InputResourceFile) for x in null_models.values()])
        logger.info(f'Running {completed[False]} null models (already found {completed[True]})...')

        use_bgen = True
        vcf_dir = f'{root_vcf}/vcf/{pop}'
        test_extension = 'bgen' if use_bgen else 'vcf.gz'
        overwrite_vcfs = args.create_vcfs
        vcfs_already_created = {}
        if not overwrite_vcfs and hl.hadoop_exists(vcf_dir):
            vcfs_already_created = {x['path'] for x in hl.hadoop_ls(vcf_dir)}
            logger.info(f'Found {len(vcfs_already_created)} VCFs in directory...')
        vcfs = {}
        for chromosome in chromosomes:
            chrom_length = chrom_lengths[chromosome]
            for start_pos in range(1, chrom_length, chunk_size):
                end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                interval = f'{chromosome}:{start_pos}-{end_pos}'
                vcf_root = f'{vcf_dir}/variants_{chromosome}_{str(start_pos).zfill(9)}'
                if f'{vcf_root}.{test_extension}' in vcfs_already_created:
                    if use_bgen:
                        vcf_file = p.read_input_group(**{'bgen': f'{vcf_root}.bgen',
                                                         'bgen.bgi': f'{vcf_root}.bgen.bgi',
                                                         'sample': f'{vcf_root}.sample'})
                    else:
                        vcf_file = p.read_input_group(**{'vcf.gz': f'{vcf_root}.vcf.gz',
                                                         'vcf.gz.tbi': f'{vcf_root}.vcf.gz.tbi'})
                else:
                    vcf_task = extract_vcf_from_mt(p, vcf_root, HAIL_DOCKER_IMAGE, 'ukbb_pan_ancestry', adj=False,
                                                   additional_args=f'{chromosome},{pop}', input_dosage=True,
                                                   reference=reference, interval=interval, export_bgen=use_bgen,
                                                   n_threads=n_threads)
                    vcf_task.attributes['pop'] = pop
                    vcf_file = vcf_task.out
                vcfs[interval] = vcf_file
                if args.local_test:
                    break
            if args.local_test:
                break

        completed = Counter([type(x) == hb.resource.InputResourceFile for x in vcfs.values()])
        logger.info(f'Creating {completed[False]} VCFs (already found {completed[True]})...')

        result_dir = f'{root}/result/{args.suffix}/{pop}'
        overwrite_results = args.overwrite_results
        for i, pheno_key_dict in enumerate(phenos_to_run):
            if stringify_pheno_key_dict(pheno_key_dict) not in null_models: continue
            model_file, variance_ratio_file = null_models[stringify_pheno_key_dict(pheno_key_dict)]

            if not i % 10:
                n_jobs = dict(Counter(map(lambda x: x.name, p.select_jobs("")))).get("run_saige", 0)
                logger.info(f'Read {i} phenotypes ({n_jobs} new to run so far)...')

            pheno_results_dir = get_pheno_output_path(result_dir, pheno_key_dict, '', legacy=False)
            results_already_created = {}

            if not overwrite_results and not args.skip_saige and hl.hadoop_exists(pheno_results_dir):
                results_already_created = {x['path'] for x in hl.hadoop_ls(pheno_results_dir)}

            saige_tasks = []
            for chromosome in chromosomes:
                if args.skip_saige: break
                chrom_length = chrom_lengths[chromosome]
                for start_pos in range(1, chrom_length, chunk_size):
                    end_pos = chrom_length if start_pos + chunk_size > chrom_length else (start_pos + chunk_size)
                    interval = f'{chromosome}:{start_pos}-{end_pos}'
                    vcf_file = vcfs[interval]
                    results_path = get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome, start_pos, legacy=False)
                    if overwrite_results or f'{results_path}.single_variant.txt' not in results_already_created:
                        samples_file = p.read_input(get_ukb_samples_file_path(pop, iteration))
                        saige_task = custom_run_saige(p, results_path, model_file, variance_ratio_file, vcf_file, samples_file,
                                                      SAIGE_DOCKER_IMAGE, trait_type=pheno_key_dict['trait_type'], use_bgen=use_bgen,
                                                      chrom=chromosome, log_pvalue=args.log_p, disable_loco=args.disable_loco,
                                                      cpu=args.n_cpu_saige, min_mac=args.min_mac)
                        saige_task.attributes.update({'interval': interval, 'pop': pop})
                        saige_task.attributes.update(copy.deepcopy(pheno_key_dict))
                        saige_tasks.append(saige_task)
                    if args.local_test:
                        break
                if args.local_test:
                    break

            res_tasks = []
            if overwrite_results or args.overwrite_hail_results or \
                    f'{pheno_results_dir}/variant_results.ht' not in results_already_created or \
                    not hl.hadoop_exists(f'{pheno_results_dir}/variant_results.ht/_SUCCESS'):
                null_glmm_root = get_pheno_output_path(null_model_dir, pheno_key_dict, f'.{analysis_type}.log',
                                                       legacy=False)

                prefix = get_results_prefix(pheno_results_dir, pheno_key_dict,
                                            f'{"chr" if reference == "GRCh38" else ""}{{chrom}}', 1,
                                            legacy=False)
                saige_log = f'{prefix}.{analysis_type}.log'

                load_task = custom_load_results_into_hail(p, pheno_results_dir, pheno_key_dict,
                                                          saige_tasks, get_ukb_vep_path(), PHENO_DOCKER_IMAGE,
                                                          analysis_type=analysis_type, saige_log=null_glmm_root,#saige_log, 
                                                          n_threads=args.n_cpu_merge, null_glmm_log=null_glmm_root,
                                                          reference=reference, legacy_annotations=True, log_pvalue=args.log_p)
                load_task.attributes['pop'] = pop
                res_tasks.append(load_task)
                qq_export, qq_plot = custom_qq_plot_results(p, pheno_results_dir, res_tasks, PHENO_DOCKER_IMAGE, QQ_DOCKER_IMAGE, n_threads=n_threads)
                qq_export.attributes.update({'pop': pop})
                qq_export.attributes.update(copy.deepcopy(pheno_key_dict))
                qq_plot.attributes.update({'pop': pop})
                qq_plot.attributes.update(copy.deepcopy(pheno_key_dict))

        logger.info(f'Setup took: {time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))}')
        logger.info(f'Submitting: {get_tasks_from_pipeline(p)}')
        logger.info(f"Total size: {sum([len(x._pretty()) for x in p.select_jobs('')])}")
        p.run(dry_run=args.dry_run, wait=False, delete_scratch_on_exit=False)
        logger.info(f'Finished: {get_tasks_from_pipeline(p)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--single_sex_only', help='Run only single sex phenotypes (experimental)', action='store_true')
    parser.add_argument('--sex_stratified', help='Run these phenotypes in a sex-stratified fashion (experimental)', choices=(None, 'all', 'only'))
    parser.add_argument('--skip_any_null_models', help='Skip running SAIGE null models', action='store_true')
    parser.add_argument('--skip_saige', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--create_null_models', help='Force creation of null models', action='store_true')
    parser.add_argument('--create_vcfs', help='Force creation of VCFs', action='store_true')
    parser.add_argument('--overwrite_pheno_data', help='Overwrite phenotype munged data and exports', action='store_true')
    parser.add_argument('--overwrite_results', help='Force run of SAIGE tests', action='store_true')
    parser.add_argument('--overwrite_hail_results', help='Force run of results loading', action='store_true')
    parser.add_argument('--local_test', help='Local test of pipeline', action='store_true')
    parser.add_argument('--non_pre_emptible', help='Local test of pipeline', action='store_true')
    parser.add_argument('--skip_case_count_filter', help='Skip running SAIGE tests', action='store_true')
    parser.add_argument('--phenos', help='Comma-separated list of trait_type-phenocode-pheno_sex-coding-modifier regexes '
                                         '(e.g. continuous-50-both_sexes--,icd10-E1.*,brain_mri-.* )')
    parser.add_argument('--pops', help='comma-searated list')
    parser.add_argument('--dry_run', help='Dry run only', action='store_true')
    parser.add_argument('--include_base_covariates', help='If true, will include the usual age, sex covariates.', action='store_true')
    parser.add_argument('--include_addl_covariates', help='Points to a .tsv with additional covariates. All will be included. Must have sample key as s.', type=str)
    parser.add_argument('--data_path', required=True, type=str, help='Path to phenotype data.')
    parser.add_argument('--sample_col', default='s', type=str, help='Sample identifier column.')
    parser.add_argument('--data_extn', required=True, type=str, help='Type of the phenotype data (ht or tsv).')
    parser.add_argument('--suffix', required=True, type=str, help='Analysis suffix.')
    parser.add_argument('--trait_type', required=True, type=str, help='Trait type; also used for munging.')
    parser.add_argument('--modifier', required=True, type=str, help='Global modifier for this dataset used for phenotype munging.')
    parser.add_argument('--source', required=True, type=str, help='Data source for tagging, used for pheontype munging.')
    parser.add_argument('--append', action='store_true', help='If enabled, will attempt to augment the pheno MT with new phenotypes.')
    parser.add_argument('--disable_loco', action='store_true', help='If enabled, will avoid LOCO. This was the original pan-ancestry pipeline behavior.')
    parser.add_argument('--log_p', action='store_true', help='Log transform p-values via SAIGE.')
    parser.add_argument('--n_cpu_saige', type=int, default=1, help='Number of threads to request for running SAIGE.')
    parser.add_argument('--n_cpu_merge', type=int, default=8, help='Number of threads to request during summary stat merging.')
    parser.add_argument('--n_cpu_pheno', default=16, type=int)
    parser.add_argument('--n_threads', default=8, type=int)
    parser.add_argument('--min_mac', default=1, type=int)
    parser.add_argument('--force_inv_normalize', action='store_true')
    parser.add_argument('--force_serial_phenotype_export', action='store_true')
    args = parser.parse_args()

    main(args)
