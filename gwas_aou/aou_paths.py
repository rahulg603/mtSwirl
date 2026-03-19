import os


# directories
BUCKET = os.getenv('WORKSPACE_BUCKET')
DATASET = os.getenv("WORKSPACE_CDR")
PCA_DIR = f'{BUCKET}/pca/'
FILTERED_SAMPLES_DIR = f'{PCA_DIR}filtered_samples/'
COVARIATES_DIR = f'{BUCKET}/covariates/'
GWAS_DIR = f'{BUCKET}/gwas/'
RESULTS_DIR = f'{BUCKET}/gwas/results/'
TEMP = f'{BUCKET}/tmp/'
CDR_v6_WGS = 'gs://fc-aou-datasets-controlled/v6/wgs/'


# paths
COVAR_PATH = f'{COVARIATES_DIR}final_covariates.ht'
ANCESTRY_INFO_PATH = f"{CDR_v6_WGS}vcf/aux/ancestry/ancestry_preds.tsv"
REL_SAMP_PATH = f"{CDR_v6_WGS}vcf/aux/relatedness/relatedness_flagged_samples.tsv"
FLAGGED_SAMP_PATH = f"{CDR_v6_WGS}vcf/aux/qc/flagged_samples.tsv"
FINAL_VARIANT_CALLSET_PATH = f'{BUCKET}/final_callset_220920/vcf/221012_filt_annotated/annotated_combined_processed_flat.tsv.bgz'
FINAL_SAMPLE_STATS = f'{BUCKET}/final_callset_220920/221012_filtered_aou_tab_per_sample_stats.tsv'
SANITY_PHENOS_PATH = f'{GWAS_DIR}sanity_phenotypes/sanity_check_demographic_phenos.ht'
SANITY_PHENOS_COVARS_PATH = f'{GWAS_DIR}sanity_phenotypes/sanity_check_demographic_phenos_covars.ht'



def get_final_case_only_hl_path(num_to_keep, type):
    base_path = f'{BUCKET}/final_callset_220920/vcf/221012_filt_annotated/munged/'
    if type == 'ht':
        return base_path + f'case_only_calls_low_hl_fullqc_N_{str(num_to_keep)}.ht'
    elif type == 'tsv':
        return base_path + f'case_only_calls_low_hl_fullqc_long_format_N_{str(num_to_keep)}.tsv'


def get_lambdas_path(suffix, pop, extn):
    return f'{RESULTS_DIR}/lambdas/{suffix}/lambda_export_{pop}.{extn}'