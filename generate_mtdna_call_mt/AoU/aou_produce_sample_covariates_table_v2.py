#!/usr/bin/env python3
import hail as hl
import pandas as pd
import argparse
import os
import re

# AOU CONSTANTS
ARRAY_PATH = os.getenv('MICROARRAY_HAIL_STORAGE_PATH')
EXOME_PATH = os.getenv('WGS_EXOME_SPLIT_HAIL_PATH')
WGS_PATH = os.getenv('WGS_ACAF_THRESHOLD_SPLIT_HAIL_PATH')
AUX_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux'

# WORKSPACE PATHS
BUCKET = os.getenv('WORKSPACE_BUCKET')
REFERENCE_PATH = os.path.join(BUCKET, 'reference_files')
GWAS_PATH = os.path.join(BUCKET, 'saige_gwas')
COVAR_PATH = os.path.join(GWAS_PATH, 'covariates')
TEMP_PATH = os.path.join(BUCKET, 'tmp')

def get_path_vat():
    return os.path.join(AUX_PATH, 'vat/vat_complete_v7.1.bgz.tsv.gz')

def get_ancestry_flat_file():
    return os.path.join(AUX_PATH, 'ancestry/ancestry_preds.tsv')

def get_genomic_metric_flat_file():
    return os.path.join(AUX_PATH, 'qc/genomic_metrics.tsv')

def get_aou_flagged_samples():
    return os.path.join(AUX_PATH, 'qc/flagged_samples.tsv')

def get_all_by_all_qc_samples():
    return os.path.join(REFERENCE_PATH, '241021_allbyall_qc_sample_ids.tsv')

def get_demographics_path(cov_folder, drc):
    drc_string = '_drc' if drc else '_axaou'
    return os.path.join(cov_folder, f'base/ht/all_covariates{drc_string}.ht')

# HELPER FUNCTIONS
def generate_indicator(ht, col, baseline_item):
    assert(ht[col].dtype == hl.tstr)

    all_values = ht.aggregate(hl.agg.collect_as_set(ht[col]))
    if baseline_item not in all_values:
        raise ValueError('Item to be used as baseline is not present in data.')
    return ht.annotate(**{f'{col}_{x}': ht[col] == x for x in all_values if x != baseline_item})

def load_aou_genotyping_qc_data():
    ht_qc_flat = hl.import_table(get_genomic_metric_flat_file(),
                                 impute=True, 
                                 types={"research_id":"tstr",
                                        "verify_bam_id2_contamination":"tfloat"},
                                 min_partitions=50,
                                 missing=['NA',''])
    ht_qc_flat = ht_qc_flat.annotate(isFemale = hl.if_else(ht_qc_flat.sex_at_birth == 'F', 1, 
                                                hl.if_else(ht_qc_flat.sex_at_birth == 'M', 0, hl.missing(hl.tint32))))
    ht_qc_flat = ht_qc_flat.rename({'research_id': 's'}).key_by('s')
    return ht_qc_flat

def load_axaou_sample_qc():
    ht = hl.import_table(get_all_by_all_qc_samples(),
                         impute=False,
                         types={"person_id":"tstr"})
    ht = ht.rename({'person_id': 's'}).key_by('s')
    return ht


def load_aou_flagged_samples():
    ht = hl.import_table(get_aou_flagged_samples(),
                         impute=True,
                         types={"s":"tstr",
                                "probabilities":hl.tarray(hl.tfloat),
                                "pca_features":hl.tarray(hl.tfloat),
                                "qc_metrics_filters":hl.tarray(hl.tstr)},
                        key='s')
    return ht


def load_ancestry_data(use_drc_data=True):
    # if use_drc_data is True, will return the first 15 covariates as produced by DRC.
    # if use_drc_data is False, will return the first 20 covariates as produced for AllxAllofUs
    if use_drc_data:
        npcs = 15
        ancestry_pred = hl.import_table(get_ancestry_flat_file(),
                                        key="research_id", 
                                        impute=True, 
                                        types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)},
                                        min_partitions=50)
        ancestry_pred = ancestry_pred.select(pop = ancestry_pred.ancestry_pred,
                                             **{f'PC{str(idx+1)}': ancestry_pred.pca_features[idx] for idx in range(0,npcs+1)})
        ancestry_pred = ancestry_pred.rename({'research_id': 's'}).key_by('s')
        return npcs, ancestry_pred
    else:
        # TODO AxAoU data
        raise NotImplementedError('AxAoU ancestry data not implemented yet.')


def get_all_demographics(overwrite=False, use_drc_ancestry_data=True):
    covar_path = get_demographics_path(COVAR_PATH, use_drc_ancestry_data)
    if not overwrite and hl.hadoop_exists(os.path.join(covar_path, '_SUCCESS')):
        sample_covariates = hl.read_table(covar_path)
    else:
        npcs, ancestry_pred = load_ancestry_data(use_drc_ancestry_data)
        dataset = os.getenv('WORKSPACE_CDR')
        person_sql = f"""
        SELECT  person.person_id,
                person.birth_datetime,
                p_location_concept.concept_name as loc,
                p_site_concept.concept_name as site
            FROM
                `{dataset}.person` person 
            LEFT JOIN
                `{dataset}.concept` p_location_concept 
                    on person.location_id = p_location_concept.CONCEPT_ID 
            LEFT JOIN
                `{dataset}.concept` p_site_concept 
                    on person.care_site_id = p_site_concept.CONCEPT_ID
            WHERE
                person.PERSON_ID IN (
                    select
                        person_id  
                    from
                        `{dataset}.cb_search_person` cb_search_person  
                    where
                        cb_search_person.person_id in (
                            select
                                person_id 
                            from
                                `{dataset}.cb_search_person` p 
                            where
                                has_whole_genome_variant = 1 
                        ) 
                    )"""

        wgs_demog = pd.read_gbq(person_sql, dialect="standard")
        wgs_demog['birth_year'] = pd.DatetimeIndex(wgs_demog.birth_datetime).year
        age_ht = hl.Table.from_pandas(wgs_demog[['person_id', 'birth_year']])
        age_ht = age_ht.annotate(s = hl.str(age_ht.person_id)).key_by('s')

        ht_qc_flat = load_aou_genotyping_qc_data()
        axaou_sample_qc = load_axaou_sample_qc()
        flagged_samples = load_aou_flagged_samples()

        # munge covariates
        sample_covariates = ht_qc_flat.select(isFemale = ht_qc_flat.isFemale, 
                                              nucdna_mean_coverage = ht_qc_flat.mean_coverage,
                                              site_id = ht_qc_flat.site_id,
                                              sample_collection_date = ht_qc_flat.biosample_collection_date,
                                              sample_source = ht_qc_flat.sample_source,
                                              ploidy = ht_qc_flat.dragen_sex_ploidy,
                                              contamination = hl.struct(dragen = ht_qc_flat.dragen_contamination,
                                                                     verify_bam_id2 = ht_qc_flat.verify_bam_id2_contamination))
        sample_covariates = sample_covariates.annotate(sample_collection_split = sample_covariates.sample_collection_date.split('-'))
        sample_covariates = sample_covariates.annotate(sample_collection_year = hl.if_else(hl.len(sample_covariates.sample_collection_split) == 3,
                                                                                           hl.int32(sample_covariates.sample_collection_split[0]),
                                                                                           hl.missing('int32'))).drop('sample_collection_split')
        
        # add birth year and compute age
        sample_covariates = sample_covariates.annotate(birth_year = age_ht[sample_covariates.s].birth_year)
        sample_covariates = sample_covariates.annotate(age = sample_covariates.sample_collection_year - sample_covariates.birth_year)

        # add GWAS covariates
        sample_covariates = sample_covariates.annotate(age_isFemale = sample_covariates.age * sample_covariates.isFemale,
                                                       age2 = sample_covariates.age**2)
        sample_covariates = sample_covariates.annotate(age2_isFemale = sample_covariates.age2 * sample_covariates.isFemale)
        sample_covariates = sample_covariates.annotate(ancestry = hl.struct(pop = ancestry_pred[sample_covariates.s].pop,
                                                                            **{f'PC{str(idx+1)}': ancestry_pred[sample_covariates.s][f'PC{str(idx+1)}'] for idx in range(0,npcs+1)}))

        # add flag for passing axaou
        sample_covariates = sample_covariates.annotate(pass_axaou_qc = hl.is_defined(axaou_sample_qc[sample_covariates.s]))

        # add flag from AoU DRC
        sample_covariates = sample_covariates.annotate(pass_aou_drc_qc = ~hl.is_defined(flagged_samples[sample_covariates.s]),
                                                       failed_flags_drc_qc = flagged_samples[sample_covariates.s].qc_metrics_filters)
        sample_covariates = sample_covariates.annotate(pass_qc = sample_covariates.pass_axaou_qc & sample_covariates.pass_aou_drc_qc)

        # output
        sample_covariates = sample_covariates.repartition(100).checkpoint(covar_path, overwrite=True)
    
    return sample_covariates


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--export_tsv_path', help='Path to tsv file where the covariates should be exported to the local file system.')
    parser.add_argument('--use_drc_ancestry_data', action='store_true')
    parser.add_argument('--overwrite', action='store_true', help='If supplied, will overwrite GCP hts')
    args = parser.parse_args()

    ht = get_all_demographics(overwrite=args.overwrite, use_drc_ancestry_data=args.use_drc_ancestry_data)
    ht = ht.select(sex = ht.isFemale,
                   age = ht.age,
                   pop = ht.ancestry.pop)
    if args.export_tsv_path:
        ht.export(f'file://{os.path.abspath(args.export_tsv_path)}')
