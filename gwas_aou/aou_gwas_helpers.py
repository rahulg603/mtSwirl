import hail as hl
import pandas as pd
import re
from .aou_paths import *
from .run_per_ancestry_pca import get_custom_pc_path, make_iteration_suffix


LIFTOVERFILTERS = set(['NoTarget','MismatchedRefAllele','IndelStraddlesMultipleIntevals'])
CUSTOMLIFTOVERFILTERS = set(['FailedPicardLiftoverVcf', 'InsertionSharesForceCalledDeletion', 'InsertionSharesForceCalledInsertion',
                             'AddGRCh38RefDeleToRefSiteIns', 'ComplexSwapField', 'NewInsertionHaplotype',
                             'SwapFirstAlleleIndel', 'ReplaceInternalBaseDeletion', 'FancyFieldInversion',
                             'DeletionSpannedHomoplasmicInsertion', 'LiftoverSuccessEntrySwap', 'ForceCalledHomoplasmy',
                             'LeftShiftedIndel', 'FailedDuplicateVariant'])
ALL_FILTERS = hl.literal(CUSTOMLIFTOVERFILTERS.union(LIFTOVERFILTERS))
HAP_CUTOFF = 30


def munge_phenotype_name(expr):
    """ Applies basic transformations to make phenotype names legible
    """
    return expr.replace('[\\[\\]\\(\\)/\s,.\\\\:\'\"]{1,}', '_').lower()


def hap_to_dummy(ht, field_name='hap', count_cutoff=HAP_CUTOFF):
    counts = ht.aggregate(hl.agg.counter(ht[field_name]))
    haps_to_keep = [hap for hap, ct in counts.items() if ct > count_cutoff]
    ht = ht.filter(hl.literal(haps_to_keep).contains(ht[field_name]))
    unique_haps = list(ht.aggregate(hl.agg.collect_as_set(ht[field_name])))
    ht = ht.annotate(**{f'{field_name}_{hap}': hl.if_else(ht[field_name] == hap, 1, 0) for hap in unique_haps})
    return ht.drop(field_name)


def get_covariates(iteration, overwrite=False, use_custom_pca=False):
    if use_custom_pca:
        covar_suff = f'final_covariates_custom{make_iteration_suffix(iteration)}.ht'
    else:
        covar_suff = f'final_covariates_original_PCs.ht'
    covar_path = f'{COVARIATES_DIR}{covar_suff}'
    
    if hl.hadoop_exists(f'{covar_path}/_SUCCESS') and not overwrite:
        sample_covariates = hl.read_table(covar_path)
    else:
        ht_samp_flat = hl.import_table(FINAL_SAMPLE_STATS,
                                       key="s", 
                                       impute=True, 
                                       types={"s":"tstr"},
                                       min_partitions=50)
        ht_samp_flat = ht_samp_flat.annotate(isFemale = hl.if_else(ht_samp_flat.sex_at_birth == 'Female', 1, 
                                                        hl.if_else(ht_samp_flat.sex_at_birth == 'Male', 0, hl.missing(hl.tint32))))

        person_sql = f"""
        SELECT  person.person_id,
                person.birth_datetime,
                p_location_concept.concept_name as loc,
                p_site_concept.concept_name as site,
                state.state_of_residence_concept_id as state_of_residence_concept_id,
                state.state_of_residence_source_value as state_of_residence_source_value
            FROM
                `{DATASET}.person` person 
            LEFT JOIN
                `{DATASET}.concept` p_location_concept 
                    on person.location_id = p_location_concept.CONCEPT_ID 
            LEFT JOIN
                `{DATASET}.concept` p_site_concept 
                    on person.care_site_id = p_site_concept.CONCEPT_ID
            LEFT JOIN
                `{DATASET}.person_ext` state 
                    on person.person_id = state.person_id
            WHERE
                person.PERSON_ID IN (
                    select
                        person_id  
                    from
                        `{DATASET}.cb_search_person` cb_search_person  
                    where
                        cb_search_person.person_id in (
                            select
                                person_id 
                            from
                                `{DATASET}.cb_search_person` p 
                            where
                                has_whole_genome_variant = 1 
                        ) 
                    )"""

        wgs_demog = pd.read_gbq(person_sql, dialect="standard")

        sample_covariates = ht_samp_flat.select(isFemale = ht_samp_flat.isFemale, 
                                                mtdna_mean_coverage = ht_samp_flat.mean_coverage, 
                                                nucdna_mean_coverage = ht_samp_flat.nuc_mean_coverage,
                                                major_haplogroup = ht_samp_flat.major_haplogroup,
                                                hap = ht_samp_flat.hap,
                                                mtcn = ht_samp_flat.mtcn, mtcn_median = ht_samp_flat.mtcn_median)
        sample_covariates = sample_covariates.annotate(mtcn_hi = hl.int32(sample_covariates.mtcn > 140))
        wgs_demog['year_of_birth'] = pd.DatetimeIndex(wgs_demog.birth_datetime).year
        wgs_demog['approx_age'] = 2021 - wgs_demog['year_of_birth']
        age_table = wgs_demog[['person_id', 'approx_age', 'year_of_birth']]
        age_ht = hl.Table.from_pandas(age_table)
        age_ht = age_ht.annotate(person_id = hl.str(age_ht.person_id)).key_by('person_id')
        sample_covariates = sample_covariates.annotate(approx_age = age_ht[sample_covariates.s].approx_age,
                                                       year_of_birth = age_ht[sample_covariates.s].year_of_birth)
        sample_covariates = sample_covariates.annotate(age_isFemale = sample_covariates.approx_age * sample_covariates.isFemale,
                                                       age2 = sample_covariates.approx_age**2)
        sample_covariates = sample_covariates.annotate(age2_isFemale = sample_covariates.age2 * sample_covariates.isFemale)

        if not use_custom_pca:
            ancestry_pred = hl.import_table(ANCESTRY_INFO_PATH,
                                            key="research_id", 
                                            impute=True, 
                                            types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)},
                                            min_partitions=50)
            sample_covariates = sample_covariates.annotate(**{f'PC{str(idx+1)}': ancestry_pred[sample_covariates.s].pca_features[idx] for idx in range(0,16)})
            sample_covariates = sample_covariates.annotate(pop = ancestry_pred[sample_covariates.s].ancestry_pred)
        else:
            ancestry_ht = hl.read_table(get_custom_pc_path(iteration))
            sample_covariates = sample_covariates.annotate(**ancestry_ht[sample_covariates.s])
        
        sample_covariates = hap_to_dummy(sample_covariates)
        sample_covariates = sample_covariates.repartition(100).checkpoint(covar_path, overwrite=True)
    
    return sample_covariates


def get_raw_mtdna_callset():
    mt_pheno = hl.read_matrix_table(f"{BUCKET}/final_callset_220920/vcf/full_aou_variant_callset.mt")
    mt_pheno = mt_pheno.annotate_entries(FT_LIFT = mt_pheno.FT.intersection(ALL_FILTERS))
    mt_pheno = mt_pheno.annotate_entries(FT = mt_pheno.FT.difference(ALL_FILTERS))    
    mt_pheno = mt_pheno.annotate_entries(FT = mt_pheno.FT.difference(hl.literal(set(['PASS']))))
    mt_pheno = mt_pheno.annotate_entries(FT = hl.if_else(hl.len(mt_pheno.FT) == 0, {"PASS"}, mt_pheno.FT))
    mt_pheno_mod = mt_pheno.filter_entries(mt_pheno.FT == {"PASS"})
    mt_pheno_mod = mt_pheno_mod.annotate_rows(ct_hl = hl.agg.count_where((mt_pheno_mod.HL > 0.05) & (mt_pheno_mod.HL < 0.95)))
    return mt_pheno_mod


def get_case_only_mtdna_callset(num_to_keep=310, overwrite=False):
    """ Here we implement the munging pipeline in Python and Pandas.
    """
    if hl.hadoop_exists(f"{get_final_case_only_hl_path(num_to_keep, 'ht')}/_SUCCESS") and not overwrite:
        ht = hl.read_table(get_final_case_only_hl_path(num_to_keep, 'ht'))
    
    else:
        df = pd.read_csv(FINAL_VARIANT_CALLSET_PATH, sep='\t', 
                         compression='gzip',
                         dtype={'locus': 'str', 'alleles': 'str', 's': 'str',
                                'rsid': 'str', 'variant':'str', 'batch':'str',
                                'common_low_heteroplasmy': 'boolean',
                                'hap_defining_variant': 'boolean',
                                'region': 'str',
                                'variant_context': 'str',
                                'dp_mean': 'float64',
                                'mq_mean': 'float64', 'tlod_mean': 'float64',
                                'AF_hom': 'float64', 'AF_het': 'float64',
                                'AC_hom': 'int64', 'AC_het': 'int64',
                                'max_hl': 'float64', 'dp_mean': 'float64',
                                'major_haplogroup':'str', 'hap':'str', 
                                'wgs_median_coverage':'float64', 'mt_mean_coverage':'float64',
                                'mito_cn':'float64', 'age':'float64', 'pop':'str',
                                'DP':'float64', 'AD_ref':'float64','AD_alt':'float64',
                                'MQ':'float64','TLOD':'float64', 'GT':'str',
                                'OriginalSelfRefAlleles':'str', 'SwappedFieldIDs':'str', 'FT':'str', 'FT_LIFT':'str',
                                'artifact_prone':'boolean', 'lifted':'boolean', 'fail_gt':'boolean', 'missing_call':'boolean'})
        df_for_enumeration = df[(~df['HL'].isna()) & (df['HL'] < 0.95) & (df['common_low_heteroplasmy'])]
        counted_variants = df_for_enumeration.groupby('variant')['HL'].count().sort_values(ascending=False)
        variants_to_extract = list(counted_variants[counted_variants > num_to_keep].index)

        variants_to_analyze = df[df['variant'].isin(variants_to_extract) & (~df['HL'].isna()) & (df['HL'] < 0.95)]
        variants_to_analyze.to_csv(get_final_case_only_hl_path(num_to_keep, 'tsv'), sep='\t', index=False)
        
        pivoted_vars = variants_to_analyze.pivot(index='s', columns='variant',values='HL').reset_index()
        alls = pd.DataFrame({'s': list(set(df['s']))})
        data_HL = pd.merge(pivoted_vars, alls, on='s', how='right')
        backbone = data_HL[['s']].copy()
        for vari in variants_to_extract:
            this_name = re.sub('[:,]','_', vari)
            backbone.loc[:,this_name] = data_HL[vari]

        backbone.to_csv(f'{TEMP}pivoted_flat_file_lowHLqc.tsv', sep='\t', index=False)

        ht = hl.import_table(f'{TEMP}pivoted_flat_file_lowHLqc.tsv', impute=True, missing='')
        ht = ht.annotate(s = hl.str(ht.s)).key_by('s')
        ht = ht.repartition(15).checkpoint(get_final_case_only_hl_path(num_to_keep, 'ht'))
    
    return ht


def get_positive_control_phenotype_names():
    measure_list = {3036277: 'centimeter', 
                    3004501: 'milligram per deciliter', 
                    903115: 'millimeter mercury column',
                    903118: 'millimeter mercury column'}
    measures = '903133, 3004501, 903115, 903118'
    allowed_operators = ['=','No matching concept']
    return measure_list, measures, allowed_operators


def grab_positive_control_phenos(measures, measure_list, allowed_operators):
    data_sql = f"""
        SELECT 
            measurement.PERSON_ID,
            measurement.MEASUREMENT_DATETIME,
            m_ext.src_id as DATA_SOURCE,
            measurement.visit_occurrence_id as VISIT_ID,
            measurement.MEASUREMENT_CONCEPT_ID, 
            m_standard_concept.concept_name as STANDARD_CONCEPT_NAME, 
            m_standard_concept.vocabulary_id as STANDARD_VOCABULARY, 
            measurement.UNIT_CONCEPT_ID, 
            m_unit.concept_name as UNIT_CONCEPT_NAME, 
            measurement.MEASUREMENT_TYPE_CONCEPT_ID, 
            m_type.concept_name as MEASUREMENT_TYPE_CONCEPT_NAME, 
            measurement.MEASUREMENT_SOURCE_CONCEPT_ID, 
            m_source_concept.concept_name as SOURCE_CONCEPT_NAME, 
            m_source_concept.vocabulary_id as SOURCE_VOCABULARY,
            measurement.VALUE_AS_NUMBER,
            measurement.VALUE_AS_CONCEPT_ID, 
            m_value.concept_name as VALUE_AS_CONCEPT_NAME, 
            m_operator.concept_name as OPERATOR_CONCEPT_NAME
        FROM `{DATASET}.measurement` measurement 
        LEFT JOIN `{DATASET}.measurement_ext` m_ext 
            ON measurement.measurement_id = m_ext.measurement_id
        LEFT JOIN `{DATASET}.concept` m_standard_concept 
            ON measurement.measurement_concept_id = m_standard_concept.concept_id 
        LEFT JOIN `{DATASET}.concept` m_unit 
            ON measurement.unit_concept_id = m_unit.concept_id 
        LEFT JOIN `{DATASET}.concept` m_type 
            ON measurement.measurement_type_concept_id = m_type.concept_id 
        LEFT JOIN `{DATASET}.concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id 
        LEFT JOIN `{DATASET}.concept` m_value 
            ON measurement.value_as_concept_id = m_value.concept_id 
        LEFT JOIN `{DATASET}.concept` m_operator 
            ON measurement.operator_concept_id = m_operator.concept_id 
        WHERE
            measurement_source_concept_id IN ({measures})"""
    
    wgs_data_measure = pd.read_gbq(data_sql, dialect="standard")
    wgs_data_measure_latest = wgs_data_measure.sort_values('MEASUREMENT_DATETIME'
                                                ).groupby(['MEASUREMENT_CONCEPT_ID', 'PERSON_ID']
                                                ).tail(1)
    wgs_data_measure_latest['correct_measure'] = wgs_data_measure_latest.MEASUREMENT_CONCEPT_ID.map(measure_list)
    wgs_data_measure_latest = wgs_data_measure_latest[wgs_data_measure_latest.correct_measure == wgs_data_measure_latest.UNIT_CONCEPT_NAME]
    wgs_data_measure_latest = wgs_data_measure_latest[wgs_data_measure_latest.OPERATOR_CONCEPT_NAME.isin(allowed_operators)]
    wgs_data_measure_latest['measurement_year'] = pd.DatetimeIndex(wgs_data_measure_latest.MEASUREMENT_DATETIME).year
    wgs_data_measure_latest['measurement_time'] = wgs_data_measure_latest.MEASUREMENT_DATETIME.astype('int') / 10**9

    return wgs_data_measure, wgs_data_measure_latest


def import_and_process_pos_control_phenos(path_latest_csv):
    ht_latest = hl.import_table(path_latest_csv, min_partitions=100, impute=True, 
                                types={'VALUE_AS_NUMBER': hl.tfloat64}, missing=['NA',''])
    ht_latest = ht_latest.annotate(phenotype = munge_phenotype_name(ht_latest.STANDARD_CONCEPT_NAME))
    ht_latest = ht_latest.annotate(s = hl.str(ht_latest.PERSON_ID)).drop('PERSON_ID').key_by('s')
    return ht_latest


def get_positive_control_phenotypes(sample_covariates, overwrite=False):

    measure_list, measures, allowed_operators = get_positive_control_phenotype_names()

    if hl.hadoop_exists(f'{SANITY_PHENOS_PATH}/_SUCCESS') and not overwrite:
        ht_pheno_out = hl.read_table(SANITY_PHENOS_PATH)
        ht_covar_out = hl.read_table(SANITY_PHENOS_COVARS_PATH)
        phenotypes = [x for x in ht_pheno_out.row if x not in ht_pheno_out.key]
    else:
        temp_phenos = f'{GWAS_DIR}sanity_phenotypes/demog_obtained_traits.tsv'
        temp_phenos_latest = f'{GWAS_DIR}sanity_phenotypes/demog_obtained_traits_latest.tsv'
        
        if not hl.hadoop_exists(temp_phenos) or overwrite:
            wgs_data_measure, wgs_data_measure_latest = grab_positive_control_phenos(measure_list=measure_list, measures=measures, allowed_operators=allowed_operators)
            wgs_data_measure.to_csv(temp_phenos, index=False, sep='\t')
            wgs_data_measure_latest.to_csv(temp_phenos_latest, index=False, sep='\t')

        #ht_all = hl.import_table(temp_phenos, min_partitions=100, impute=True)
        ht_latest = import_and_process_pos_control_phenos(temp_phenos_latest)
        ht_latest = ht_latest.annotate(year_of_birth = sample_covariates[ht_latest.s].year_of_birth,
                                       isFemale = sample_covariates[ht_latest.s].isFemale)
        ht_latest = ht_latest.annotate(approx_age = ht_latest.measurement_year - ht_latest.year_of_birth)
        ht_latest = ht_latest.filter(ht_latest.approx_age > 0)
        ht_latest = ht_latest.checkpoint(f'{TEMP}sanity_check_alldata_checkpoint.ht', overwrite=True)
        ht_pheno_out = ht_latest.select().distinct()
        ht_covar_out = ht_latest.select().distinct()
        phenotypes = list(ht_latest.aggregate(hl.agg.collect_as_set(ht_latest.phenotype)))

        for pheno in phenotypes:
            htf = ht_latest.filter(ht_latest.phenotype == pheno)
            ht_pheno_out = ht_pheno_out.annotate(**{pheno: htf[ht_pheno_out.s].VALUE_AS_NUMBER})

            struct_item = hl.Struct(**{'measurement_time': htf[ht_covar_out.s].measurement_time,
                                       'collection_site': htf[ht_covar_out.s].DATA_SOURCE,
                                       'visit_id': htf[ht_covar_out.s].VISIT_ID,
                                       'approx_age': htf[ht_covar_out.s].approx_age,
                                       'age2': htf[ht_covar_out.s].approx_age**2,
                                       'age_isFemale':htf[ht_covar_out.s].isFemale * htf[ht_covar_out.s].approx_age,
                                       'age2_isFemale': (htf[ht_covar_out.s].approx_age**2)*htf[ht_covar_out.s].isFemale})
            ht_covar_out = ht_covar_out.annotate(**{pheno: struct_item})
            ht_covar_out = ht_covar_out.checkpoint(f'{TEMP}sanity_check_covars_after_{pheno}_checkpoint.ht', overwrite=True)

        ht_pheno_out = ht_pheno_out.checkpoint(SANITY_PHENOS_PATH, overwrite=overwrite)
        ht_covar_out = ht_covar_out.checkpoint(SANITY_PHENOS_COVARS_PATH, overwrite=overwrite)

    return ht_pheno_out, ht_covar_out, phenotypes


def apply_irnt(ht, cols):
    # similar to Neale lab round 2 approach:
    # https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2/irnt.biomarkers.py
    
    df = ht.select(*cols).to_pandas()
    df.index = df['s']
    df = df.drop('s', axis=1)

    dfp = df.rank()
    dfp = (dfp - 0.5) / (~dfp.isnull()).sum()
    dfp.columns = [x + '_prob' for x in dfp.columns]
    dfp.loc[:, 's'] = dfp.index

    ht = hl.Table.from_pandas(dfp, key='s')
    ht = ht.annotate(**{x.replace('_prob', ''): hl.qnorm(ht[x])
                        for x in ht.row_value if x.endswith('_prob')})
    ht = ht.annotate(**{x: hl.or_missing(~hl.is_nan(ht[x]), ht[x]) for x in ht.row_value})
    ht = ht.drop(*[x for x in ht.row_value if x.endswith('_prob')])
    
    return ht


def get_genotypes(AF_cutoff=0.001, n_partitions=5000, overwrite=False, use_global_hwe=True):
    intermediate_path = f'{GWAS_DIR}filtered_gt_2.mt'
    final_path = f'{GWAS_DIR}filtered_gt_3.mt'

    if hl.hadoop_exists(f'{final_path}/_SUCCESS') and not overwrite:
        mt = hl.read_matrix_table(final_path)

    else:
        if hl.hadoop_exists(f'{intermediate_path}/_SUCCESS') and not overwrite:
            mt = hl.read_matrix_table(intermediate_path)
        else:
            mt = hl.read_matrix_table(os.getenv("WGS_HAIL_STORAGE_PATH"))
            mt = mt.select_entries('GT')

            mt = hl.split_multi_hts(mt, permit_shuffle=False)
            #bi = mt.filter_rows(hl.len(mt.alleles) == 2)
            #bi = bi.annotate_rows(a_index=1, was_split=False)
            #multi = mt.filter_rows(hl.len(mt.alleles) > 2)
            #split = hl.split_multi_hts(multi, permit_shuffle=False)
            #mt = split.union_rows(bi)

            mt = mt.annotate_rows(MAF = mt.info.AF[mt.a_index-1])
            mt = mt.annotate_rows(MAF = hl.min([mt.MAF, 1-mt.MAF]))
            mt = mt.filter_rows(mt.MAF > AF_cutoff, keep = True)
            mt = mt.naive_coalesce(n_partitions*2).checkpoint(f'{GWAS_DIR}filtered_gt_1.mt', overwrite=True)

            related_remove = hl.import_table(REL_SAMP_PATH,
                                             types={"sample_id.s":"tstr"},
                                             key="sample_id.s", min_partitions=20)
            mt = mt.anti_join_cols(related_remove)

            flagged_s_remove = hl.import_table(FLAGGED_SAMP_PATH,
                                               types={"s":"tstr"},
                                               key="s", min_partitions=20)
            mt = mt.anti_join_cols(flagged_s_remove)
            
            mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt.filters) == 0) | (mt.filters == {'PASS'}))
            mt = mt.checkpoint(intermediate_path, overwrite=True)

        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(hwe = hl.agg.hardy_weinberg_test(mt.GT, one_sided=True))
        if use_global_hwe:
            mt = mt.filter_rows(mt.hwe.p_value > 1e-10, keep=True)
            # mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-10, keep = True)
        mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95, keep = True)
        mt = mt.drop('variant_qc', 'hwe').repartition(n_partitions).checkpoint(final_path, overwrite=True)
    
    return mt


def run_regressions(mt, phenos, covars, pass_through, gwas_name, thresh=0, model='additive', overwrite=False, n_partition=1000):
    """
    Much of the heavy lifiting for converting the ht to mt is pulled from UKB round 2.
    Performs covariate correction for all elements in the covariates struct.
    """
    mt_dir = f'{RESULTS_DIR}mt/'
    cutoff = '' if thresh == 0 else f'_filt_prob_{thresh}'
    filename = f'{gwas_name}_gwas_{model}{cutoff}.mt'
    filename_raw = f'{gwas_name}_gwas_{model}{cutoff}_raw.mt'
    if (not overwrite) & hl.hadoop_exists(mt_dir+filename):
        mt = hl.read_matrix_table(mt_dir+filename)

    else:
        if model == 'additive':
            entry = mt.GT.n_alt_alleles()
        elif model == 'recessive':
            entry = hl.if_else(mt.GT.n_alt_alleles() == 2, 1, 0)

        ht_dir = TEMP + f'{gwas_name}_gwas_{model}{cutoff}_res.ht'
        if hl.hadoop_exists(ht_dir) and not overwrite:
            ht_out = hl.read_table(ht_dir)
        else:
            ht_out = hl.linear_regression_rows(
                y=[[mt['phenotypes'][y]] for y in phenos],
                x=entry,
                covariates=[1, *[mt['covariates'][item] for item in covars]],
                pass_through=pass_through)
            ht_out = ht_out.repartition(n_partition)
            ht_out = ht_out.checkpoint(ht_dir, overwrite=True)

        ht_res = ht_out.annotate_globals(columns=hl.map(lambda i: hl.struct(phenotype=hl.literal(phenos)[i]), 
                                                        hl.range(0, hl.len(phenos))))
        ht_res = ht_res.annotate(entries=hl.map(
            lambda i: hl.struct(
                n=ht_res['n'][i],
                sum_x=ht_res['sum_x'][i],
                y_transpose_x=ht_res['y_transpose_x'][i][0],
                beta=ht_res['beta'][i][0],
                standard_error=ht_res['standard_error'][i][0],
                t_stat=ht_res['t_stat'][i][0],
                p_value=ht_res['p_value'][i][0]),
            hl.range(0, hl.len(phenos))))
        ht_res = ht_res.select(*(pass_through + ['entries']))
        mt = ht_res._unlocalize_entries('entries', 'columns', ['phenotype'])
        mt = mt.repartition(n_partition)
        mt = mt.checkpoint(mt_dir+filename_raw, overwrite=True)

        mt = mt.select_entries(N = mt.n,
                               AC = mt.sum_x,
                               ytx = mt.y_transpose_x,
                               BETA = mt.beta,
                               SE = mt.standard_error,
                               tstat = mt.t_stat,
                               Pvalue = mt.p_value)
        mt = mt.annotate_entries(AF = mt.AC / (2 * mt.N))
        mt = mt.annotate_entries(minor_AF = hl.cond(mt.AF <= 0.5, mt.AF, 1.0-mt.AF),
                                 minor_AC = hl.cond(mt.AF <= 0.5, mt.AC, (2 * mt.N)-mt.AC))
        mt = mt.annotate_entries(low_confidence = mt.minor_AC <= 20)

        mt = mt.checkpoint(mt_dir+filename, overwrite=True)
    
    return mt


def export_for_manhattan(mt, phenos, entry_keep, model, fold, suffix, overwrite, include_cols_for_mung):
    if type(entry_keep) == str:
        entry_keep = [entry_keep]

    for pheno in phenos:
        file_out = f'{GWAS_DIR}/sumstats/{model}/{fold}/{pheno}{suffix}'
        if include_cols_for_mung:
            extra_cols = ['rsid']
        else:
            extra_cols = ['rsid']
        extra_cols = [x for x in extra_cols if x in mt.row]
        if overwrite or not hl.hadoop_exists(file_out):
            ht_f = mt.filter_cols(mt.phenotype == pheno).entries().select(*(entry_keep + extra_cols)).repartition(100)
            ht_f = ht_f.key_by()
            if include_cols_for_mung:
                ht_f = ht_f.rename({'minor_AF':'MAF', 'rsid':'SNP', 'n':'N', 
                                    'p_value': 'P', 'beta':'BETA', 
                                    'standard_error':'BETA_SE'})
                ht_f = ht_f.annotate(A1 = ht_f.alleles[1], A2 = ht_f.alleles[0])
                ht_f = ht_f.drop(ht_f.locus, ht_f.alleles, ht_f.phenotype).key_by('SNP')
            else:
                ht_f = ht_f.annotate(variant = hl.str(ht_f.locus)+ ':' + hl.delimit(ht_f.alleles, ':'))
                ht_f = ht_f.drop(ht_f.locus, ht_f.alleles, ht_f.phenotype).key_by('variant')
            ht_f.export(file_out)


def aou_generate_final_lambdas(mt, suffix, overwrite):
    mt = mt.annotate_cols(
        pheno_data=hl.zip(mt.pheno_data, hl.agg.array_agg(
            lambda ss: hl.agg.filter(~ss.low_confidence,
                hl.struct(lambda_gc=hl.methods.statgen._lambda_gc_agg(ss.Pvalue),
                          n_variants=hl.agg.count_where(hl.is_defined(ss.Pvalue)),
                          n_sig_variants=hl.agg.count_where(ss.Pvalue < 5e-8))),
            mt.summary_stats)).map(lambda x: x[0].annotate(**x[1]))
    )
    ht = mt.cols()
    ht = ht.checkpoint(get_lambdas_path(suffix, 'full', 'ht'), overwrite=overwrite, _read_if_exists=not overwrite)
    ht.explode('pheno_data').flatten().export(get_lambdas_path(suffix, 'full', 'txt.bgz'))
    return mt


def get_snv_count_phenotype():
    ht = hl.import_table(FINAL_VARIANT_CALLSET_PATH,
                        types={'locus': hl.tstr, 'alleles': hl.tstr, 's': hl.tstr,
                            'rsid': hl.tstr, 'variant':hl.tstr, 'batch':hl.tstr,
                            'common_low_heteroplasmy': hl.tbool,
                            'hap_defining_variant': hl.tbool,
                            'region': hl.tstr,
                            'variant_context': hl.tstr,
                            'dp_mean': hl.tfloat64,
                            'mq_mean': hl.tfloat64, 'tlod_mean': hl.tfloat64,
                            'AF_hom': hl.tfloat64, 'AF_het': hl.tfloat64,
                            'AC_hom': hl.tint64, 'AC_het': hl.tint64,
                            'max_hl': hl.tfloat64, 'dp_mean': hl.tfloat64,
                            'major_haplogroup':hl.tstr, 'hap':hl.tstr, 
                            'wgs_median_coverage':hl.tfloat64, 'mt_mean_coverage':hl.tfloat64,
                            'mito_cn':hl.tfloat64, 'age':hl.tfloat64, 'pop':hl.tstr,
                            'DP':hl.tfloat64, 'AD_ref':hl.tfloat64,'AD_alt':hl.tfloat64,
                            'MQ':hl.tfloat64,'TLOD':hl.tfloat64, 'GT':hl.tstr, 'HL':hl.tfloat64,
                            'OriginalSelfRefAlleles':hl.tstr, 'SwappedFieldIDs':hl.tstr, 'FT':hl.tstr, 'FT_LIFT':hl.tstr,
                            'artifact_prone':hl.tbool, 'lifted':hl.tbool, 'fail_gt':hl.tbool, 'missing_call':hl.tbool}, min_partitions=50)
    all_s_table = ht.select('s').key_by('s').distinct()
    ht_heteroplasmies = ht.filter(hl.is_defined(ht.HL) & (ht.HL < 0.95)).repartition(40)
    ht_heteroplasmies = ht_heteroplasmies.annotate(alleles = ht_heteroplasmies.alleles.split(','))
    ht_heteroplasmies = ht_heteroplasmies.filter(~hl.is_indel(ht_heteroplasmies.alleles[0], ht_heteroplasmies.alleles[1]))
    ht_snv_count = ht_heteroplasmies.group_by(ht_heteroplasmies.s).aggregate(N = hl.agg.count())
    s_not_found = all_s_table.anti_join(ht_snv_count).annotate(N = hl.int64(hl.literal(0)))
    ht_snv_count = hl.Table.union(ht_snv_count, s_not_found)
    return ht_snv_count.rename({'N':'snv_count_qcpass'})


def filter_mt_per_pop_maf(mt, pop, cutoff, overwrite_gt, perform_per_pop_hwe=False):
    """ Expects that MT has a GT and .covariates.pop field.
    Also filters based on p_hwe (two sided) per-population.
    """
    this_pop_path = f'{TEMP}mt/genotype_mt_filtered_{pop}_maf_pass.mt'
    if hl.hadoop_exists(f'{this_pop_path}/_SUCCESS') and not overwrite_gt:
        mt_af_filt = hl.read_matrix_table(this_pop_path)
    else:
        mt_af_filt = mt.filter_cols(mt.covariates.pop == pop)
        
        if perform_per_pop_hwe:
            print('Running per-pop HWE filtering!')
            mt_af_filt = mt_af_filt.annotate_rows(hwe = hl.agg.hardy_weinberg_test(mt_af_filt.GT))
            mt_af_filt = mt_af_filt.filter_rows(mt_af_filt.hwe.p_value > 1e-10, keep = True).drop('hwe')
        
        mt_af_filt = mt_af_filt.annotate_rows(gt_stats = hl.agg.call_stats(mt_af_filt.GT, mt_af_filt.alleles))
        mt_af_filt = mt_af_filt.annotate_rows(minor_AF = hl.min(mt_af_filt.gt_stats.AF))
        mt_af_filt = mt_af_filt.filter_rows(mt_af_filt.minor_AF > cutoff).drop('gt_stats')
        mt_af_filt = mt_af_filt.checkpoint(this_pop_path, overwrite=True)
    return mt_af_filt