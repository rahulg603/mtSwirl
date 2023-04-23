import os, re
from tempfile import tempdir
import pandas as pd
import hail as hl
from gnomad_mitochondria.aou_gwas.aou_gwas_helpers import get_covariates, get_case_only_mtdna_callset, \
     filter_mt_per_pop_maf, apply_irnt, get_genotypes, run_regressions, export_for_manhattan, aou_generate_final_lambdas
from gnomad_mitochondria.aou_gwas.aou_paths import *


ANALYSIS_POP = ['afr','amr','eur','sas','eas']
MIN_CASES = 100
IRNT = True
OVERWRITE_SUMSTATS = True
EXPORT_SUMSTATS = True
overwrite_full_gt = True


def apply_qc_continuous(mt, min_case_count: int = 50):
    mt = mt.filter_cols(mt.n_cases >= min_case_count)
    return mt


def run_full_gwas(sample_covariates, mt, ht_pheno, num_PC, overwrite_gt, naming_insert, fold, pheno, min_cases):
    
    # Prepare mt for regressions
    mt_a = mt.annotate_cols(phenotypes = ht_pheno[mt.s])
    mt_a = mt_a.annotate_cols(covariates = sample_covariates[mt_a.s])
    covar_full_list = list(mt_a.covariates.keys())
    covars_base = ['isFemale','approx_age','age_isFemale','age2','age2_isFemale'] + \
        [f'PC{str(x)}'for x in range(1,num_PC+1)]
    covars = covars_base + [x for x in covar_full_list if re.search('^hap_', x)]
    print('Using covariates:')
    print(covars)
    irntsuff = '_irnt' if IRNT else ''
    this_suffix = lambda pop: f'221206_{pop}_{naming_insert}_mtdna_variant_qc_hl_case_only{irntsuff}'

    # Run per-population GWAS
    mts = []
    for pop in ANALYSIS_POP:

        # filter table by per-pop MAF (and to specific pop)
        print(f'For pop {pop}, started with {str(mt_a.count_rows())} rows.')
        mt_af_filt = filter_mt_per_pop_maf(mt_a, pop, 0.0005, overwrite_gt)
        print(f'After filtering, now have {str(mt_af_filt.count_rows())} rows.')
        pheno_f = [x for x in pheno if mt_af_filt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt_af_filt.phenotypes[x]))) > min_cases]


        # Run variant HL GWAS and export sumstats
        res = run_regressions(mt_af_filt, pheno_f, covars, ['rsid'], 
                              this_suffix(pop), overwrite=OVERWRITE_SUMSTATS)

        if EXPORT_SUMSTATS:
            export_for_manhattan(mt=res, phenos=pheno_f, entry_keep=['N','Pvalue','BETA','SE','tstat','ytx','AC','minor_AC','AF', 'minor_AF', 'low_confidence'], 
                                 model='additive', fold=fold, suffix=f'_{this_suffix(pop)}_geno_af_0.001.tsv.bgz', 
                                 overwrite=OVERWRITE_SUMSTATS, include_cols_for_mung=False)
        
        # finish formatting MT
        res = res.annotate_cols(n_cases = hl.array(hl.agg.collect_as_set(res.N)).filter(lambda x: hl.is_defined(x))[0],
                                n_controls = hl.missing(hl.tint32),
                                pop = pop,
                                inv_normalized = IRNT,
                                log_pvalue = False)
        res = apply_qc_continuous(res)
        res = res.select_cols(pheno_data=res.col_value)
        res = res.select_entries(summary_stats=res.entry)
        mts.append(res)
    
    # Collapse into single MT
    full_mt = mts[0]
    for this_mt in mts[1:]:
        full_mt = full_mt.union_cols(this_mt, row_join_type='outer')
    full_mt = full_mt.checkpoint(f'{TEMP}mt/staging_full.mt', overwrite=True)
    full_mt = full_mt.collect_cols_by_key()
    full_mt = full_mt.checkpoint(f'{TEMP}mt/staging_lambdas.mt', overwrite=True)
    full_mt = aou_generate_final_lambdas(full_mt, this_suffix('full'), overwrite=True)
    full_mt.write(f'{RESULTS_DIR}/all_pop_mt/{this_suffix("full")}.mt', overwrite=True)


# Import covariates
sample_covariates_iter0 = get_covariates(0, overwrite=False, use_custom_pca=True)

# Import phenotypes and IRNT
ht_pheno = get_case_only_mtdna_callset(num_to_keep=300, overwrite=False)
pheno = [x for x in ht_pheno.row if x not in ht_pheno.key]
if IRNT:
    ht_pheno = apply_irnt(ht_pheno, pheno)
    ht_pheno = ht_pheno.checkpoint(f'{TEMP}phenotypes_after_irnt_checkpoint.ht', overwrite=True)
    ht_pheno.export(f'{BUCKET}/analyses/221206_final_gwas_heteroplasmies_hwefix/Data/phenotypes_post_irnt.tsv')

# Grab genotypes
mt = get_genotypes(0.0001, 10000, overwrite_full_gt, use_global_hwe=True).drop('MAF')

# Run GWAS using new PCs, no iterations (raw)
fold = '221206_final_gwas_heteroplasmies_hwefix'
run_full_gwas(sample_covariates_iter0, mt, ht_pheno, num_PC=20, overwrite_gt=True, 
              naming_insert='newPCs_iter0_final', fold=fold, pheno=pheno, min_cases=MIN_CASES)

# AFR: 40832109
# AMR: 34642693
# EUR: 22513978
# SAS: 27266217
# EAS: 17856139

