from gnomad_methods.gnomad.sample_qc.ancestry import pc_project
from .aou_paths import *
import hail as hl
import os

# globals
GNOMAD_VAR_PATH = f'{PCA_DIR}/gnomad_hq_variants_without_washu_preldprune.ht'
POPS = ('eur', 'amr', 'afr', 'eas', 'sas') # dropping MID as sample size is too small

# parameters
n_partitions = 2000
min_maf_hq = 0.001
ld_r2 = 0.1
iteration = 1
k = 20
global_overwrite = True

# initialize
hl.init(default_reference="GRCh38", tmp_dir=TEMP)


def get_custom_pc_path(iteration):
    return f'{COVARIATES_DIR}recomputed_pca_gnomad{make_iteration_suffix(iteration)}.ht'


def make_iteration_suffix(iteration):
    if iteration == 0:
        return ''
    else:
        return f'_iter{str(iteration)}'


def remove_ancestry_outliers(mt: hl.MatrixTable, iteration):
    """ Filters ancestry outliers according to iteration number.
    """
    if iteration == 0:
        return mt
    else:
        ht_list = [hl.read_table(f'{FILTERED_SAMPLES_DIR}manual_ancestry_outliers{make_iteration_suffix(x)}.ht') for x in range(1, iteration+1)]
        ht_for_removal = hl.Table.union(*ht_list).distinct()
        return mt.anti_join_cols(ht_for_removal)


def run_pca(mt: hl.MatrixTable, out_prefix: str, k=20, overwrite: bool = False):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=k, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = add_pcs(pca_scores, k)
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', overwrite)  # PCA loadings


def project_individuals(pca_loadings, project_mt, k):
    """
    Project samples into predefined PCA space
    :param pca_loadings: existing PCA space
    :param project_mt: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = add_pcs(ht_projections, k)
    return ht_projections


def import_ld_pruned_hq_variants(mt, overwrite=False):
    """ Perform per-population LD pruning of the original high quality variant set from gnomAD.
    """
    final_pruned_ht_path = f'{PCA_DIR}gnomad_hq_variants_without_washu_perpop_pruned.ht'
    if hl.hadoop_exists(f'{final_pruned_ht_path}/_SUCCESS') and not overwrite:
        ht = hl.read_table(final_pruned_ht_path)

    else:
        # munge gnomAD pre-pruned hq variants
        pre_pruning_ht_path = f'{PCA_DIR}gnomad_hq_variants_without_washu_prepruning.ht'
        if hl.hadoop_exists(f'{pre_pruning_ht_path}/_SUCCESS') and not overwrite:
            ht_pre = hl.read_table(pre_pruning_ht_path)
            
        else:
            ht_pre = hl.import_table(f'{PCA_DIR}resources/pre_ld_pruning_combined_variants_without_washu_gnomad.tsv.bgz', 
                                     types={'locus':hl.tstr, 'alleles':hl.tarray(hl.tstr)}, impute=True)
            ht_pre = ht_pre.annotate(locus = hl.parse_locus(ht_pre.locus, reference_genome='GRCh38')).key_by('locus','alleles')
            ht_pre = ht_pre.repartition(100).checkpoint(pre_pruning_ht_path, overwrite=overwrite)

        # filter the genotype mt to those variants found in the HQ variant table
        ld_clump_gt_path = f'{PCA_DIR}/filtered_gt_for_LD_clump_gnomadqc.mt'
        if hl.hadoop_exists(f'{ld_clump_gt_path}/_SUCCESS') and not overwrite:
            mt_pre_clump = hl.read_matrix_table(ld_clump_gt_path)

        else:
            mt_pre_clump = mt.semi_join_rows(ht_pre)
            mt_pre_clump = mt_pre_clump.naive_coalesce(n_partitions*2).checkpoint(ld_clump_gt_path, overwrite=overwrite)
            # this saved successfully with 258457 rows, a subset of the starting 259482 and 98k samples.
        
        ht_post_clump = mt_pre_clump.rows()
        ht_post_clump = ht_post_clump.annotate_globals(pop = hl.empty_array(hl.tstr), starting_vars = hl.empty_array(hl.tint32), samples = hl.empty_array(hl.tint32), ld_r2 = hl.empty_array(hl.tfloat64))
        ht_post_clump = ht_post_clump.annotate(found_in_clumped_set = hl.empty_array(hl.tbool))
        
        for pop in POPS:
            this_pruned_var_path = f'{PCA_DIR}/pruning/ld_pruned_results_maf{str(min_maf_hq)}_gnomadqc_{pop}.ht'
            if hl.hadoop_exists(f'{this_pruned_var_path}/_SUCCESS') and not overwrite:
                ht_pruned = hl.read_table(this_pruned_var_path)

            else:
                # keep only variants with MAF > min MAF within the pop of interest
                mtf = mt_pre_clump.filter_cols(mt_pre_clump.pop == pop)
                mtf = mtf.annotate_rows(call_stats = hl.agg.call_stats(mtf.GT, mtf.alleles))
                mtf = mtf.filter_rows(hl.min(mtf.call_stats.AF) >= min_maf_hq, keep=True)
                this_shape = mtf.count()

                # perform LD pruning
                ht_pruned = hl.ld_prune(mtf.GT, r2=ld_r2)
                ht_pruned = ht_pruned.annotate_globals(pop=pop, starting_vars=this_shape[0], samples=this_shape[1], ld_r2=ld_r2)
                ht_pruned = ht_pruned.checkpoint(this_pruned_var_path, overwrite=overwrite)
            
            # append true if a variant was retained after filtering and pruning per population; else false
            ht_post_clump = ht_post_clump.annotate_globals(**{x: ht_post_clump[x].append(ht_pruned[x].collect()[0]) for x in ['pop','starting_vars','samples','ld_r2']})
            ht_post_clump = ht_post_clump.annotate(found_in_clumped_set = ht_post_clump.found_in_clumped_set.append(hl.is_defined(ht_pruned[ht_post_clump.key])))

        ht = ht_post_clump.filter(hl.any(ht_post_clump.found_in_clumped_set)) # has 137394 variants that passed in >1 ancestry
        ht = ht.checkpoint(final_pruned_ht_path, overwrite=overwrite)
    
    return ht


def import_aou_gt_for_pca(overwrite=False):
    """ Generates full genotype file for PCA. If MT is not found, will
    draw this file from the raw genotype data.
    """
    final_mt_path = f'{PCA_DIR}final_filtered_gt_for_PCA_gnomadqc.mt'
    
    if not overwrite and hl.hadoop_exists(f'{final_mt_path}/_SUCCESS'):
        mt = hl.read_matrix_table(final_mt_path)
    
    else:
        prelim_filt_mt = f'{TEMP}preliminary_filt_gt_for_PCA.mt'
        if hl.hadoop_exists(f'{prelim_filt_mt}/_SUCCESS') and not overwrite:
            mt = hl.read_matrix_table(prelim_filt_mt)
        else:
            # in step 1, split multi allelics, remove variants with 0 MAF, remove variants which don't pass
            mt = hl.read_matrix_table(os.getenv("WGS_HAIL_STORAGE_PATH"))
            mt = mt.select_entries('GT')
            mt = hl.split_multi_hts(mt, permit_shuffle=False)
            mt = mt.annotate_rows(AF = mt.info.AF[mt.a_index-1])
            mt = mt.filter_rows(hl.min([mt.AF, 1-mt.AF]) > 0, keep = True)
            mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt.filters) == 0) | (mt.filters == {'PASS'}))
            mt = mt.naive_coalesce(n_partitions*10).checkpoint(prelim_filt_mt, overwrite=overwrite)

        # filter full MT and add ancestry information
        ancestry_pred = hl.import_table(ANCESTRY_INFO_PATH,
                                        key="research_id", 
                                        impute=True, 
                                        types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)},
                                        min_partitions=50)
        related_remove = hl.import_table(REL_SAMP_PATH,
                                        types={"sample_id.s":"tstr"},
                                        key="sample_id.s", min_partitions=20)
        flagged_s_remove = hl.import_table(FLAGGED_SAMP_PATH,
                                            types={"s":"tstr"},
                                            key="s", min_partitions=20)
        mt = mt.anti_join_cols(flagged_s_remove) # remove flagged samples
        # NOTE bug here -- ~is_defined means that a sample is tagged as related if it is NOT found in the related list
        mt = mt.annotate_cols(pop = ancestry_pred[mt.col_key].ancestry_pred, 
                              related = ~hl.is_defined(related_remove[mt.col_key]))

        # import hq site data
        # we do not use the AoU provided list because this includes only chromosome 20 and 21
        # we check if the LD-pruned set of variants has been produced; if not we do LD-pruning
        ht_hq = import_ld_pruned_hq_variants(mt, overwrite=overwrite)

        mt = mt.semi_join_rows(ht_hq) # filter to HQ variants and add per-ancestry inclusion variables
        mt = mt.annotate_globals(pop_index=ht_hq.pop.collect()[0])
        mt = mt.annotate_rows(include_in_pca_by_anc = ht_hq[mt.row_key].found_in_clumped_set)
        mt = mt.naive_coalesce(n_partitions).checkpoint(final_mt_path, overwrite=True)

    return mt


def add_pcs(ht, k, goal=20):
    """ If ht contains fewer PCs than k, adds columns with 0s to get up to k.
    """
    ht = ht.annotate(**{f'PC{i}': ht.scores[i - 1] for i in range(1, k+1)})
    if k < goal:
        ht = ht.annotate(**{f'PC{i}': 0 for i in range(1, goal+1) if f'PC{i}' not in ht.row})
    ht = ht.select(*[f'PC{i}' for i in range(1, goal+1)])
    return ht


def filter_to_clumped_variants(mt, pop):
    popidx = hl.eval(mt.pop_index.index(pop))
    if popidx is None:
        raise ValueError('pop must be found in the globals of the imported mt!')
    return mt.filter_rows(mt.include_in_pca_by_anc[popidx])


if __name__ == '__main__':
    # import full MT
    # this function will produce this PCA set if it isn't found
    mt = import_aou_gt_for_pca(overwrite=False)
    mt = mt.annotate_cols(related = ~mt.related) # band aid relatedness fix
    mt = remove_ancestry_outliers(mt, iteration)
    iter_suff = make_iteration_suffix(iteration)

    # for each population, remove related samples, run PCA, and then project unrelated samples onto PCs
    hts = []
    for pop in POPS:
        related_prefix = f'{PCA_DIR}{pop}_results_related_gnomadqc{iter_suff}'
        unrelated_prefix = f'{PCA_DIR}{pop}_results_unrelated_gnomadqc{iter_suff}'

        file_found = hl.hadoop_exists(f'{related_prefix}.scores_projected.ht/_SUCCESS') and \
            hl.hadoop_exists(f'{unrelated_prefix}.scores.ht/_SUCCESS')
        
        if global_overwrite or not file_found:
            mtf = mt.filter_cols(mt.pop == pop)
            mtf = filter_to_clumped_variants(mtf, pop)
            mtf_unrel = mtf.filter_cols(~mtf.related)
            mtf_rel = mtf.filter_cols(mtf.related)
            print(f'For pop {pop}, {str(mtf_unrel.count_cols())} unrelated samples found with {str(mtf_unrel.count_rows())} variants.')
            run_pca(mtf_unrel, f'{unrelated_prefix}.', k, True)

            pca_loadings = hl.read_table(f'{unrelated_prefix}.loadings.ht')
            ht = project_individuals(pca_loadings, mtf_rel, k)
            ht.write(f'{related_prefix}.scores_projected.ht', overwrite=True)
            hl.read_table(f'{related_prefix}.scores_projected.ht').export(
                f'{related_prefix}.scores_projected.txt.bgz')

        ht_score_rel = hl.read_table(f'{related_prefix}.scores_projected.ht')
        hts.append(ht_score_rel.annotate(pop=pop, related=True))
        ht_score_unrel = hl.read_table(f'{unrelated_prefix}.scores.ht')
        ht_score_unrel = add_pcs(ht_score_unrel, k)
        hts.append(ht_score_unrel.annotate(pop=pop, related=False))

    # merged table now contains PCs, pop, and relatedness information
    ht = hts[0].union(*hts[1:])
    ht.write(get_custom_pc_path(iteration), overwrite=True)

    # filtering iteration 0
    # EUR: 34576 variants
    # AMR: 48101 variants
    # AFR: 65611 variants
    # EAS: 27044 variants
    # SAS: 46434 variants

    # filtering iteration 1
    # EUR: 34576 variants
    # AMR: 48101 variants
    # AFR: 65611 variants
    # EAS: 27037 variants
    # SAS: 46346 variants