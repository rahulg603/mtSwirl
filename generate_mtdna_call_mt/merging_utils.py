import math
import os
import hail as hl
import logging

from hail.utils.java import info


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int, min_partitions: int, check_from_disk: bool, prefix: str) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.

    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    if check_from_disk:
        staging = [x for x in mts]
    else:
        staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        
        next_stage = []
        
        if check_from_disk:
            all_exists = True
            for idx in range(n_jobs):
                path = os.path.join(temp_dir, f"stage_{stage}_job_{idx}.ht")
                exists = hl.hadoop_is_file(f'{path}/_SUCCESS')
                if not exists:
                    print(path + ' is missing.')
                    if stage == 0:
                        raise ValueError('ERROR: --check-from-disk was enabled but not all stage 0 MTs were found. This is unsupported.')
                    all_exists = False
                    break
            
            if all_exists:
                info(f"Reading stage {stage} from disk...")
                staging.clear()
                for idx in range(n_jobs):
                    staging.append(hl.read_table(os.path.join(temp_dir, f"stage_{stage}_job_{idx}.ht")))
                info(f"Stage {stage} imported from disk.")
                stage += 1
                continue

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            if min_partitions > 10:
                merged = merged.checkpoint(f"{prefix}stage_{stage}_job_{i}_pre.ht", overwrite=True)
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"{prefix}stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def chunks(items, binsize):
    lst = []
    for item in items:
        lst.append(item)
        if len(lst) == binsize:
            yield lst
            lst = []
    if len(lst) > 0:
        yield lst


def coverage_merging(paths, num_merges, chunk_size, check_from_disk, 
                     temp_dir, n_read_partitions, n_final_partitions, 
                     keep_targets, logger, no_batch_mode=False):
    
    if no_batch_mode:
        pairs_for_coverage = paths.annotate(pairs = (paths.s, paths.coverage)).pairs.collect()
    else:
        pairs_for_coverage = paths.annotate(pairs = (paths.batch, paths.coverage)).pairs.collect()
    
    if num_merges > 1:
        # check_from_disk is not compatible with multiple merges and will not be used
        merged_prefix = f'coverage_merging_final_{str(num_merges)}subsets/'
        this_merged_mt = os.path.join(temp_dir, f"{merged_prefix}final_merged.mt")
        if hl.hadoop_is_file(this_merged_mt + '/_SUCCESS'):
            cov_mt = hl.read_matrix_table(this_merged_mt)
        else:
            subsets = chunks(pairs_for_coverage, len(pairs_for_coverage) // num_merges)
            mt_list_subsets = []
            for subset_number, subset in enumerate(subsets):
                print(f'Importing subset {str(subset_number)}...')
                this_prefix = f'coverage_merging_subset{str(subset_number)}_{str(num_merges)}subsets/'
                this_subset_mt = os.path.join(temp_dir, f"{this_prefix}final_merged.mt")
                if hl.hadoop_is_file(f'{this_subset_mt}/_SUCCESS'):
                    mt_list_subsets.append(hl.read_matrix_table(this_subset_mt))
                    print(f'Subset {str(subset_number)} already processed and imported with {str(mt_list_subsets[len(mt_list_subsets)-1].count_cols())} samples.')
                else:
                    mt_list = []
                    idx = 0
                    for batch, base_level_coverage_metrics in subset:
                        idx+=1
                        mt = hl.import_matrix_table(
                            base_level_coverage_metrics,
                            delimiter="\t",
                            row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                            row_key=["chrom", "pos"],
                            min_partitions=n_read_partitions,
                        )
                        if not keep_targets:
                            mt = mt.drop("target")
                        else:
                            mt = mt.key_rows_by(*["chrom", "pos", "target"])

                        if no_batch_mode:
                            mt = mt.key_cols_by().annotate_cols(col_id = batch)
                            mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                        else:
                            mt = mt.key_cols_by().rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                            mt = mt.annotate_cols(batch = batch)

                        mt_list.append(mt)
                        if idx % 10 == 0:
                            logger.info(f"Imported batch {str(idx)}, subset {str(subset_number)}...")

                    logger.info(f"Joining individual coverage mts for subset {str(subset_number)}...")
                    cov_mt_this = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=False, prefix=this_prefix)
                    cov_mt_this = cov_mt_this.repartition(n_final_partitions // num_merges).checkpoint(this_subset_mt, overwrite=True)
                    mt_list_subsets.append(cov_mt_this)
            cov_mt = multi_way_union_mts(mt_list_subsets, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=False, prefix=merged_prefix)
            cov_mt = cov_mt.repartition(n_final_partitions).checkpoint(this_merged_mt, overwrite=True)
    else:
        mt_list = []
        idx = 0
        if check_from_disk:
            logger.info("NOTE: Skipping reading individual coverage MTs since --check-from-disk was enabled.")
            n_append = len(pairs_for_coverage)-1
            pairs_for_coverage = pairs_for_coverage[0:1]
        
        for batch, base_level_coverage_metrics in pairs_for_coverage:
            idx+=1
            mt = hl.import_matrix_table(
                base_level_coverage_metrics,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
                min_partitions=n_read_partitions,
            )
            if not keep_targets:
                mt = mt.drop("target")
            else:
                mt = mt.key_rows_by(*["chrom", "pos", "target"])
            
            if no_batch_mode:
                mt = mt.key_cols_by().annotate_cols(col_id = batch)
                mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
            else:
                mt = mt.key_cols_by().rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
                mt = mt.annotate_cols(batch = batch)

            mt_list.append(mt)
            if idx % 10 == 0:
                logger.info(f"Imported batch {str(idx)}...")

        if check_from_disk:
            mt_list.extend([None for x in range(n_append)])

        logger.info("Joining individual coverage mts...")
        cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=n_read_partitions, check_from_disk=check_from_disk, prefix='')

    cov_mt = cov_mt.annotate_rows(locus = hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"))
    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")
    return cov_mt


def join_two_mts(mt1, mt2, row_keep, col_keep, temp_dir, partitions):

    mt1 = mt1.select_rows(*row_keep).select_cols(*col_keep)
    mt2 = mt2.select_rows(*row_keep).select_cols(*col_keep)

    if (mt1.row.dtype != mt2.row.dtype) or (mt1.col.dtype != mt2.col.dtype) or (mt1.entry.dtype != mt2.entry.dtype):
        raise ValueError('ERROR: when joining MatrixTables, schemas must be the same.')
    
    mt_append = multi_way_union_mts(mts=[mt1,mt2], temp_dir=temp_dir, chunk_size=2, min_partitions=partitions, check_from_disk=False, prefix='appending_')
    return mt_append


def append_coverage_to_old(cov_mt, old_mt_path, n_final_partitions, temp_dir):

    this_merged_mt = os.path.join(temp_dir, 'coverage_tmp_appended_to_old_dataset_final.mt')
    cov_mt = cov_mt.checkpoint(os.path.join(temp_dir, 'coverage_mt_new_keyed_pre_merge_with_old.mt'), overwrite=True)

    if hl.hadoop_is_file(this_merged_mt + '/_SUCCESS'):
        cov_mt = hl.read_matrix_table(this_merged_mt)
    else:
        old_mt = hl.read_matrix_table(old_mt_path)
        print(f'Second database imported with {str(old_mt.count()[1])} samples.')
        cov_mt = join_two_mts(mt1 = old_mt, mt2 = cov_mt, row_keep = [], col_keep = ['batch'], temp_dir=temp_dir, partitions=n_final_partitions)
        cov_mt = cov_mt.repartition(n_final_partitions).checkpoint(this_merged_mt, overwrite=True) 

    return cov_mt


def add_coverage_annotations(cov_mt):
    n_samples = cov_mt.count_cols()
    cov_mt = cov_mt.annotate_rows(
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )
    return cov_mt

