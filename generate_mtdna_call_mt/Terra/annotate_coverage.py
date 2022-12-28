import argparse
import logging
import math
import os
import re
import sys

import hail as hl

from os.path import dirname
from gnomad.utils.slack import slack_notifications
from hail.utils.java import info

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Annotate coverage")
logger.setLevel(logging.INFO)


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


# def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int, min_partitions: int) -> hl.MatrixTable:
#     """
#     Hierarchically join together MatrixTables in the provided list.

#     :param mts: List of MatrixTables to join together
#     :param temp_dir: Path to temporary directory for intermediate results
#     :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
#     :return: Joined MatrixTable
#     """
#     # Convert the MatrixTables to tables where entries are an array of structs
#     staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
#     stage = 0
#     while len(staging) > 1:
#         # Calculate the number of jobs to run based on the chunk size
#         n_jobs = int(math.ceil(len(staging) / chunk_size))
#         info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
#         next_stage = []

#         for i in range(n_jobs):
#             # Grab just the tables for the given job
#             to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
#             info(
#                 f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
#             )

#             # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
#             merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
#             if min_partitions > 4:
#                 merged = merged.checkpoint(f"stage_{stage}_job_{i}_pre.ht", overwrite=True)
#             # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
#             merged = merged.annotate(
#                 __entries=hl.flatten(
#                     hl.range(hl.len(merged.__entries)).map(
#                         # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
#                         lambda i: hl.coalesce(
#                             merged.__entries[i].__entries,
#                             hl.range(hl.len(merged.__cols[i].__cols)).map(
#                                 lambda j: hl.null(
#                                     merged.__entries.__entries.dtype.element_type.element_type
#                                 )
#                             ),
#                         )
#                     )
#                 )
#             )

#             # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
#             merged = merged.annotate_globals(
#                 __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
#             )

#             next_stage.append(
#                 merged.checkpoint(
#                     os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
#                 )
#             )
#         info(f"Completed stage {stage}")
#         stage += 1
#         staging.clear()
#         staging.extend(next_stage)

#     # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
#     return (
#         staging[0]
#         ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
#         .unfilter_entries()
#     )


def chunks(items, binsize):
    lst = []
    for item in items:
        lst.append(item)
        if len(lst) == binsize:
            yield lst
            lst = []
    if len(lst) > 0:
        yield lst


def main(args):  # noqa: D103
    input_tsv = args.input_tsv
    output_ht = args.output_ht
    temp_dir = args.temp_dir
    chunk_size = args.chunk_size
    overwrite = args.overwrite
    keep_targets = args.keep_targets
    num_merges = args.split_merging
    hl.init(tmp_dir=temp_dir)
    hl._set_flags(no_whole_stage_codegen='1')

    if args.overwrite == False and hl.hadoop_exists(output_ht):
        logger.warning(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_ht,
        )
    # Ensure that user supplied ht extension for output_ht
    if not output_ht.endswith(".ht"):
        sys.exit("Path supplied as output_ht must end with .ht extension")

    logger.info(
        "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
    )

    paths = hl.import_table(input_tsv)
    pairs_for_coverage = paths.annotate(pairs = (paths.s, paths.coverage)).pairs.collect()

    if num_merges > 1:
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
                    for s, base_level_coverage_metrics in subset:
                        idx+=1
                        mt = hl.import_matrix_table(
                            base_level_coverage_metrics,
                            delimiter="\t",
                            row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                            row_key=["chrom", "pos"],
                            min_partitions=args.n_read_partitions,
                        )
                        if not keep_targets:
                            mt = mt.drop("target")
                        else:
                            mt = mt.key_rows_by(*["chrom", "pos", "target"])
                        mt = mt.key_cols_by().annotate_cols(col_id = s)
                        mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')

                        mt_list.append(mt)
                        if idx % 10 == 0:
                            logger.info(f"Imported batch {str(idx)}, subset {str(subset_number)}...")

                    logger.info(f"Joining individual coverage mts for subset {str(subset_number)}...")
                    cov_mt_this = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=args.n_read_partitions, check_from_disk=False, prefix=this_prefix)
                    cov_mt_this = cov_mt_this.repartition(args.n_final_partitions // num_merges).checkpoint(this_subset_mt, overwrite=True)
                    mt_list_subsets.append(cov_mt_this)
            cov_mt = multi_way_union_mts(mt_list_subsets, temp_dir, chunk_size, min_partitions=args.n_read_partitions, check_from_disk=False, prefix=merged_prefix)
            cov_mt = cov_mt.repartition(args.n_final_partitions).checkpoint(this_merged_mt, overwrite=True)
    else:
        mt_list = []
        idx = 0
        
        for s, base_level_coverage_metrics in pairs_for_coverage:
            idx+=1
            mt = hl.import_matrix_table(
                base_level_coverage_metrics,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
                min_partitions=args.n_read_partitions,
            )
            if not keep_targets:
                mt = mt.drop("target")
            else:
                mt = mt.key_rows_by(*["chrom", "pos", "target"])
            mt = mt.key_cols_by().annotate_cols(col_id = s)
            mt = mt.rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')

            mt_list.append(mt)
            if idx % 10 == 0:
                logger.info(f"Imported batch {str(idx)}...")

        logger.info("Joining individual coverage mts...")
        cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=args.n_read_partitions, check_from_disk=False, prefix='')
    
    n_samples = cov_mt.count_cols()

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_rows(
        locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )
    cov_mt.show()

    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

    output_mt = re.sub(r"\.ht$", ".mt", output_ht)
    output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
    output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

    if not args.hail_only:
        logger.info("Writing sample level coverage...")
        sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
        sample_mt.coverage.export(output_samples)

    logger.info("Writing coverage mt and ht...")
    cov_mt.repartition(args.n_final_partitions).write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.rows()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    if not args.hail_only:
        cov_ht.export(output_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "-i",
        "--input-tsv",
        help="Input file with coverage files to combine in tab-delimited format of participant_id, base_level_coverage_metrics, sample",
        required=True,
    )
    parser.add_argument(
        "-o", "--output-ht", help="Name of ht to write output", required=True
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    parser.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--keep-targets", help="Will add an annotation for target from the coverage file", action="store_true"
    )
    parser.add_argument(
        "--n-read-partitions", type=int, help="The number of partitions to use when reading tsvs. This should be 1 if the files are small.", default=1
    )
    parser.add_argument(
        "--hail-only", action='store_true', help='Skip generating flat files.'
    )    
    parser.add_argument(
        "--n-final-partitions", type=int, default=1000, help='Number of partitions for final mt.'
    )    
    parser.add_argument(
        '--split-merging', type=int, default=1, help='Will split the merging into this many jobs which will be merged at the end. Uses the same order each time such that if it fails we can read from previous files.'
    )

    args = parser.parse_args()
    main(args)
