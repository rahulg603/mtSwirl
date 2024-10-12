import argparse
import logging
import re
import sys, os
import dxpy
import pyspark

import hail as hl

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from merging_utils import coverage_merging, append_coverage_to_old, add_coverage_annotations

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Annotate coverage")
logger.setLevel(logging.INFO)


def main(args):  # noqa: D103
    # start SQL session and initialize constants
    my_database = dxpy.find_one_data_object(name=args.dx_init.lower())["id"]
    input_ht = f'dnax://{my_database}/{args.input_ht}/'
    output_ht = f'dnax://{my_database}/{args.output_ht}'
    temp_dir = f'dnax://{my_database}/{args.temp_dir}/'
    chunk_size = args.chunk_size
    overwrite = args.overwrite
    keep_targets = args.keep_targets
    check_from_disk = args.check_from_disk
    num_merges = args.split_merging
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    hl.init(sc=sc, tmp_dir=temp_dir)
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
    paths = hl.read_table(input_ht)

    cov_mt = coverage_merging(paths=paths, num_merges=num_merges, chunk_size=chunk_size, check_from_disk=check_from_disk, 
                              n_read_partitions=args.n_read_partitions, n_final_partitions=args.n_final_partitions, 
                              keep_targets=keep_targets, logger=logger, temp_dir=temp_dir)
    
    if args.append_to_existing is not None:
        logger.info('Appending coverage table to existing MT.')

        if args.existing_dx_init is not None:
            existing_database = dxpy.find_one_data_object(name=args.existing_dx_init.lower())["id"]
        else:
            existing_database = my_database
        
        old_mt_path = f'dnax://{existing_database}/{args.append_to_existing}/'

        cov_mt = append_coverage_to_old(cov_mt, old_mt_path, col_keep=['batch'],
                                        n_final_partitions=args.n_final_partitions, temp_dir=temp_dir)
        logger.info('Coverage table successfully appended.')

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = add_coverage_annotations(cov_mt)
    cov_mt.show()

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
        "-i", "--input-ht", help="Input ht with paths to coverage file.", required=True,
    )
    parser.add_argument(
        "-o", "--output-ht", help="Name of ht to write output", required=True
    )
    parser.add_argument(
        "-t", "--temp-dir", help="Temporary directory to use for intermediate outputs", required=True,
    )
    parser.add_argument(
        "--chunk-size", type=int, default=5,
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
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
        "--check-from-disk", action='store_true', help='If enabled and if the correct number of temporary files are found for a given stage, will read from those files.'
    )
    parser.add_argument(
        "--dx-init", type=str, required=True, help='SQL database path for use in DNAnexus.'
    )
    parser.add_argument(
        "--n-final-partitions", type=int, default=1000, help='Number of partitions for final mt.'
    )
    parser.add_argument(
        '--split-merging', type=int, default=1, help='Will split the merging into this many jobs which will be merged at the end. Uses the same order each time such that if it fails we can read from previous files.'
    )
    parser.add_argument(
        '--append-to-existing', type=str, help='Optional: specify relative path to existing coverage MatrixTable to merge this dataset into.'
    )
    parser.add_argument(
        '--existing-dx-init', type=str, help='Optional: specify dx database for existing MT. If not provided, will assume it is --dx-init.'
    )
    
    args = parser.parse_args()
    main(args)
