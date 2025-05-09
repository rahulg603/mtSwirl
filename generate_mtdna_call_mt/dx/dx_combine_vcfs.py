import argparse
from curses import pair_content
import logging
import math
import os, sys
import dxpy
import pyspark

import hail as hl

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from merging_utils import collect_vcf_paths, vcf_merging_and_processing


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("combine_mitochondria_vcfs_into_mt")
logger.setLevel(logging.INFO)


def main(args):  # noqa: D103
    # start SQL session and initialize constants
    my_database = dxpy.find_one_data_object(name=args.dx_init.lower())["id"]
    participant_data = f'dnax://{my_database}/{args.input_ht}'
    coverage_mt_path = f'dnax://{my_database}/{args.coverage_mt_path}/'
    output_bucket = f'dnax://{my_database}/{args.output_bucket}'
    temp_dir = f'dnax://{my_database}/{args.temp_dir}/'
    participants_to_subset = None if args.participants_to_subset is None else f'dnax://{my_database}/{args.participants_to_subset}'
    chunk_size = args.chunk_size
    artifact_prone_sites_path = args.artifact_prone_sites_path
    artifact_prone_sites_reference = args.artifact_prone_sites_reference
    include_extra_v2_fields = args.include_extra_v2_fields
    file_name = args.file_name
    minimum_homref_coverage = args.minimum_homref_coverage
    num_merges = args.split_merging
    vcf_col_name = args.vcf_col_name
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    hl.init(sc=sc, tmp_dir=temp_dir)

    if int(hl.version().split('-')[0].split('.')[2]) >= 75: # only use this if using hail 0.2.75 or greater
        logger.info("Setting hail flag to avoid array index out of bounds error...")
        # Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
        # TODO: reassess if this flag is still needed for future versions of Hail
        hl._set_flags(no_whole_stage_codegen="1")

    logger.info("Collecting VCF paths for samples to subset...")
    vcf_paths = collect_vcf_paths(
        participant_data, vcf_col_name, participants_to_subset
    )
    vcf_paths = {k: f'file://{v}' for k, v in vcf_paths.items()}

    if args.append_to_existing is not None:
        logger.info('Will search for existing mt to join new data into.')
        if args.existing_dx_init is not None:
            existing_database = dxpy.find_one_data_object(name=args.existing_dx_init.lower())["id"]
        else:
            existing_database = my_database
        old_mt_path = f'dnax://{existing_database}/{args.append_to_existing}/'
    else:
        old_mt_path = None

    combined_mt, meta = vcf_merging_and_processing(vcf_paths=vcf_paths, 
                                                   coverage_mt_path=coverage_mt_path,
                                                   include_extra_v2_fields=include_extra_v2_fields, 
                                                   single_sample=False,
                                                   old_mt_path=old_mt_path,
                                                   artifact_prone_sites_path=artifact_prone_sites_path,
                                                   artifact_prone_sites_reference=artifact_prone_sites_reference,
                                                   minimum_homref_coverage=minimum_homref_coverage,
                                                   logger=logger, 
                                                   chunk_size=chunk_size, 
                                                   num_merges=num_merges,
                                                   n_final_partitions=args.n_final_partitions,
                                                   output_bucket=output_bucket,
                                                   temp_dir=temp_dir,
                                                   overwrite=args.overwrite)

    logger.info("Writing combined MT...")

    # Set the file names for output files
    out_vcf = f"file://{os.getcwd()}/{file_name}.vcf.bgz"
    out_mt = f"{output_bucket}/{file_name}.mt"
    out_tsv = f"file://{os.getcwd()}/{file_name}.tsv.bgz"

    combined_mt = combined_mt.repartition(args.n_final_partitions).checkpoint(out_mt, overwrite=args.overwrite, _read_if_exists=not args.overwrite)

    logger.info("Writing trimmed variants table...")
    logger.info("We export missing variants; all others are considered homozygous reference.")
    ht_for_tsv = combined_mt.entries()
    ht_for_tsv = ht_for_tsv.filter(hl.is_missing(ht_for_tsv.HL) | (ht_for_tsv.HL > 0))
    ht_for_tsv.repartition(300).export(out_tsv)

    logger.info("Writing combined VCF...")
    # For the VCF output, join FT values by semicolon
    combined_mt = combined_mt.annotate_entries(
        FT=hl.str(";").join(hl.array(combined_mt.FT))
    )
    hl.export_vcf(combined_mt.repartition(300), out_vcf, metadata=meta)



if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="This script combines individual mitochondria VCF files into one MatrixTable, determines homoplasmic reference sites, and applies an artifact_prone_site filter"
    )
    p.add_argument(
        "-i",
        "--input-ht",
        help="Input ht with paths to vcf file.",
        required=True,
    )
    p.add_argument(
        "-c",
        "--coverage-mt-path",
        help="Path to MatrixTable of sample-level coverage (per-sample and per-base, can be generated by running annotate_coverage.py)",
        required=True,
    )
    p.add_argument(
        "-v",
        "--vcf-col-name",
        help="Name of column in participant data file that contains the path to the VCF output by Mutect2",
        required=True,
    )
    p.add_argument(
        "-a",
        "--artifact-prone-sites-path",
        help="Path to BED file of artifact-prone sites to flag in the FILTER column",
        required=True,
    )
    p.add_argument(
        "-a-ref",
        "--artifact-prone-sites-reference",
        help="Reference genome for artifact-prone sites BED file. If specified can be GRCh37 or GRCh38",
        default='default',
        type=str,
    )
    p.add_argument(
        "-o",
        "--output-bucket",
        help="Path to bucket to which results should be written",
        required=True,
    )
    p.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    p.add_argument(
        "-f",
        "--file-name",
        help="File name to use for output files (will be used the the .vcf.bgz and .mt outputs)",
        required=True,
    )
    p.add_argument(
        "-s",
        "--participants-to-subset",
        help="Path to txt file of participant_ids to which the data should be subset (file should contain header (named 'participant') and one line for each participant_id matching the 'entity:participant_id's supplied in Terra",
    )
    p.add_argument(
        "--minimum-homref-coverage",
        help="Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing",
        type=int,
        default=100,
    )
    p.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=5,
    )
    p.add_argument("--overwrite", help="Overwrites existing files", action="store_true")
    p.add_argument("--include-extra-v2-fields", help="Loads in exåtra fields from MitochondriaPipeline v2.1, specifically AD, OriginalSelfRefAlleles, and SwappedFieldIDs. If missing will fill with missing.", action="store_true")
    p.add_argument(
        "--dx-init", type=str, required=True, help='SQL database path for use in DNAnexus.'
    )
    p.add_argument(
        "--n-final-partitions", type=int, default=1000, help='Number of partitions for final mt.'
    )
    p.add_argument(
        '--split-merging', type=int, default=1, help='Will split the merging into this many jobs which will be merged at the end. Uses the same order each time such that if it fails we can read from previous files.'
    )
    p.add_argument(
        '--append-to-existing', type=str, help='Optional: specify relative path to existing coverage MatrixTable to merge this dataset into.'
    )
    p.add_argument(
        '--existing-dx-init', type=str, help='Optional: specify dx database for existing MT. If not provided, will assume it is --dx-init.'
    )

    args = p.parse_args()

    main(args)
