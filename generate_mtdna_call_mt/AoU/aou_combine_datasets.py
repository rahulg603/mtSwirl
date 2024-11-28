import argparse
from curses import pair_content
import logging
import math
import re
import os, sys

import hail as hl

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from merging_utils import append_coverage_to_old, add_coverage_annotations, \
                          append_vcf_to_old, get_vcf_metadata, apply_mito_artifact_filter

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("combine_tables_into_mt")
logger.setLevel(logging.INFO)

def main(args):
    if not args.coverage and not args.variants:
        raise ValueError('ERROR: enable either coverage mode or variants mode.')
    if args.coverage and args.variants:
        raise ValueError('ERROR: cannot enable both coverage and variants mode.')

    hl.init(tmp_dir=args.temp_dir, log='combine_datasets.log')
    if int(hl.version().split('-')[0].split('.')[2]) >= 75: # only use this if using hail 0.2.75 or greater
        logger.info("Setting hail flag to avoid array index out of bounds error...")
        # Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
        # TODO: reassess if this flag is still needed for future versions of Hail
        hl._set_flags(no_whole_stage_codegen="1")

    mt1 = hl.read_matrix_table(args.t1)

    if args.coverage:
        mt = append_coverage_to_old(mt1, args.t2, col_keep=['batch'], 
                                    n_final_partitions=args.n_final_partitions,
                                    temp_dir=args.temp_dir)
        logger.info('Coverage table successfully appended.')
        logger.info("Adding coverage annotations...")
        # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
        cov_mt = add_coverage_annotations(mt)

        output_ht = re.sub(r"\.mt$", ".ht", args.output_mt)
        output_tsv = re.sub(r"\.mt$", ".tsv", args.output_mt)
        output_samples = re.sub(r"\.mt$", "_sample_level.txt", args.output_mt)

        logger.info("Writing coverage mt and ht...")
        cov_mt.repartition(args.n_final_partitions).write(args.output_mt, overwrite=args.overwrite)
        cov_ht = cov_mt.rows()
        cov_ht = cov_ht.checkpoint(output_ht, overwrite=args.overwrite)

        if args.output_coverage_flat:
            logger.info("Writing sample level coverage...")
            sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
            sample_mt.coverage.export(output_samples)
            cov_ht.export(output_tsv)


    if args.variants:
        logger.info("Appending new VCF to old VCF database...")
        col_keep = [] if args.single_sample_analysis else ['batch']
        combined_mt = append_vcf_to_old(mt1, args.t2, col_keep, 
                                        args.n_final_partitions, args.temp_dir)

        logger.info("Applying artifact_prone_site filter...")
        combined_mt = apply_mito_artifact_filter(combined_mt, args.artifact_prone_sites_path, args.artifact_prone_sites_reference)

        logger.info("Writing combined MT...")
        # Set the file names for output files
        output_vcf = re.sub(r"\.mt$", ".vcf.bgz", args.output_mt)
        output_tsv = re.sub(r"\.mt$", ".tsv.bgz", args.output_mt)

        combined_mt = combined_mt.repartition(args.n_final_partitions).checkpoint(args.output_mt, overwrite=args.overwrite)

        logger.info("Writing trimmed variants table...")
        ht_for_tsv = combined_mt.entries()
        ht_for_tsv = ht_for_tsv.filter(hl.is_missing(ht_for_tsv.HL) | (ht_for_tsv.HL > 0))
        ht_for_tsv.naive_coalesce(300).export(output_tsv)

        if args.output_vcf:
            logger.info("Writing combined VCF...")
            # For the VCF output, join FT values by semicolon
            combined_mt = combined_mt.annotate_entries(
                FT=hl.str(";").join(hl.array(combined_mt.FT))
            )
            meta = get_vcf_metadata(args.include_extra_v2_fields)
            hl.export_vcf(combined_mt.naive_coalesce(300), output_vcf, metadata=meta)


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="This script combines individual mitochondria VCF files into one MatrixTable, determines homoplasmic reference sites, and applies an artifact_prone_site filter"
    )
    p.add_argument("--t1", help="Path to first MT", required=True)
    p.add_argument("--t2", help="Path to second MT", required=True)
    p.add_argument("--coverage", action='store_true', help='Coverage merging mode')
    p.add_argument("--variants", action='store_true', help='VCF merging mode')
    p.add_argument("--overwrite", help="Overwrites existing files", action="store_true")
    p.add_argument(
        "-o",
        "--output-mt",
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
        "-a",
        "--artifact-prone-sites-path",
        help="Path to BED file of artifact-prone sites to flag in the FILTER column"
    )
    p.add_argument(
        "-a-ref",
        "--artifact-prone-sites-reference",
        help="Reference genome for artifact-prone sites BED file. If specified can be GRCh37 or GRCh38",
        default='default',
        type=str,
    )
    p.add_argument("--n-final-partitions", type=int, default=1000, help='Number of partitions for final mt.')
    p.add_argument('--output-vcf', type=str, help='If specified, will generate the VCF (size explodes with sample size).')
    p.add_argument('--single-sample-analysis', type=str, help='If specified, uses single-sample mode.')
    p.add_argument('--output-coverage-flat', type=str, help='If specified, will export per-base coverage flat file..')
    p.add_argument('--include-extra-v2-fields', help='Includes new v2 fields if enabled.', action='store_true')

    args = p.parse_args()

    main(args)