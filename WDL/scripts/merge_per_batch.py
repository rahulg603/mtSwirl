import argparse
import math
import os

import hail as hl


META_DICT = {
    "filter": {
        "artifact_prone_site": {
            "Description": "Variant overlaps an artifact-prone site"
        }
    },
    "format": {
        "DP": {"Description": "Depth of coverage", "Number": "1", "Type": "Integer"},
        "FT": {
            "Description": "Sample-level genotype filters",
            "Number": ".",
            "Type": "String",
        },
        "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
        "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
        "TLOD": {
            "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
            "Number": "1",
            "Type": "Float",
        },
        "MPOS": {
            "Description": "median distance from end of read",
            "Number": "1",
            "Type": "Integer",
        },
        "AS_SB_TABLE": {
            "Description": "Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.",
            "Number": "1",
            "Type": "String",
        },
        "STR": {
            "Description": "Variant is a short tandem repeat. 1 if True, 0 if False.",
            "Number": "1",
            "Type": "Integer",
        },
        "STRQ": {
            "Description": "Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors",
            "Number": "1",
            "Type": "Integer",
        },
        "RPA": {
            "Description": "Number of times tandem repeat unit is repeated, for each allele (including reference)",
            "Number": "R",
            "Type": "Integer",
        },
        'AD': {"Description": "Allelic depth of REF and ALT", "Number": "R", "Type": "Integer"},
        'OriginalSelfRefAlleles': {
            'Description':'Original self-reference alleles (only if alleles were changed in Liftover repair pipeline)', 
            'Number':'R', 
            'Type':'String'},
        'SwappedFieldIDs': {
            'Description':'Fields remapped during liftover (only if alleles were changed in Liftover repair pipeline)', 
            'Number':'1', 
            'Type':'String'
        },
        'F2R1': {"Description": "Count of reads in F2R1 pair orientation supporting each allele", "Number": "R", "Type": "Integer"},
        'F1R2': {"Description": "Count of reads in F1R2 pair orientation supporting each allele", "Number": "R", "Type": "Integer"}
    },
}


def read_input_data(path):
    vcf_paths = {}

    ht = hl.import_table(path)
    ht = ht.filter(ht.path != "")
    
    # Add the vcf path to a dictionary with sample name as key
    df = ht.to_pandas()

    for _, row in df.iterrows():
        vcf_paths[row["s"]] = row['path']

    return vcf_paths


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int, min_partitions: int) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.

    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        next_stage = []

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            if min_partitions > 10:
                merged = merged.checkpoint(f"stage_{stage}_job_{i}_pre.ht", overwrite=True)
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
                    os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def run_coverage(args):
    input_paths = read_input_data(args.input_tsv)
    mt_list = []
    for sample, file_path in input_paths.items():
        mt = hl.import_matrix_table(
            file_path,
            delimiter="\t",
            row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
            row_key=["chrom", "pos"],
            min_partitions=args.n_read_partitions,
        )
        mt = mt.key_rows_by(*["chrom", "pos", "target"])
        mt = mt.rename({"x": "coverage"})
        mt = mt.key_cols_by(s=sample)
        mt_list.append(mt)

    cov_mt = multi_way_union_mts(mt_list, args.temp_dir, args.chunk_size, min_partitions=args.n_read_partitions)
    
    cov_mt.describe()
    cov_mt.show()
    
    cov_mt.coverage.export(args.output_flat_file)


def run_variants(args):
    input_paths = read_input_data(args.input_tsv)
    mt_list = []
    for sample, file_path in input_paths.items():
        mt = hl.import_vcf(file_path, reference_genome="GRCh38")
        fields_of_interest = {'OriginalSelfRefAlleles':'array<str>', 'SwappedFieldIDs':'str',
                              'F2R1':'array<int32>', 'F1R2':'array<int32>'}
        if 'GT' in mt.entry:
            mt = mt.drop('GT')
        for x, item_type in fields_of_interest.items():
            if x not in mt.entry:
                mt = mt.annotate_entries(**{x: hl.missing(item_type)})
        mt = mt.select_entries("DP", "AD", *list(fields_of_interest.keys()), HL=mt.AF[0])
        # tuples are if we should index to the 0th element of an array
        info_of_interest = {'MPOS':('int32',True), 'AS_SB_TABLE':('str',False), 
                            'STR':('int32',False), 'STRQ':('int32',False), 'RPA':('array<int32>',False)}
        mt = mt.annotate_entries(
            MQ=hl.float(mt.info["MMQ"][1]),
            TLOD=mt.info["TLOD"][0],
            FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters),
        )
        for x, (item_type, to_index) in info_of_interest.items():
            if x not in mt.info:
                mt = mt.annotate_entries(**{x: hl.missing(item_type)})
            else:
                if to_index:
                    mt = mt.annotate_entries(**{x: mt.info[x][0]})
                else:
                    mt = mt.annotate_entries(**{x: mt.info[x]})
                if mt.info[x].dtype == hl.dtype('tbool'):
                    # flag is not supported in FORMAT
                    mt = mt.annotate_entries(**{x: hl.if_else(mt.info[x], 1, 0)})
        mt = mt.key_cols_by(s=sample)
        mt = mt.select_rows()
        mt_list.append(mt)

    cov_mt = multi_way_union_mts(mt_list, args.temp_dir, args.chunk_size, min_partitions=args.n_read_partitions)

    cov_mt.describe()
    cov_mt.show()
    cov_mt.entries().show()

    hl.export_vcf(cov_mt, args.output_flat_file, metadata=META_DICT)


parser = argparse.ArgumentParser()
parser.add_argument('--run-coverage', action='store_true', 
                    help='If enabled, expects to merge coverage tsvs.')
parser.add_argument('--run-variants', action='store_true', 
                    help='If true, expects to merge variants. Superceded by --run-coverage.')
parser.add_argument("--n-read-partitions", type=int, 
                    help="The number of partitions to use when reading tsvs. This should be 1 if the files are small.", default=1)
parser.add_argument("--chunk-size", 
                    help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
                    type=int,default=100)
parser.add_argument("--input-tsv", 
                    help="Table with a 's' column for sample name and 'path' column for files to import.", required=True)
parser.add_argument("--output-flat-file", 
                    help="Name of flat file to write output", required=True)
parser.add_argument("--temp-dir", 
                    help="Temporary directory to use for intermediate outputs", required=True)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.run_coverage:
        run_coverage(args)
    else:
        run_variants(args)
