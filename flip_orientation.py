"""
author: @rahulg

This simple function loads a VCF, performs basic checks, and swaps reference sequence.
The input VCF should be the original input into `bcftools consensus` -- i.e. the sites that were
  flipped in the new sequence. The purpose of this script is to produce a new VCF in
  target sequence coordinates with variants (relative to reference) that need to be force-called.
"""

import hail as hl
import argparse


def fai_to_len(fai):
    with open(fai) as f:
        line = f.readline()
    return int(line.split('\t')[1])


def check_vcf_integrity(mt):
    # check that locus, alleles are the two key fields
    if sorted(list(mt.row_key)) != ['alleles', 'locus']:
        raise ValueError('VCFs must always be keyed by locus, alleles.')
    
    # check that all sites are bi-allelic
    if mt.aggregate_rows(~hl.agg.all(hl.len(mt.alleles) == 2)):
        raise ValueError('This function only supports biallelic sites (run SplitMultiAllelics!)')
    
    # check that there is no missingness in locus
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.locus))):
        raise ValueError('ERROR: locus must always be defined, both before and after Liftover. This should be a reversible operation, thus finding missing loci after reverse liftover is very concerning.')

    # check that there is no missingness in alleles
    if mt.aggregate_rows(~hl.agg.all(hl.map(hl.is_defined, mt.alleles))):
        raise ValueError('ERROR: alleles should always be defined.')


parser = argparse.ArgumentParser()
parser.add_argument('--target-name', required=True)
parser.add_argument('--target-fasta', required=True)
parser.add_argument('--target-fai', required=True)
parser.add_argument('--source-name', required=True)
parser.add_argument('--source-fasta', required=True)
parser.add_argument('--source-fai', required=True)
parser.add_argument('--chain', required=True)
parser.add_argument('--source-reference-vcf', required=True)
parser.add_argument('--output', required=True)


if __name__ == "__main__":
    args = parser.parse_args()

    target = hl.ReferenceGenome(args.target_name, ['chrM'], {'chrM':fai_to_len(args.target_fai)}, mt_contigs=['chrM'])
    source = hl.ReferenceGenome(args.source_name, ['chrM'], {'chrM':fai_to_len(args.source_fai)}, mt_contigs=['chrM'])
    target.add_sequence(args.target_fasta, args.target_fai)
    source.add_sequence(args.source_fasta, args.source_fai)
    source.add_liftover(args.chain, args.target_name)

    mt_new = hl.import_vcf(args.source_reference_vcf, reference_genome=args.source_name)
    check_vcf_integrity(mt_new)
    mt_new = mt_new.select_rows().select_entries()
    mt_new = mt_new.annotate_rows(new_locus = hl.liftover(mt_new.locus, args.target_name))
    mt_new = mt_new.annotate_rows(allele_flip = hl.reversed(mt_new.alleles))
    mt_new = mt_new.key_rows_by().rename({'locus':'locus_orig', 'alleles':'alleles_orig'}).rename({'new_locus':'locus', 'allele_flip': 'alleles'}).key_rows_by('locus','alleles')
    mt_new = mt_new.drop('locus_orig', 'alleles_orig')
    check_vcf_integrity(mt_new)

    hl.export_vcf(mt_new, args.output, tabix=True)