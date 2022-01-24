version 1.0

#import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates
#import "https://api.firecloud.org/ga4gh/v1/tools/mitochondria:AlignmentPipeline/versions/1/plain-WDL/descriptor" as AlignAndMarkDuplicates
import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/AlignmentPipeline.wdl" as AlignAndMarkDuplicates
import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/AlignAndCall.wdl" as AlignAndCallFull

workflow AlignAndCall {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File unmapped_bam
    String base_name

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    #Shifted reference is used for calling the control region (edge of mitochondria reference).
    #This solves the problem that BWA doesn't support alignment to circular contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    File shift_back_chain

    # Optional arguments override hardcoded definitions of the (shifted) control region and non-control region
    # Intervals are built using the first row from the following interval files in M2:
    File? non_control_interval
    File? control_shifted

    File? gatk_override
    String? gatk_docker_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Boolean compress_output_vcf

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    #Optional runtime arguments
    Int? preemptible_tries
  }

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallFull.CollectWgsMetrics as CollectWgsMetrics{
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  Int M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7
  Boolean defined_custom_noncntrl = defined(non_control_interval)
  String noncntrl_args_suffix = if defined_custom_noncntrl then "" else " -L chrM:576-16024 "

  call AlignAndCallFull.M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + noncntrl_args_suffix,
      custom_interval = non_control_interval,
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  Boolean defined_custom_cntrl = defined(control_shifted)
  String cntrl_args_suffix = if defined_custom_cntrl then "" else " -L chrM:8025-9144 "

  call AlignAndCallFull.M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fai = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      # Only the control region is now called.
      m2_extra_args = select_first([m2_extra_args, ""]) + cntrl_args_suffix,
      custom_interval = control_shifted,
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallFull.LiftoverAndCombineVcfs as LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallShiftedMt.raw_vcf,
      vcf = CallMt.raw_vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain,
      preemptible_tries = preemptible_tries
  }

  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = LiftoverAndCombineVcfs.merged_vcf
    File out_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index
  }
}