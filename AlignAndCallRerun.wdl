version 1.0

#import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates
#import "https://api.firecloud.org/ga4gh/v1/tools/mitochondria:AlignmentPipeline/versions/1/plain-WDL/descriptor" as AlignAndMarkDuplicates
import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/AlignmentPipeline.wdl" as AlignAndMarkDuplicates
import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/AlignAndCall.wdl" as AlignAndCallFull

workflow AlignAndCallRerun {
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
    File blacklisted_sites
    File blacklisted_sites_index

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
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf

    Float? verifyBamID
    File? force_call_vcf
    File? force_call_vcf_index
    File? force_call_vcf_shifted
    File? force_call_vcf_shifted_index
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    #Optional runtime arguments
    Int? preemptible_tries
    Int? n_cpu
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
      preemptible_tries = preemptible_tries,
      force_call_vcf = force_call_vcf,
      force_call_vcf_index = force_call_vcf_index,
      n_cpu = n_cpu
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
      preemptible_tries = preemptible_tries,
      force_call_vcf = force_call_vcf_shifted,
      force_call_vcf_index = force_call_vcf_shifted_index,
      n_cpu = n_cpu
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

  call AlignAndCallFull.MergeStats as MergeStats {
    input:
      shifted_stats = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallFull.Filter as FilterContamination {
    input:
      raw_vcf = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      run_contamination = true,
      hasContamination = hasContamination,
      contamination_major = contamination_major,
      contamination_minor = contamination_minor,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      preemptible_tries = preemptible_tries
 }

  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_idx
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    Float median_coverage = CollectWgsMetrics.median_coverage
  }
}