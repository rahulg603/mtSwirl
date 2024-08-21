version 1.0

import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/Parallel_MongoTasks_v2_6_Multi.wdl" as MongoTasks_Multi

workflow ParallelAlignAndCallR2 {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    Array[File] unmapped_bam
    Array[String] sample_name
    String suffix

    File selfref_bundle
    Array[String] mt_interval_list

    Array[String] mt_self
    Array[String] mt_self_index
    Array[String] mt_self_dict

    Array[String] self_cat
    Array[String] self_cat_index
    Array[String] self_cat_dict

    Array[String] mt_self_shifted
    Array[String] mt_self_shifted_index
    Array[String] mt_self_shifted_dict

    Array[String] self_shifted_cat
    Array[String] self_shifted_cat_index
    Array[String] self_shifted_cat_dict

    Array[File] blacklisted_sites
    Array[File] blacklisted_sites_index

    Array[String] force_call_vcf
    Array[String] force_call_vcf_idx
    Array[String] force_call_vcf_shifted
    Array[String] force_call_vcf_shifted_idx

    Array[String] shift_back_chain

    Array[String] non_control_interval
    Array[String] control_shifted

    Array[Float]? verifyBamID
    Array[String] hasContamination
    Array[Float] contamination_major
    Array[Float] contamination_minor

    File? gatk_override
    String? gatk_docker_override
    File JsonTools
    String gatk_version = "4.2.6.0"
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    #Optional runtime arguments
    Int? preemptible_tries
    Int? n_cpu
    Int? n_cpu_bwa
  }

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call MongoTasks_Multi.ParallelMongoAlignToMtRegShiftedAndMetrics as AlignToMtRegShiftedAndMetrics {
    input:
      input_bam = unmapped_bam,
      sample_base_name = sample_name,
      suffix = suffix,

      selfref_bundle = selfref_bundle,
      mt = mt_self,
      mt_index = mt_self_index,
      mt_dict = mt_self_dict,

      mt_cat = self_cat,
      mt_cat_index = self_cat_index,
      mt_cat_dict = self_cat_dict,

      mt_shifted = mt_self_shifted,
      mt_shifted_index = mt_self_shifted_index,
      mt_shifted_dict = mt_self_shifted_dict,

      mt_shifted_cat = self_shifted_cat,
      mt_shifted_cat_index = self_shifted_cat_index,
      mt_shifted_cat_dict = self_shifted_cat_dict,
      
      mt_interval_list = mt_interval_list,

      read_length = max_read_length,
      coverage_cap = 100000,
      
      preemptible_tries = preemptible_tries,
      JsonTools = JsonTools,
      n_cpu = 16
  }

  Int M2_mem = if AlignToMtRegShiftedAndMetrics.max_mean_coverage > 25000 then 14 else 7

  call MongoTasks_Multi.MongoCallMtAndShifted as CallMtAndShifted {
    input:
      sample_base_name = AlignToMtRegShiftedAndMetrics.samples,
      suffix = suffix,
      selfref_bundle = selfref_bundle,
      # Everything is called except the control region.
      input_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam,
      input_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai,
      mt_self = mt_self,
      mt_self_index = mt_self_index,
      mt_self_dict = mt_self_dict,
      mt_interval_list = non_control_interval,
      force_call_vcf = force_call_vcf,
      force_call_vcf_idx = force_call_vcf_idx,

      # Only the control region is now called.
      shifted_input_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam,
      #shifted_input_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai,
      shifted_mt_self = mt_self_shifted,
      shifted_mt_self_index = mt_self_shifted_index,
      shifted_mt_self_dict = mt_self_shifted_dict,
      shifted_mt_interval_list = control_shifted,
      shifted_force_call_vcf = force_call_vcf_shifted,
      shifted_force_call_vcf_idx = force_call_vcf_shifted_idx,

      m2_extra_args = select_first([m2_extra_args, ""]),
      shifted_m2_extra_args = select_first([m2_extra_args, ""]),

      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      mem = M2_mem,
      preemptible_tries = preemptible_tries,
      JsonTools = JsonTools,
      n_cpu = n_cpu
  }

  call MongoTasks_Multi.MongoLiftoverCombineMergeFilterContamSplit as LiftoverCombineMergeFilterContamSplit {
    input:
      sample_base_name = sample_name,
      suffix = suffix,
      selfref_bundle = selfref_bundle,

      mt_self = mt_self,
      mt_self_index = mt_self_index,
      mt_self_dict = mt_self_dict,
      shifted_vcf = CallMtAndShifted.shifted_raw_vcf,
      shifted_vcf_idx = CallMtAndShifted.shifted_raw_vcf_idx,
      non_shifted_vcf = CallMtAndShifted.raw_vcf,
      non_shifted_vcf_idx = CallMtAndShifted.raw_vcf_idx,
      shifted_stats = CallMtAndShifted.shifted_stats,
      non_shifted_stats = CallMtAndShifted.stats,
      shift_back_chain = shift_back_chain, 
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,

      hasContamination = hasContamination,
      contamination_major = contamination_major,
      contamination_minor = contamination_minor,
      verifyBamID = verifyBamID,
      
      compress = compress_output_vcf,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries
  }

  output {
    Array[String] samples = LiftoverCombineMergeFilterContamSplit.samples
    Array[File] mt_aligned_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam
    Array[File] mt_aligned_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai
    Array[File] mt_aligned_shifted_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam
    #Array[File] mt_aligned_shifted_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai
    Array[File] nuc_mt_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bam
    #Array[File] nuc_mt_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bai
    Array[File] nuc_mt_shifted_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bam
    #Array[File] nuc_mt_shifted_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bai
    Array[File] out_vcf = LiftoverCombineMergeFilterContamSplit.filtered_vcf
    Array[File] out_vcf_idx = LiftoverCombineMergeFilterContamSplit.filtered_vcf_idx
    Array[File] split_vcf = LiftoverCombineMergeFilterContamSplit.split_vcf
    Array[File] split_vcf_idx = LiftoverCombineMergeFilterContamSplit.split_vcf_idx
    Array[File] duplicate_metrics = AlignToMtRegShiftedAndMetrics.duplicate_metrics
    Array[File] coverage_metrics = AlignToMtRegShiftedAndMetrics.wgs_metrics
    Array[File] theoretical_sensitivity_metrics = AlignToMtRegShiftedAndMetrics.theoretical_sensitivity
    Array[Int] mean_coverage = AlignToMtRegShiftedAndMetrics.mean_coverage
    Array[Float] median_coverage = AlignToMtRegShiftedAndMetrics.median_coverage
  }
}

# workflow AlignAndCallR2 {
#   meta {
#     description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
#   }

#   input {
#     Array[File] unmapped_bam
#     Array[String] sample_name
#     String suffix

#     File selfref_bundle
#     Array[String] mt_interval_list

#     Array[String] mt_self
#     Array[String] mt_self_index
#     Array[String] mt_self_dict

#     Array[String] self_cat
#     Array[String] self_cat_index
#     Array[String] self_cat_dict

#     Array[String] mt_self_shifted
#     Array[String] mt_self_shifted_index
#     Array[String] mt_self_shifted_dict

#     Array[String] self_shifted_cat
#     Array[String] self_shifted_cat_index
#     Array[String] self_shifted_cat_dict

#     Array[File] blacklisted_sites
#     Array[File] blacklisted_sites_index

#     Array[String] force_call_vcf
#     Array[String] force_call_vcf_idx
#     Array[String] force_call_vcf_shifted
#     Array[String] force_call_vcf_shifted_idx

#     Array[String] shift_back_chain

#     Array[String] non_control_interval
#     Array[String] control_shifted

#     Array[Float]? verifyBamID
#     Array[String] hasContamination
#     Array[Float] contamination_major
#     Array[Float] contamination_minor

#     File? gatk_override
#     String? gatk_docker_override
#     File JsonTools
#     String gatk_version = "4.2.6.0"
#     String? m2_extra_args
#     String? m2_filter_extra_args
#     Float? vaf_filter_threshold
#     Float? f_score_beta
#     Boolean compress_output_vcf

#     # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
#     # affected by this number. Default is 151.
#     Int? max_read_length

#     #Optional runtime arguments
#     Int? preemptible_tries
#     Int? n_cpu
#     Int? n_cpu_bwa
#   }

#   parameter_meta {
#     unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
#   }

#   call MongoTasks_Multi.MongoAlignToMtRegShiftedAndMetrics as AlignToMtRegShiftedAndMetrics {
#     input:
#       input_bam = unmapped_bam,
#       sample_base_name = sample_name,
#       suffix = suffix,

#       selfref_bundle = selfref_bundle,
#       mt = mt_self,
#       mt_index = mt_self_index,
#       mt_dict = mt_self_dict,

#       mt_cat = self_cat,
#       mt_cat_index = self_cat_index,
#       mt_cat_dict = self_cat_dict,

#       mt_shifted = mt_self_shifted,
#       mt_shifted_index = mt_self_shifted_index,
#       mt_shifted_dict = mt_self_shifted_dict,

#       mt_shifted_cat = self_shifted_cat,
#       mt_shifted_cat_index = self_shifted_cat_index,
#       mt_shifted_cat_dict = self_shifted_cat_dict,
      
#       mt_interval_list = mt_interval_list,

#       read_length = max_read_length,
#       coverage_cap = 100000,
      
#       preemptible_tries = preemptible_tries,
#       JsonTools = JsonTools,
#       n_cpu = n_cpu_bwa
#   }

#   Int M2_mem = if AlignToMtRegShiftedAndMetrics.max_mean_coverage > 25000 then 14 else 7

#   call MongoTasks_Multi.MongoCallMtAndShifted as CallMtAndShifted {
#     input:
#       sample_base_name = AlignToMtRegShiftedAndMetrics.samples,
#       suffix = suffix,
#       selfref_bundle = selfref_bundle,
#       # Everything is called except the control region.
#       input_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam,
#       input_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai,
#       mt_self = mt_self,
#       mt_self_index = mt_self_index,
#       mt_self_dict = mt_self_dict,
#       mt_interval_list = non_control_interval,
#       force_call_vcf = force_call_vcf,
#       force_call_vcf_idx = force_call_vcf_idx,

#       # Only the control region is now called.
#       shifted_input_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam,
#       #shifted_input_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai,
#       shifted_mt_self = mt_self_shifted,
#       shifted_mt_self_index = mt_self_shifted_index,
#       shifted_mt_self_dict = mt_self_shifted_dict,
#       shifted_mt_interval_list = control_shifted,
#       shifted_force_call_vcf = force_call_vcf_shifted,
#       shifted_force_call_vcf_idx = force_call_vcf_shifted_idx,

#       m2_extra_args = select_first([m2_extra_args, ""]),
#       shifted_m2_extra_args = select_first([m2_extra_args, ""]),

#       compress = compress_output_vcf,
#       gatk_override = gatk_override,
#       gatk_docker_override = gatk_docker_override,
#       gatk_version = gatk_version,
#       mem = M2_mem,
#       preemptible_tries = preemptible_tries,
#       JsonTools = JsonTools,
#       n_cpu = n_cpu
#   }

#   call MongoTasks_Multi.MongoLiftoverCombineMergeFilterContamSplit as LiftoverCombineMergeFilterContamSplit {
#     input:
#       sample_base_name = sample_name,
#       suffix = suffix,
#       selfref_bundle = selfref_bundle,

#       mt_self = mt_self,
#       mt_self_index = mt_self_index,
#       mt_self_dict = mt_self_dict,
#       shifted_vcf = CallMtAndShifted.shifted_raw_vcf,
#       shifted_vcf_idx = CallMtAndShifted.shifted_raw_vcf_idx,
#       non_shifted_vcf = CallMtAndShifted.raw_vcf,
#       non_shifted_vcf_idx = CallMtAndShifted.raw_vcf_idx,
#       shifted_stats = CallMtAndShifted.shifted_stats,
#       non_shifted_stats = CallMtAndShifted.stats,
#       shift_back_chain = shift_back_chain, 
#       blacklisted_sites = blacklisted_sites,
#       blacklisted_sites_index = blacklisted_sites_index,

#       hasContamination = hasContamination,
#       contamination_major = contamination_major,
#       contamination_minor = contamination_minor,
#       verifyBamID = verifyBamID,
      
#       compress = compress_output_vcf,
#       m2_extra_filtering_args = m2_filter_extra_args,
#       max_alt_allele_count = 4,
#       vaf_filter_threshold = vaf_filter_threshold,
#       f_score_beta = f_score_beta,
      
#       gatk_override = gatk_override,
#       gatk_docker_override = gatk_docker_override,
#       gatk_version = gatk_version,
#       JsonTools = JsonTools,
#       preemptible_tries = preemptible_tries
#   }

#   output {
#     Array[String] samples = LiftoverCombineMergeFilterContamSplit.samples
#     Array[File] mt_aligned_bam = AlignToMtRegShiftedAndMetrics.mt_aligned_bam
#     Array[File] mt_aligned_bai = AlignToMtRegShiftedAndMetrics.mt_aligned_bai
#     Array[File] mt_aligned_shifted_bam = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bam
#     #Array[File] mt_aligned_shifted_bai = AlignToMtRegShiftedAndMetrics.shifted_mt_aligned_bai
#     Array[File] nuc_mt_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bam
#     #Array[File] nuc_mt_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_mt_aligned_bai
#     Array[File] nuc_mt_shifted_aligned_bam = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bam
#     #Array[File] nuc_mt_shifted_aligned_bai = AlignToMtRegShiftedAndMetrics.nuc_and_shifted_mt_aligned_bai
#     Array[File] out_vcf = LiftoverCombineMergeFilterContamSplit.filtered_vcf
#     Array[File] out_vcf_idx = LiftoverCombineMergeFilterContamSplit.filtered_vcf_idx
#     Array[File] split_vcf = LiftoverCombineMergeFilterContamSplit.split_vcf
#     Array[File] split_vcf_idx = LiftoverCombineMergeFilterContamSplit.split_vcf_idx
#     Array[File] duplicate_metrics = AlignToMtRegShiftedAndMetrics.duplicate_metrics
#     Array[File] coverage_metrics = AlignToMtRegShiftedAndMetrics.wgs_metrics
#     Array[File] theoretical_sensitivity_metrics = AlignToMtRegShiftedAndMetrics.theoretical_sensitivity
#     Array[Int] mean_coverage = AlignToMtRegShiftedAndMetrics.mean_coverage
#     Array[Float] median_coverage = AlignToMtRegShiftedAndMetrics.median_coverage
#   }
# }