version 1.0

import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/AlignAndCallR1_v2_6_Multi.wdl" as AlignAndCallR1_Multi
import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/AlignAndCallR2_v2_6_Multi.wdl" as AlignAndCallR2_Multi
import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/LiftoverTools_v2_6_Multi.wdl" as LiftoverTools_Multi
import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/ProduceSelfReferenceFiles_v2_6_Multi.wdl" as ProduceSelfReferenceFiles_Multi
import "https://raw.githubusercontent.com/gnchau/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/MongoTasks_v2_6_Multi.wdl" as MongoTasks_Multi

workflow MitochondriaPipeline {

  meta {
    description: "Takes in an hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    Array[File] wgs_aligned_input_bam_or_cram
    Array[File] wgs_aligned_input_bam_or_cram_index
    Array[String] sample_name

    File mt_interval_list
    File nuc_interval_list

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File blacklisted_sites
    File blacklisted_sites_index

    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list

    File HailLiftover
    File FaRenamingScript
    File CheckVariantBoundsScript
    File CheckHomOverlapScript
    File JsonTools
    File MergePerBatch

    Boolean force_manual_download
    String? requester_pays_project
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Array[Float]+? verifyBamID
    Boolean compress_output_vcf = false
    Boolean compute_numt_coverage = false
    Boolean use_haplotype_caller_nucdna = true
    Int haplotype_caller_nucdna_dp_lower_bound = 10    
    
    # Some CRAMs (e.g., AoU) contain an XQ tag that doesn't play well with gatk RevertSam.
    # Enable this flag to avoid using this tag by skipping the restore hardclips step.
    Boolean skip_restore_hardclips = false

    #Docker and version arguments
    String gatk_version = "4.2.6.0"
    File? gatk_override
    String? gatk_docker_override
    String ucsc_docker
    String genomes_cloud_docker
    String haplochecker_docker
    String gatk_samtools_docker

    #Optional runtime arguments
    Int? printreads_mem
    Int? lift_coverage_mem
    Int? n_cpu_subsetbam
    Int? n_cpu_m2_hc_lift
    Int? n_cpu_bwa
    Int? preemptible_tries
  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    mt_interval_list: "Picard style interval list file, with header and single interval representing chrM, eg chrM 1 16569 + ., and putative NUMT intervals."
  }

  String self_ref_suffix = ".self.ref"

  call MongoTasks_Multi.MongoSubsetBam as SubsetBam {
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      input_bai = wgs_aligned_input_bam_or_cram_index,
      sample_name = sample_name,

      mt_interval_list = mt_interval_list,
      nuc_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      requester_pays_project = requester_pays_project,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_samtools_docker,
      gatk_version = gatk_version,
      force_manual_download = force_manual_download,
      mem = printreads_mem,
      JsonTools = JsonTools,
      n_cpu = n_cpu_subsetbam,
      preemptible_tries = preemptible_tries
  }

  call MongoTasks_Multi.MongoProcessBamAndRevert as ProcessBam {
    input:
      subset_bam = SubsetBam.subset_bam,
      subset_bai = SubsetBam.subset_bai,
      flagstat_pre_metrics = SubsetBam.flagstat_pre_metrics,
      sample_name = SubsetBam.samples,

      mt_interval_list = mt_interval_list,
      nuc_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      read_length = max_read_length,
      skip_restore_hardclips = skip_restore_hardclips,
      coverage_cap = 100000,
      JsonTools = JsonTools,
      n_cpu = n_cpu_bwa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallR1_Multi.AlignAndCallR1 as AlignAndCallR1 {
    input:
      input_bam = ProcessBam.output_bam,
      input_bai = ProcessBam.output_bai,
      sample_name = ProcessBam.samples,

      mt_interval_list = mt_interval_list,
      nuc_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_mean_coverage = ProcessBam.max_mean_coverage,
      mt_mean_coverage_array = ProcessBam.mean_coverage,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_args = m2_extra_args,
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      verifyBamID = verifyBamID,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      use_haplotype_caller_nucdna = use_haplotype_caller_nucdna,
      hc_dp_lower_bound = haplotype_caller_nucdna_dp_lower_bound,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries,
      haplochecker_docker = haplochecker_docker,
      n_cpu = n_cpu_m2_hc_lift
  }

  call ProduceSelfReferenceFiles_Multi.ProduceSelfReferenceFiles as ProduceSelfRefFiles {
    input:
      sample_name = AlignAndCallR1.samples,
      mtdna_variants = AlignAndCallR1.split_vcf,
      nuc_variants = AlignAndCallR1.split_nuc_vcf,
      suffix = self_ref_suffix,
      
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_dict = mt_dict,
      mt_interval_list = mt_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      nuc_interval_list = nuc_interval_list,
      reference_name = "reference",
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      n_shift = 8000,
      compute_numt_coverage = compute_numt_coverage,
      FaRenamingScript = FaRenamingScript,
      CheckVariantBoundsScript = CheckVariantBoundsScript,
      CheckHomOverlapScript = CheckHomOverlapScript,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries,
      genomes_cloud_docker = genomes_cloud_docker,
      ucsc_docker = ucsc_docker
  }

  call AlignAndCallR2_Multi.AlignAndCallR2 as AlignAndCallR2 {
    input:
      unmapped_bam = ProcessBam.unmapped_bam,
      sample_name = ProduceSelfRefFiles.samples,
      suffix = self_ref_suffix,

      selfref_bundle = ProduceSelfRefFiles.selfref_bundle,
      mt_interval_list = ProduceSelfRefFiles.mt_interval_list_self,

      mt_self = ProduceSelfRefFiles.mt_self,
      mt_self_index = ProduceSelfRefFiles.mt_self_index,
      mt_self_dict = ProduceSelfRefFiles.mt_self_dict,
      self_cat = ProduceSelfRefFiles.mt_andNuc_self,
      self_cat_index = ProduceSelfRefFiles.mt_andNuc_self_index,
      self_cat_dict = ProduceSelfRefFiles.mt_andNuc_self_dict,
      mt_self_shifted = ProduceSelfRefFiles.mt_shifted_self,
      mt_self_shifted_index = ProduceSelfRefFiles.mt_shifted_self_index,
      mt_self_shifted_dict = ProduceSelfRefFiles.mt_shifted_self_dict,
      self_shifted_cat = ProduceSelfRefFiles.mt_andNuc_shifted_self,
      self_shifted_cat_index = ProduceSelfRefFiles.mt_andNuc_shifted_self_index,
      self_shifted_cat_dict = ProduceSelfRefFiles.mt_andNuc_shifted_self_dict,
      shift_back_chain = ProduceSelfRefFiles.self_shift_back_chain,

      force_call_vcf = ProduceSelfRefFiles.force_call_vcf,
      force_call_vcf_idx = ProduceSelfRefFiles.force_call_vcf_idx,
      force_call_vcf_shifted = ProduceSelfRefFiles.force_call_vcf_shifted,
      force_call_vcf_shifted_idx = ProduceSelfRefFiles.force_call_vcf_shifted_idx,

      non_control_interval = ProduceSelfRefFiles.non_control_interval_self,
      control_shifted = ProduceSelfRefFiles.control_shifted_self,
      blacklisted_sites = ProduceSelfRefFiles.blacklisted_sites_self,
      blacklisted_sites_index = ProduceSelfRefFiles.blacklisted_sites_index_self,
      hasContamination = AlignAndCallR1.hasContamination,
      contamination_major = AlignAndCallR1.contamination_major,
      contamination_minor = AlignAndCallR1.contamination_minor,
      
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_args = select_first([m2_extra_args," "]),
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      verifyBamID = verifyBamID,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries,
      n_cpu_bwa = n_cpu_bwa,
      n_cpu = n_cpu_m2_hc_lift
  }

  call MongoTasks_Multi.MongoLiftoverVCFAndGetCoverage as LiftOverAfterSelf {
    input:
      sample_name = AlignAndCallR2.samples,
      selfref_bundle = ProduceSelfRefFiles.selfref_bundle,
      ref_homoplasmies_vcf = ProduceSelfRefFiles.ref_homoplasmies_vcf,
      r2_self_ref_vcf = AlignAndCallR2.split_vcf,
      self_homoplasmies_vcf = ProduceSelfRefFiles.force_call_vcf_filters,

      mt_self = ProduceSelfRefFiles.mt_self,
      mt_self_index = ProduceSelfRefFiles.mt_self_index,
      mt_self_dict = ProduceSelfRefFiles.mt_self_dict,
      mt_self_shifted = ProduceSelfRefFiles.mt_shifted_self,
      mt_self_shifted_index = ProduceSelfRefFiles.mt_shifted_self_index,
      mt_self_shifted_dict = ProduceSelfRefFiles.mt_shifted_self_dict,
      chain_self_to_ref = ProduceSelfRefFiles.self_to_ref_chain,
      chain_ref_to_self = ProduceSelfRefFiles.ref_to_self_chain,

      input_bam_regular_ref = AlignAndCallR2.mt_aligned_bam,
      input_bam_regular_ref_index = AlignAndCallR2.mt_aligned_bai,
      input_bam_shifted_ref = AlignAndCallR2.mt_aligned_shifted_bam,
      #input_bam_shifted_ref_index = AlignAndCallR2.mt_aligned_shifted_bai,
      self_control_region_shifted_reference_interval_list = ProduceSelfRefFiles.control_shifted_self,
      self_non_control_region_interval_list = ProduceSelfRefFiles.non_control_interval_self,
      
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      HailLiftover = HailLiftover,
      JsonTools = JsonTools,
      self_suffix = self_ref_suffix,

      n_cpu = n_cpu_m2_hc_lift,
      genomes_cloud_docker = genomes_cloud_docker,
      preemptible_tries = preemptible_tries
  }

  call MongoTasks_Multi.MongoLiftoverSelfAndCollectOutputs as LiftoverSelfCoverage {
    input:
      sample_name = LiftOverAfterSelf.samples,
      self_ref_table = LiftOverAfterSelf.self_coverage_table,
      chain = ProduceSelfRefFiles.self_to_ref_chain,
      homoplasmic_deletions_coverage = LiftOverAfterSelf.gap_coverage,
      
      liftover_table = LiftOverAfterSelf.liftoverStats,
      mean_coverage = AlignAndCallR2.mean_coverage,
      median_coverage = AlignAndCallR2.median_coverage,
      major_haplogroup = AlignAndCallR1.major_haplogroup,
      contamination = AlignAndCallR1.contamination,
      nuc_variants_pass = AlignAndCallR1.nuc_variants_pass,
      n_reads_unpaired_dropped = ProcessBam.reads_dropped,
      nuc_variants_dropped = ProduceSelfRefFiles.nuc_variants_dropped,
      mtdna_consensus_overlaps = ProduceSelfRefFiles.mtdna_consensus_overlaps,
      nuc_consensus_overlaps = ProduceSelfRefFiles.nuc_consensus_overlaps,
      JsonTools = JsonTools,
      ucsc_docker = ucsc_docker,
      preemptible_tries = preemptible_tries
  }

  if (compute_numt_coverage) {
    call NucCoverageAtEveryBase {
      input:
        sample_name = AlignAndCallR2.samples,
        selfref_bundle = ProduceSelfRefFiles.selfref_bundle,
        input_bam_regular_ref = ProcessBam.output_bam,
        input_bam_regular_ref_index = ProcessBam.output_bai,
        input_bam_self_ref = AlignAndCallR2.nuc_mt_aligned_bam,
        #input_bam_self_ref_index = AlignAndCallR2.nuc_mt_aligned_bai,
        input_bam_self_ref_shifted = AlignAndCallR2.nuc_mt_shifted_aligned_bam,
        #input_bam_self_ref_shifted_index = AlignAndCallR2.nuc_mt_shifted_aligned_bai,
        nuc_interval_list = nuc_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        self_nuc_interval_list = ProduceSelfRefFiles.nuc_interval_list_self,
        self_fasta = ProduceSelfRefFiles.mt_andNuc_self,
        self_fasta_index = ProduceSelfRefFiles.mt_andNuc_self_index,
        self_dict = ProduceSelfRefFiles.mt_andNuc_self_dict,
        self_nuc_interval_list_shifted = ProduceSelfRefFiles.nuc_interval_list_shifted_self,
        self_shifted_fasta = ProduceSelfRefFiles.mt_andNuc_shifted_self,
        self_shifted_fasta_index = ProduceSelfRefFiles.mt_andNuc_shifted_self_index,
        self_shifted_dict = ProduceSelfRefFiles.mt_andNuc_shifted_self_dict,
        JsonTools = JsonTools,
        preemptible_tries = preemptible_tries
    }

    call LiftoverTools_Multi.LiftOverAndJoinCoverage as LiftOverAndJoinCoverage {
      input:
        sample_name = NucCoverageAtEveryBase.samples,
        ref_table = NucCoverageAtEveryBase.table_old_ref,
        self_table = NucCoverageAtEveryBase.table_new_self,
        self_table_shifted = NucCoverageAtEveryBase.table_new_self_shifted,
        chain = ProduceSelfRefFiles.nuc_self_to_ref_chain,
        mem = lift_coverage_mem,
        JsonTools = JsonTools,
        ucsc_docker = ucsc_docker,
        preemptible_tries = preemptible_tries
    }
  }

  call ExtractProduceSelfRefOutputs {
    input:
      sample_name = LiftoverSelfCoverage.samples,
      selfref_bundle = ProduceSelfRefFiles.selfref_bundle,
      mt_self = ProduceSelfRefFiles.mt_self,
      ref_to_self_chain = ProduceSelfRefFiles.ref_to_self_chain,
      control_shifted_self = ProduceSelfRefFiles.control_shifted_self,
      non_control_interval_self = ProduceSelfRefFiles.non_control_interval_self,
      JsonTools = JsonTools,
      haplochecker_docker = haplochecker_docker,
      preemptible_tries = preemptible_tries
  }

  call MergeMitoMultiSampleOutputs {
    input:
      sample_name = ExtractProduceSelfRefOutputs.samples,
      variant_vcf = LiftOverAfterSelf.liftover_r2_final_vcf,
      coverage_table = LiftoverSelfCoverage.reference_coverage,
      statistics = LiftoverSelfCoverage.table,
      idxstats_metrics = SubsetBam.idxstats_metrics, # changed from SubsetBamToChrMAndRevert
      yield_metrics = ProcessBam.yield_metrics,
      MergePerBatch = MergePerBatch,
      genomes_cloud_docker = genomes_cloud_docker,
      preemptible_tries = preemptible_tries
  }

  output { # need to edit
    Array[File] subset_bam = ProcessBam.output_bam
    Array[File] subset_bai = ProcessBam.output_bai
    Array[File] r1_vcf = AlignAndCallR1.out_vcf
    Array[File] r1_vcf_index = AlignAndCallR1.out_vcf_index
    Array[File] r1_nuc_vcf = AlignAndCallR1.nuc_vcf
    Array[File] r1_nuc_vcf_index = AlignAndCallR1.nuc_vcf_index
    Array[File] r1_nuc_vcf_unfiltered = AlignAndCallR1.nuc_vcf_unfiltered
    Array[File] r1_split_vcf = AlignAndCallR1.split_vcf
    Array[File] r1_split_vcf_index = AlignAndCallR1.split_vcf_index

    Array[File] self_mt_aligned_bam = AlignAndCallR2.mt_aligned_bam
    Array[File] self_mt_aligned_bai = AlignAndCallR2.mt_aligned_bai
    Array[File] self_ref_vcf = AlignAndCallR2.out_vcf
    Array[File] self_ref_vcf_index = AlignAndCallR2.out_vcf_idx
    Array[File] self_ref_split_vcf = AlignAndCallR2.split_vcf
    Array[File] self_ref_split_vcf_index = AlignAndCallR2.split_vcf_idx
    Array[File] self_base_level_coverage_metrics = LiftOverAfterSelf.self_coverage_table
    Array[File] self_reference_fasta = ExtractProduceSelfRefOutputs.mt_self_out
    Array[File] reference_to_self_ref_chain = ExtractProduceSelfRefOutputs.ref_to_self_chain_out
    Array[File] self_control_region_shifted = ExtractProduceSelfRefOutputs.control_shifted_self_out
    Array[File] self_non_control_region = ExtractProduceSelfRefOutputs.non_control_interval_self_out

    Array[File] liftover_fix_pipeline_log = LiftOverAfterSelf.liftover_r2_log
    Array[File] stats_outputs = LiftoverSelfCoverage.table

    Array[File] final_vcf = LiftOverAfterSelf.liftover_r2_final_vcf
    Array[File] final_rejected_vcf = LiftOverAfterSelf.liftover_r2_rejected_vcf
    Array[File] final_base_level_coverage_metrics = LiftoverSelfCoverage.reference_coverage
    Array[File]? numt_base_level_coverage = LiftOverAndJoinCoverage.reference_coverage
    
    Array[File] input_vcf_for_haplochecker = AlignAndCallR1.input_vcf_for_haplochecker
    Array[File] duplicate_metrics = AlignAndCallR2.duplicate_metrics
    Array[File] coverage_metrics = AlignAndCallR2.coverage_metrics
    Array[File] theoretical_sensitivity_metrics = AlignAndCallR2.theoretical_sensitivity_metrics
    Array[File] contamination_metrics = AlignAndCallR1.contamination_metrics

    # liftover stats – the full set are outputted in the stats outputs file above
    Array[Int] success_liftover_variants = LiftOverAfterSelf.n_liftover_r2_pass
    Array[Int] failed_liftover_variants = LiftOverAfterSelf.n_liftover_r2_failed
    Array[Int] fixed_liftover_variants = LiftOverAfterSelf.n_liftover_r2_fixed
    Array[Int] n_liftover_r2_left_shift = LiftOverAfterSelf.n_liftover_r2_left_shift
    Array[Int] n_liftover_r2_injected_from_success = LiftOverAfterSelf.n_liftover_r2_injected_from_success
    Array[Int] n_liftover_r2_ref_insertion_new_haplo = LiftOverAfterSelf.n_liftover_r2_ref_insertion_new_haplo
    Array[Int] n_liftover_r2_failed_het_dele_span_insertion_boundary = LiftOverAfterSelf.n_liftover_r2_failed_het_dele_span_insertion_boundary
    Array[Int] n_liftover_r2_failed_new_dupes_leftshift = LiftOverAfterSelf.n_liftover_r2_failed_new_dupes_leftshift
    Array[Int] n_liftover_r2_het_ins_sharing_lhs_hom_dele = LiftOverAfterSelf.n_liftover_r2_het_ins_sharing_lhs_hom_dele
    Array[Int] n_liftover_r2_spanning_complex = LiftOverAfterSelf.n_liftover_r2_spanning_complex
    Array[Int] n_liftover_r2_spanningfixrhs_sharedlhs = LiftOverAfterSelf.n_liftover_r2_spanningfixrhs_sharedlhs
    Array[Int] n_liftover_r2_spanningfixlhs_upstream = LiftOverAfterSelf.n_liftover_r2_spanningfixlhs_upstream
    Array[Int] n_liftover_r2_repaired_success = LiftOverAfterSelf.n_liftover_r2_repaired_success

    # other stats and constants; also in stats_outputs
    Array[Int] mean_coverage = AlignAndCallR2.mean_coverage
    Array[Float] median_coverage = AlignAndCallR2.median_coverage
    Array[String] major_haplogroup = AlignAndCallR1.major_haplogroup
    Array[Float] contamination = AlignAndCallR1.contamination
    Array[Int] nuc_variants_pass = AlignAndCallR1.nuc_variants_pass
    Array[Int] n_reads_unpaired_dropped = ProcessBam.reads_dropped
    Array[Int] nuc_variants_dropped = ProduceSelfRefFiles.nuc_variants_dropped
    Array[Int] mtdna_consensus_overlaps = ProduceSelfRefFiles.mtdna_consensus_overlaps
    Array[Int] nuc_consensus_overlaps = ProduceSelfRefFiles.nuc_consensus_overlaps

    # merged files
    File merged_statistics = MergeMitoMultiSampleOutputs.merged_statistics
    File merged_coverage = MergeMitoMultiSampleOutputs.merged_coverage
    File merged_calls = MergeMitoMultiSampleOutputs.merged_calls
    File merged_idxstats = MergeMitoMultiSampleOutputs.merged_idxstats
    File merged_yield = MergeMitoMultiSampleOutputs.merged_yield
  }
}

task NucCoverageAtEveryBase {
  input {
    Array[String] sample_name
    File selfref_bundle
    Array[File] input_bam_regular_ref
    Array[File] input_bam_regular_ref_index
    Array[File] input_bam_self_ref
    #Array[File] input_bam_self_ref_index
    Array[File] input_bam_self_ref_shifted
    #Array[File] input_bam_self_ref_shifted_index
    Array[String] self_fasta
    Array[String] self_fasta_index
    Array[String] self_dict
    Array[File]? self_nuc_interval_list
    Array[String] self_shifted_fasta
    Array[String] self_shifted_fasta_index
    Array[String] self_shifted_dict
    Array[File]? self_nuc_interval_list_shifted

    File nuc_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File JsonTools

    Int? preemptible_tries

    # never specify the below input!
    Array[File]? tmp
  }

  # This function is only called when compute_numt_coverage is True thus
  # self_nuc_interval_list/self_nuc_interval_list_shifted should never be missing.
  Array[File] self_nuc_interval_list_internal = select_first([self_nuc_interval_list, tmp])
  Array[File] self_nuc_interval_list_shifted_internal = select_first([self_nuc_interval_list_shifted, tmp])

  Int disk_size = ceil(size(input_bam_regular_ref, "GB") + size(input_bam_self_ref, "GB") + size(ref_fasta, "GB") * 2) + 20
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  meta {
    description: "Mainly for QC to understand how NUMT coverage changes with remapping."
  }
  command <<<
    set -e

    tar xf "~{selfref_bundle}"

    sampleNames=('~{sep="' '" sample_name}')
    ref_bams=('~{sep="' '" input_bam_regular_ref}')
    self_fasta=('~{sep="' '" self_fasta}')
    self_fasta_shifted=('~{sep="' '" self_shifted_fasta}')
    self_bams=('~{sep="' '" input_bam_self_ref}')
    self_shifted_bams=('~{sep="' '" input_bam_self_ref_shifted}')
    self_intervals=('~{sep="' '" self_nuc_interval_list_internal}')
    self_intervals_shifted=('~{sep="' '" self_nuc_interval_list_shifted_internal}')

    mkdir out

    for i in "~{d}{!sampleNames[@]}"; do

      this_sample=out/"~{d}{sampleNames[i]}"
      this_bam="~{d}{ref_bams[i]}"
      this_self_bam="~{d}{self_bams[i]}"
      this_self_shifted_bam="~{d}{self_shifted_bams[i]}"
      this_self_fasta="~{d}{self_fasta[i]}"
      this_self_fasta_shifted="~{d}{self_fasta_shifted[i]}"
      this_self_intervals="~{d}{self_intervals[i]}"
      this_self_intervals_shifted="~{d}{self_intervals_shifted[i]}"

      samtools index "~{d}{this_self_bam}"
      samtools index "~{d}{this_self_shifted_bam}"

      java -jar /usr/gitc/picard.jar CollectHsMetrics \
        I="~{d}{this_bam}" \
        R=~{ref_fasta} \
        PER_BASE_COVERAGE="~{d}{this_sample}.nuc.reference_pre_realignment.tsv" \
        O=reference_pre_realignment.metrics \
        TI=~{nuc_interval_list} \
        BI=~{nuc_interval_list} \
        COVMAX=20000 \
        SAMPLE_SIZE=1

      java -jar /usr/gitc/picard.jar CollectHsMetrics \
        I="~{d}{this_self_bam}" \
        R="~{d}{this_self_fasta}" \
        PER_BASE_COVERAGE="~{d}{this_sample}.nuc.self_post_realignment.tsv" \
        O=self_post_realignment.metrics \
        TI="~{d}{this_self_intervals}" \
        BI="~{d}{this_self_intervals}" \
        COVMAX=20000 \
        SAMPLE_SIZE=1

      java -jar /usr/gitc/picard.jar CollectHsMetrics \
        I="~{d}{this_self_shifted_bam}" \
        R="~{d}{this_self_fasta_shifted}" \
        PER_BASE_COVERAGE="~{d}{this_sample}.nuc.self_post_realignment_shifted.tsv" \
        O=self_post_realignment_shifted.metrics \
        TI="~{d}{this_self_intervals_shifted}" \
        BI="~{d}{this_self_intervals_shifted}" \
        COVMAX=20000 \
        SAMPLE_SIZE=1

      python ~{JsonTools} \
      --path out/jsonout.json \
      --set samples="~{d}{sampleNames[i]}" \
        table_old_ref="~{d}{this_sample}.nuc.reference_pre_realignment.tsv" \
        table_new_self="~{d}{this_sample}.nuc.self_post_realignment.tsv" \
        table_new_self_shifted="~{d}{this_sample}.nuc.self_post_realignment_shifted.tsv"

    done
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2000 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    Object obj_out = read_json('out/jsonout.json')
    Array[String] samples = obj_out.samples
    Array[File] table_old_ref = obj_out.table_old_ref
    Array[File] table_new_self = obj_out.table_new_self
    Array[File] table_new_self_shifted = obj_out.table_new_self_shifted
  }
}

task TransposeTable {
  input {
    File input_table
    String genomes_cloud_docker
  }

  String output_table = basename(input_table, '.txt') + '_t.txt'

  command <<<
    set -e

    R --vanilla <<CODE
      df = read.csv("~{input_table}", header=F, sep='\t')
      swapped_df = as.data.frame(t(df))
      write.table(swapped_df, sep ='\t', row.names = F, col.names = F, file = "~{output_table}", quote = F)
    CODE
  >>>

  output {
    File table = "~{output_table}"
  }

  runtime {
    disks: "local-disk 10 HDD"
    memory: "1 GB"
    docker: genomes_cloud_docker
    preemptible: 5
  }
}

task MergeMitoMultiSampleOutputs {
  # this task merges arrayed inputs for:
  # (1) output stats
  # (2) lifted-over variant calls
  # (3) per-base coverage
  input {
    Array[String] sample_name
    Array[File] variant_vcf
    Array[File] coverage_table
    Array[File] statistics
    Array[File] idxstats_metrics
    Array[File] yield_metrics

    File MergePerBatch
    Int? preemptible_tries
    String genomes_cloud_docker
  }

  Int disk_size = ceil(size(variant_vcf, "GB") + size(coverage_table, "GB") + size(statistics, "GB") * 2) + 20

  command <<<
    set -e

    # start by merging statistics
    R --vanilla <<CODE
      files_of_interest <- read.csv("~{write_lines(statistics)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      dfs <- lapply(files_of_interest, function(x)read.csv(x, sep='\t', stringsAsFactors=F))
      df <- do.call("rbind", dfs)
      write.table(df, sep ='\t', row.names = F, file = "batch_analysis_statistics.tsv", quote = F)
    CODE

    # now produce the relevant inputs for the per-batch MT script
    R --vanilla <<CODE
      sample_ids <- read.csv("~{write_lines(sample_name)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      paths_coverage <- read.csv("~{write_lines(coverage_table)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      paths_vcf <- read.csv("~{write_lines(variant_vcf)}", sep='\t', stringsAsFactors=F, header=F)[[1]]

      write.table(data.frame(s=sample_ids, path=paths_coverage), sep='\t', row.names=F, file='coverage_paths.tsv', quote=F)
      write.table(data.frame(s=sample_ids, path=paths_vcf), sep='\t', row.names=F, file='vcf_paths.tsv', quote=F)
    CODE

    # merge idxstats stuff
    R --vanilla <<CODE
      files_of_interest <- read.csv("~{write_lines(idxstats_metrics)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      samples <- read.csv("~{write_lines(sample_name)}", stringsAsFactors=F, header=F)[[1]]
      dfs <- lapply(1:length(files_of_interest), function(idx) {
        tab <- read.csv(files_of_interest[idx], sep='\t', stringsAsFactors=F, col.names=c('chr','len','mapped_reads','unmapped_reads'), header=F)
        tab[,'s'] <- samples[idx]
        return(tab)
      })
      df <- do.call("rbind", dfs)
      gz1 <- gzfile("batch_idxstats_metrics.tsv.gz", 'w')
      write.table(df, sep ='\t', row.names = F, file = gz1, quote = F)
      close(gz1)
    CODE

    # merge yield stuff
    R --vanilla <<CODE
      files_of_interest <- read.csv("~{write_lines(yield_metrics)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      samples <- read.csv("~{write_lines(sample_name)}", stringsAsFactors=F, header=F)[[1]]
      dfs <- lapply(1:length(files_of_interest), function(idx) {
        tab <- read.csv(files_of_interest[idx], sep='\t', stringsAsFactors=F)
        tab[,'s'] <- samples[idx]
        return(tab)
      })
      df <- do.call("rbind", dfs)
      gz1 <- gzfile("batch_yield_metrics.tsv.gz", 'w')
      write.table(df, sep ='\t', row.names = F, file = gz1, quote = F)
      close(gz1)
    CODE

    mkdir tmp

    # now merge coverage
    python3.7 ~{MergePerBatch} \
    --run-coverage \
    --input-tsv coverage_paths.tsv \
    --temp-dir tmp/ \
    --output-flat-file batch_merged_mt_coverage.tsv.bgz

    # finally merge variant calls
    python3.7 ~{MergePerBatch} \
    --run-variants \
    --input-tsv vcf_paths.tsv \
    --temp-dir tmp/ \
    --output-flat-file batch_merged_mt_calls.vcf.bgz
  >>>

  output {
    File merged_statistics = "batch_analysis_statistics.tsv"
    File merged_coverage = "batch_merged_mt_coverage.tsv.bgz"
    File merged_calls = "batch_merged_mt_calls.vcf.bgz"
    File merged_idxstats = "batch_idxstats_metrics.tsv.gz"
    File merged_yield = "batch_yield_metrics.tsv.gz"
  }

  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: genomes_cloud_docker
    preemptible: select_first([preemptible_tries, 5])
  }
}

task ExtractProduceSelfRefOutputs {
  input {
    Array[String] sample_name
    File selfref_bundle
    Array[String] mt_self
    Array[String] ref_to_self_chain
    Array[String] control_shifted_self
    Array[String] non_control_interval_self
    
    File JsonTools
    String haplochecker_docker
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(selfref_bundle, "GB") * 5)
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  command <<<
    set -e

    mkdir out

    tar xf "~{selfref_bundle}"

    sampleNames=('~{sep="' '" sample_name}')
    self_fasta=('~{sep="' '" mt_self}')
    chains=('~{sep="' '" ref_to_self_chain}')
    control_shifteds=('~{sep="' '" control_shifted_self}')
    non_controls=('~{sep="' '" non_control_interval_self}')

    for i in "~{d}{!sampleNames[@]}"; do

      this_self_fasta="~{d}{self_fasta[i]}"
      this_self_chain="~{d}{chains[i]}"
      this_self_control_shifted="~{d}{control_shifteds[i]}"
      this_self_noncontrol="~{d}{non_controls[i]}"

      python ~{JsonTools} \
      --path out/jsonout.json \
      --set samples="~{d}{sampleNames[i]}" \
        mt_self_out="~{d}{this_self_fasta}" \
        ref_to_self_chain_out="~{d}{this_self_chain}" \
        control_shifted_self_out="~{d}{this_self_control_shifted}" \
        non_control_interval_self_out="~{d}{this_self_noncontrol}"

    done
  >>>

  output {
    Object obj_out = read_json('out/jsonout.json')
    Array[String] samples = obj_out.samples
    Array[File] mt_self_out = obj_out.mt_self_out
    Array[File] ref_to_self_chain_out = obj_out.ref_to_self_chain_out
    Array[File] control_shifted_self_out = obj_out.control_shifted_self_out
    Array[File] non_control_interval_self_out = obj_out.non_control_interval_self_out
  }

  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: haplochecker_docker
  }
}
