version 1.0

import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/WDL/v2.5_MongoSwirl_Single/fullMitoPipeline_v2_5_Single.wdl" as MitochondriaPipeline_v2_5

workflow MitochondriaPipelineWrapper {

  meta {
    description: "Takes in a list of hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    File wgs_aligned_input_bam_or_cram_list
    File wgs_aligned_input_bam_or_cram_index_list
    File sample_name_list

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
    File JsonTools
    File CheckHomOverlapScript
    File MergePerBatch

    Boolean force_manual_download
    String? requester_pays_project
    String? m2_extra_args
    String? m2_filter_extra_args
    String? printreads_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float? verifyBamID
    Boolean compress_output_vcf = false
    Boolean compute_numt_coverage = false
    Boolean use_haplotype_caller_nucdna = true
    Boolean skip_restore_hardclips = false
    Int haplotype_caller_nucdna_dp_lower_bound = 10

    #Docker and version arguments
    String gatk_version = "4.2.6.0"
    File? gatk_override
    String? gatk_docker_override
    String ucsc_docker = "docker.io/rahulg603/ucsc_genome_toolkit"
    String genomes_cloud_docker = "docker.io/rahulg603/genomes_cloud_bcftools"
    String haplochecker_docker = "docker.io/rahulg603/haplochecker"
    String gatk_samtools_docker = "docker.io/rahulg603/gatk46_samtools"

    #Optional runtime arguments
    Int? printreads_mem
    Int? lift_coverage_mem
    Int? n_cpu_subsetbam
    Int? n_cpu_m2_hc_lift
    Int? n_cpu_bwa
    Int? preemptible_tries
  }

  Array[File] input_cram = read_lines(wgs_aligned_input_bam_or_cram_list)
  Array[File] input_crai = read_lines(wgs_aligned_input_bam_or_cram_index_list)
  Array[String] sample_name = read_lines(sample_name_list)

  scatter (cram_pair in zip(zip(input_cram, input_crai), sample_name)) {
    call MitochondriaPipeline_v2_5.MitochondriaPipeline as MitochondriaPipeline_v2_5 {
        input:
          wgs_aligned_input_bam_or_cram = cram_pair.left.left,
          wgs_aligned_input_bam_or_cram_index = cram_pair.left.right,
          sample_name = cram_pair.right,
          mt_interval_list = mt_interval_list,
          nuc_interval_list = nuc_interval_list,

          max_read_length = max_read_length,
          skip_restore_hardclips = skip_restore_hardclips,

          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,

          mt_dict = mt_dict,
          mt_fasta = mt_fasta,
          mt_fasta_index = mt_fasta_index,
          blacklisted_sites = blacklisted_sites,
          blacklisted_sites_index = blacklisted_sites_index,

          control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
          non_control_region_interval_list = non_control_region_interval_list,

          HailLiftover = HailLiftover,
          FaRenamingScript = FaRenamingScript,
          CheckVariantBoundsScript = CheckVariantBoundsScript,
          CheckHomOverlapScript = CheckHomOverlapScript,
          JsonTools = JsonTools,

          force_manual_download = force_manual_download,
          requester_pays_project = requester_pays_project,
          m2_extra_args = m2_extra_args,
          m2_filter_extra_args = m2_filter_extra_args,
          printreads_extra_args = printreads_extra_args,
          vaf_filter_threshold = vaf_filter_threshold,
          f_score_beta = f_score_beta,
          verifyBamID = verifyBamID,
          compress_output_vcf = compress_output_vcf,
          compute_numt_coverage = compute_numt_coverage,
          use_haplotype_caller_nucdna = use_haplotype_caller_nucdna,
          haplotype_caller_nucdna_dp_lower_bound = haplotype_caller_nucdna_dp_lower_bound,

          gatk_version = gatk_version,
          gatk_override = gatk_override,
          gatk_docker_override = gatk_docker_override,
          ucsc_docker = ucsc_docker,
          genomes_cloud_docker = genomes_cloud_docker,
          haplochecker_docker = haplochecker_docker,
          gatk_samtools_docker = gatk_samtools_docker,

          printreads_mem = printreads_mem,
          lift_coverage_mem = lift_coverage_mem,
          n_cpu_subsetbam = n_cpu_subsetbam,
          n_cpu_m2_hc_lift = n_cpu_m2_hc_lift,
          n_cpu_bwa = n_cpu_bwa,
          preemptible_tries = preemptible_tries
    }
  }

  call MergeMitoMultiSampleOutputsInternal {
    input:
      sample_name = sample_name,
      variant_vcf = MitochondriaPipeline_v2_5.final_vcf,
      coverage_table = MitochondriaPipeline_v2_5.final_base_level_coverage_metrics,
      statistics = MitochondriaPipeline_v2_5.stats_outputs,
      MergePerBatch = MergePerBatch,
      preemptible_tries = preemptible_tries,
      genomes_cloud_docker = genomes_cloud_docker
  }

  output {
    Array[File] subset_bam = MitochondriaPipeline_v2_5.subset_bam
    Array[File] subset_bai = MitochondriaPipeline_v2_5.subset_bai
    Array[File] r1_vcf = MitochondriaPipeline_v2_5.r1_vcf
    Array[File] r1_vcf_index = MitochondriaPipeline_v2_5.r1_vcf_index
    Array[File] r1_nuc_vcf = MitochondriaPipeline_v2_5.r1_nuc_vcf
    Array[File] r1_nuc_vcf_index = MitochondriaPipeline_v2_5.r1_nuc_vcf_index
    Array[File] r1_nuc_vcf_unfiltered = MitochondriaPipeline_v2_5.r1_nuc_vcf_unfiltered
    Array[File] r1_split_vcf = MitochondriaPipeline_v2_5.r1_split_vcf
    Array[File] r1_split_vcf_index = MitochondriaPipeline_v2_5.r1_split_vcf_index

    Array[File] self_mt_aligned_bam = MitochondriaPipeline_v2_5.self_mt_aligned_bam
    Array[File] self_mt_aligned_bai = MitochondriaPipeline_v2_5.self_mt_aligned_bai
    Array[File] self_ref_vcf = MitochondriaPipeline_v2_5.self_ref_vcf
    Array[File] self_ref_vcf_index = MitochondriaPipeline_v2_5.self_ref_vcf_index
    Array[File] self_ref_split_vcf = MitochondriaPipeline_v2_5.self_ref_split_vcf
    Array[File] self_ref_split_vcf_index = MitochondriaPipeline_v2_5.self_ref_split_vcf_index
    Array[File] self_base_level_coverage_metrics = MitochondriaPipeline_v2_5.self_base_level_coverage_metrics
    Array[File] self_reference_fasta = MitochondriaPipeline_v2_5.self_reference_fasta
    Array[File] reference_to_self_ref_chain = MitochondriaPipeline_v2_5.reference_to_self_ref_chain
    Array[File] self_control_region_shifted = MitochondriaPipeline_v2_5.self_control_region_shifted
    Array[File] self_non_control_region = MitochondriaPipeline_v2_5.self_non_control_region

    Array[File] liftover_fix_pipeline_log = MitochondriaPipeline_v2_5.liftover_fix_pipeline_log
    Array[File] stats_outputs = MitochondriaPipeline_v2_5.stats_outputs

    Array[File] final_vcf = MitochondriaPipeline_v2_5.final_vcf
    Array[File] final_rejected_vcf = MitochondriaPipeline_v2_5.final_rejected_vcf
    Array[File] final_base_level_coverage_metrics = MitochondriaPipeline_v2_5.final_base_level_coverage_metrics
    Array[File?] numt_base_level_coverage = MitochondriaPipeline_v2_5.numt_base_level_coverage
    
    Array[File] input_vcf_for_haplochecker = MitochondriaPipeline_v2_5.input_vcf_for_haplochecker
    Array[File] duplicate_metrics = MitochondriaPipeline_v2_5.duplicate_metrics
    Array[File] coverage_metrics = MitochondriaPipeline_v2_5.coverage_metrics
    Array[File] theoretical_sensitivity_metrics = MitochondriaPipeline_v2_5.theoretical_sensitivity_metrics
    Array[File] contamination_metrics = MitochondriaPipeline_v2_5.contamination_metrics

    # liftover stats
    Array[Int] success_liftover_variants = MitochondriaPipeline_v2_5.success_liftover_variants
    Array[Int] failed_liftover_variants = MitochondriaPipeline_v2_5.failed_liftover_variants
    Array[Int] fixed_liftover_variants = MitochondriaPipeline_v2_5.fixed_liftover_variants
    Array[Int] n_liftover_r2_left_shift = MitochondriaPipeline_v2_5.n_liftover_r2_left_shift
    Array[Int] n_liftover_r2_injected_from_success = MitochondriaPipeline_v2_5.n_liftover_r2_injected_from_success
    Array[Int] n_liftover_r2_ref_insertion_new_haplo = MitochondriaPipeline_v2_5.n_liftover_r2_ref_insertion_new_haplo
    Array[Int] n_liftover_r2_failed_het_dele_span_insertion_boundary = MitochondriaPipeline_v2_5.n_liftover_r2_failed_het_dele_span_insertion_boundary
    Array[Int] n_liftover_r2_failed_new_dupes_leftshift = MitochondriaPipeline_v2_5.n_liftover_r2_failed_new_dupes_leftshift
    Array[Int] n_liftover_r2_het_ins_sharing_lhs_hom_dele = MitochondriaPipeline_v2_5.n_liftover_r2_het_ins_sharing_lhs_hom_dele
    Array[Int] n_liftover_r2_spanning_complex = MitochondriaPipeline_v2_5.n_liftover_r2_spanning_complex
    Array[Int] n_liftover_r2_spanningfixrhs_sharedlhs = MitochondriaPipeline_v2_5.n_liftover_r2_spanningfixrhs_sharedlhs
    Array[Int] n_liftover_r2_spanningfixlhs_upstream = MitochondriaPipeline_v2_5.n_liftover_r2_spanningfixlhs_upstream
    Array[Int] n_liftover_r2_repaired_success = MitochondriaPipeline_v2_5.n_liftover_r2_repaired_success

    # other stats
    Array[Int] mean_coverage = MitochondriaPipeline_v2_5.mean_coverage
    Array[Float] median_coverage = MitochondriaPipeline_v2_5.median_coverage
    Array[String] major_haplogroup = MitochondriaPipeline_v2_5.major_haplogroup
    Array[Float?] contamination = MitochondriaPipeline_v2_5.contamination
    Array[Int] nuc_variants_pass = MitochondriaPipeline_v2_5.nuc_variants_pass
    Array[Int] n_reads_unpaired_dropped = MitochondriaPipeline_v2_5.n_reads_unpaired_dropped
    Array[Int] nuc_variants_dropped = MitochondriaPipeline_v2_5.nuc_variants_dropped
    Array[Int] mtdna_consensus_overlaps = MitochondriaPipeline_v2_5.mtdna_consensus_overlaps
    Array[Int] nuc_consensus_overlaps = MitochondriaPipeline_v2_5.nuc_consensus_overlaps

    # merged files
    File merged_statistics = MergeMitoMultiSampleOutputsInternal.merged_statistics
    File merged_coverage = MergeMitoMultiSampleOutputsInternal.merged_coverage
    File merged_calls = MergeMitoMultiSampleOutputsInternal.merged_calls
    #File merged_idxstats = MergeMitoMultiSampleOutputs.merged_idxstats
    #File merged_yield = MergeMitoMultiSampleOutputs.merged_yield
  }
}

task MergeMitoMultiSampleOutputsInternal {
  # this task merges arrayed inputs for:
  # (1) output stats
  # (2) lifted-over variant calls
  # (3) per-base coverage
  input {
    Array[String] sample_name
    Array[File] variant_vcf
    Array[File] coverage_table
    Array[File] statistics

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
  }

  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: genomes_cloud_docker
    preemptible: select_first([preemptible_tries, 5])
  }
}
