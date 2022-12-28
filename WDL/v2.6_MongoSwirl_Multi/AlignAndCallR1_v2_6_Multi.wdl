version 1.0

import "https://raw.githubusercontent.com/rahulg603/testing-mito-wdl/master/WDL/v2.6_MongoSwirl_Multi/MongoTasks_v2_6_Multi.wdl" as MongoTasks_Multi

workflow AlignAndCallR1 {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    Array[File] input_bam
    Array[File] input_bai
    Array[String] sample_name

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File blacklisted_sites
    File blacklisted_sites_index

    File nuc_interval_list
    File mt_interval_list

    Int mt_mean_coverage
    Array[Int] mt_mean_coverage_array

    Boolean use_haplotype_caller_nucdna
    Int hc_dp_lower_bound
    File? gatk_override
    String? gatk_docker_override
    String gatk_version = "4.2.6.0"
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf
    Array[Float]+? verifyBamID    
    File JsonTools

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length
    String haplochecker_docker

    #Optional runtime arguments
    Int? preemptible_tries
    Int? n_cpu
  }

  parameter_meta {
  }

  if (use_haplotype_caller_nucdna) {
    call MongoTasks_Multi.MongoHC as CallNucHCIntegrated {
      input:
        input_bam = input_bam,
        input_bai = input_bai,
        sample_name = sample_name,
        nuc_interval_list = nuc_interval_list,
        ref_fasta = ref_fasta,
        ref_fai = ref_fasta_index,
        ref_dict = ref_dict,
        suffix = '.nuc',
        compress = compress_output_vcf,
        gatk_override = gatk_override,
        gatk_docker_override = gatk_docker_override,
        gatk_version = gatk_version, 
        hc_dp_lower_bound = hc_dp_lower_bound,
        mem = 4,
        preemptible_tries = preemptible_tries,
        JsonTools = JsonTools,
        n_cpu = n_cpu
    }
  }
  if (!use_haplotype_caller_nucdna) {
    call MongoTasks_Multi.MongoNucM2 as CallNucM2Integrated {
      input:
        input_bam = input_bam,
        input_bai = input_bai,
        sample_name = sample_name,

        ref_fasta = ref_fasta,
        ref_fai = ref_fasta_index,
        ref_dict = ref_dict,
        suffix = '.nuc',
        mt_interval_list = nuc_interval_list,

        compress = compress_output_vcf,
        m2_extra_args = select_first([m2_extra_args, ""]),

        max_alt_allele_count = 4,
        vaf_filter_threshold = 0.95,
        verifyBamID = verifyBamID,
        compress = compress_output_vcf,

        gatk_override = gatk_override,
        gatk_docker_override = gatk_docker_override,
        gatk_version = gatk_version,
        mem = 4,
        preemptible_tries = preemptible_tries,
        JsonTools = JsonTools,
        n_cpu = n_cpu
    }
  }
  
  Int M2_mem = if mt_mean_coverage > 25000 then 14 else 7

  call MongoTasks_Multi.MongoRunM2InitialFilterSplit as CallMt {
    input:
      #sample_name = sample_name,
      sample_name = select_first([CallNucHCIntegrated.samples, CallNucM2Integrated.samples]),
      input_bam = input_bam,
      input_bai = input_bai,
      verifyBamID = verifyBamID,
      mt_interval_list = mt_interval_list,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      suffix = "",
      compress = compress_output_vcf,
      m2_extra_filtering_args = select_first([m2_filter_extra_args, ""]) + " --min-median-mapping-quality 0",
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_args = select_first([m2_extra_args, ""]),
      mem = M2_mem,
      preemptible_tries = preemptible_tries,
      JsonTools = JsonTools,
      n_cpu = n_cpu
  }

  call GetContamination {
    input:
      input_vcf = CallMt.vcf_for_haplochecker,
      sample_name = CallMt.samples,
      mean_coverage = mt_mean_coverage_array,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries,
      haplochecker_docker = haplochecker_docker
  }

  call MongoTasks_Multi.MongoM2FilterContaminationSplit as FilterContamination {
    input:
      raw_vcf = CallMt.filtered_vcf,
      raw_vcf_index = CallMt.filtered_vcf_idx,
      raw_vcf_stats = CallMt.stats,
      sample_name = GetContamination.samples,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      suffix = "",
      run_contamination = true,
      verifyBamID = verifyBamID,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_filtering_args = select_first([m2_filter_extra_args, ""]) + " --min-median-mapping-quality 0",
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      JsonTools = JsonTools,
      preemptible_tries = preemptible_tries
  }

  output {
    Array[String] samples = FilterContamination.samples
    Array[File] out_vcf = FilterContamination.filtered_vcf
    Array[File] out_vcf_index = FilterContamination.filtered_vcf_idx
    Array[File] split_vcf = FilterContamination.split_vcf
    Array[File] split_vcf_index = FilterContamination.split_vcf_index
    Array[File] nuc_vcf = select_first([CallNucHCIntegrated.full_pass_vcf, CallNucM2Integrated.full_pass_vcf])
    Array[File] nuc_vcf_index = select_first([CallNucHCIntegrated.full_pass_vcf_index, CallNucM2Integrated.full_pass_vcf_index])
    Array[File] nuc_vcf_unfiltered = select_first([CallNucHCIntegrated.filtered_vcf, CallNucM2Integrated.filtered_vcf])
    Array[File] split_nuc_vcf = select_first([CallNucHCIntegrated.split_vcf, CallNucM2Integrated.split_vcf])
    Array[File] split_nuc_vcf_index = select_first([CallNucHCIntegrated.split_vcf_index, CallNucM2Integrated.split_vcf_index])
    Array[Int] nuc_variants_pass = select_first([CallNucHCIntegrated.post_filt_vars, CallNucM2Integrated.post_filt_vars])
    Array[File] input_vcf_for_haplochecker = CallMt.vcf_for_haplochecker
    Array[File] contamination_metrics = GetContamination.contamination_file
    Array[String] major_haplogroup = GetContamination.major_hg
    Array[Float] contamination = FilterContamination.contamination
    Array[String] hasContamination = GetContamination.hasContamination
    Array[Float] contamination_major = GetContamination.major_level
    Array[Float] contamination_minor = GetContamination.minor_level
  }
}

task GetContamination {
  input {
    Array[File] input_vcf
    Array[String] sample_name
    Array[Int] mean_coverage
    # runtime    
    File JsonTools
    Int? preemptible_tries
    String haplochecker_docker
  }

  Int disk_size = ceil(size(input_vcf, "GB")) + 20
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }
  command <<<
  set -e

  mkdir out

  sampleNames=('~{sep="' '" sample_name}')
  vcfs=('~{sep="' '" input_vcf}')
  mean_covs=('~{sep="' '" mean_coverage}')

  for i in "~{d}{!sampleNames[@]}"; do

    this_basename=out/"~{d}{sampleNames[i]}"
    this_mean_cov="~{d}{mean_covs[i]}"
    this_vcf="~{d}{vcfs[i]}"

    this_vcf_nvar=$(cat "~{d}{this_vcf}" | grep ^chrM | wc -l | sed 's/^ *//g')
    echo "~{d}{sampleNames[i]} has VCF with ~{d}{this_vcf_nvar} variants for contamination."

    mkdir thisvcf
    cp "~{d}{this_vcf}" thisvcf
    PARENT_DIR=thisvcf
    
    java -jar /usr/mtdnaserver/haplocheckCLI.jar "~{d}{PARENT_DIR}"
    rm -rf thisvcf

    sed 's/\"//g' output > output-noquotes
    cp 'output-noquotes' "~{d}{this_basename}_output_noquotes"

    grep "SampleID" output-noquotes > headers
    FORMAT_ERROR="Bad contamination file format"
    if [ `awk '{print $2}' headers` != "Contamination" ]; then
      echo $FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print $6}' headers` != "HgMajor" ]; then
      echo $FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print $8}' headers` != "HgMinor" ]; then
      echo $FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
      echo $FORMAT_ERROR; exit 1
    fi
    if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
      echo $FORMAT_ERROR; exit 1
    fi

    if [ "~{d}{this_mean_cov}" -gt 0 ] && [ "~{d}{this_vcf_nvar}" -gt 0 ]; then
      grep -v "SampleID" output-noquotes > output-data
      awk -F "\t" '{print $2}' output-data > "~{d}{this_basename}.contamination.txt"
      awk -F "\t" '{print $6}' output-data > "~{d}{this_basename}.major_hg.txt"
      awk -F "\t" '{print $8}' output-data > "~{d}{this_basename}.minor_hg.txt"
      awk -F "\t" '{print $14}' output-data > "~{d}{this_basename}.mean_het_major.txt"
      awk -F "\t" '{print $15}' output-data > "~{d}{this_basename}.mean_het_minor.txt"
    else
    
      if [ "~{d}{this_mean_cov}" -eq 0 ]; then
        echo "Sample ~{d}{sampleNames[i]} has a mean coverage of 0."
      elif [ "~{d}{this_vcf_nvar}" -eq 0 ]; then
        echo "Sample ~{d}{sampleNames[i]} has no variants."
      else
        echo "Unsupported contamination exception for ~{d}{sampleNames[i]}"; exit 1
      fi
      
      echo "NO" > "~{d}{this_basename}.contamination.txt"
      echo "NONE" > "~{d}{this_basename}.major_hg.txt"
      echo "NONE" > "~{d}{this_basename}.minor_hg.txt"
      echo "0.000" > "~{d}{this_basename}.mean_het_major.txt"
      echo "0.000" > "~{d}{this_basename}.mean_het_minor.txt"
    fi

    python ~{JsonTools} \
    --path out/jsonout.json \
    --set-float major_level="$(cat ~{d}{this_basename}.mean_het_major.txt)" \
      minor_level="$(cat ~{d}{this_basename}.mean_het_minor.txt)" \
    --set samples="~{d}{sampleNames[i]}" \
      contamination_file="~{d}{this_basename}_output_noquotes" \
      hasContamination="$(cat ~{d}{this_basename}.contamination.txt)" \
      major_hg="$(cat ~{d}{this_basename}.major_hg.txt)" \
      minor_hg="$(cat ~{d}{this_basename}.minor_hg.txt)"

  done
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: haplochecker_docker
  }
  output {
    Object obj_out = read_json("out/jsonout.json")
    Array[String] samples = obj_out.samples
    Array[File] contamination_file = obj_out.contamination_file
    Array[String] hasContamination = obj_out.hasContamination
    Array[String] major_hg = obj_out.major_hg
    Array[String] minor_hg = obj_out.minor_hg
    Array[Float] major_level = obj_out.major_level
    Array[Float] minor_level = obj_out.minor_level
  }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    String base_name
    Int max_reads_per_alignment_start = 75
    String? m2_extra_args
    Boolean make_bamout = false
    Boolean nucdna = false
    Boolean disable_filters = false
    Boolean compress
    File? gatk_override
    String gatk_version

    File? mt_interval_list
    File? force_call_vcf
    File? force_call_vcf_index
    # runtime
    String? gatk_docker_override
    Int mem
    Int? preemptible_tries
    Int? n_cpu
  }

  String nucstr = if nucdna then ".nuc" else ""
  String output_vcf = base_name + nucstr + ".raw" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    echo "Extra arguments for mutect2: ""~{m2_extra_args}""$cust_interval"
    ~{"echo 'Obtaining force calls for specified VCF: '" + force_call_vcf}

    # We need to create these files regardless, even if they stay empty
    touch bamout.bam

    # Fix for DNANexus weirdness
    ~{"gatk IndexFeatureFile -I " + force_call_vcf}

    gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      ~{"-L " + mt_interval_list} \
      -O ~{output_vcf} \
      ~{true='--bam-output bamout.bam' false='' make_bamout} \
      ~{"--genotype-filtered-alleles --alleles " + force_call_vcf} \
      ~{m2_extra_args} \
      --annotation StrandBiasBySample \
      ~{true='' false='--read-filter MateOnSameContigOrNoMappedMateReadFilter' disable_filters} \
      ~{true='' false='--read-filter MateUnmappedAndUnmappedReadFilter' disable_filters} \
      ~{true='' false='--mitochondria-mode' nucdna} \
      --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
      --max-mnp-distance 0
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
      memory: machine_mem + " MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: select_first([n_cpu,2])
  }
  output {
      File raw_vcf = "~{output_vcf}"
      File raw_vcf_idx = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task Filter {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    Float? vaf_cutoff
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination 
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor
    Float? verifyBamID
     
    File? blacklisted_sites
    File? blacklisted_sites_index

    File? gatk_override
    String? gatk_docker_override
    String gatk_version

  # runtime
    Int? preemptible_tries
  }

  String output_vcf = base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20
  
  # hc_contamination will be None if hasContamination is not defined (I think) OR contamination_major not defined OR contamination_minor not defined
  String hasContamination_2 = select_first([hasContamination,"NOT FOUND"])
  Float? hc_contamination = if run_contamination && hasContamination_2 == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float hc_contamination_2 = select_first([hc_contamination, 0.0])
  Float? max_contamination = if defined(verifyBamID) then (if verifyBamID > hc_contamination_2 then verifyBamID else hc_contamination_2) else hc_contamination_2

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx2500m" FilterMutectCalls -V ~{raw_vcf} \
        -R ~{ref_fasta} \
        -O filtered.vcf \
        --stats ~{raw_vcf_stats} \
        ~{m2_extra_filtering_args} \
        --max-alt-allele-count ~{max_alt_allele_count} \
        --mitochondria-mode \
        ~{"--min-allele-fraction " + vaf_filter_threshold} \
        ~{"--f-score-beta " + f_score_beta} \
        ~{"--contamination-estimate " + max_contamination}

      ~{"gatk IndexFeatureFile -I " + blacklisted_sites}

      gatk VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        ~{"--mask-name 'blacklisted_site' --mask " + blacklisted_sites}
        
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
      memory: "4 MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_idx = "~{output_vcf_index}"
      Float? contamination = hc_contamination # now an optional output, producing UNDEFINED if not computed
  }
}