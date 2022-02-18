version 1.0

workflow AlignAndCallR1 {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File input_bam
    File input_bai
    String base_name

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

    File? gatk_override
    String? gatk_docker_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf
    Float? verifyBamID

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    #Optional runtime arguments
    Int? preemptible_tries
    Int? n_cpu
  }

  parameter_meta {
  }

  call CollectWgsMetrics {
    # get info from just the mtDNA
    input:
      input_bam = input_bam,
      input_bam_index = input_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      mt_interval_list = mt_interval_list,
      preemptible_tries = preemptible_tries
  }

  call PreProcessBam {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  # call CallNuc {
  #   input:
  #     input_bam = input_bam,
  #     input_bai = input_bai,
  #     nuc_interval_list = nuc_interval_list,
  #     ref_fasta = mt_fasta,
  #     ref_fai = mt_fasta_index,
  #     ref_dict = mt_dict,
  #     compress = compress_output_vcf,
  #     gatk_override = gatk_override,
  #     gatk_docker_override = gatk_docker_override,
  #     preemptible_tries = preemptible_tries,
  #     n_cpu = n_cpu
  # }

  call M2 as CallNucM2 {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      mt_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      nucdna = true,
      disable_filters = true,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_args = select_first([m2_extra_args, ""]) + " --minimum-allele-fraction 0.95",
      mem = M2_mem,
      preemptible_tries = preemptible_tries,
      n_cpu = n_cpu
  }

  call Filter as FilterNuc {
    input:
      raw_vcf = CallNucM2.raw_vcf,
      raw_vcf_index = CallNucM2.raw_vcf_idx,
      raw_vcf_stats = CallNucM2.stats,
      base_name = base_name + ".nuc",
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0.95,
      run_contamination = false,
      preemptible_tries = preemptible_tries
  }

  Int M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call M2 as CallMt {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      mt_interval_list = mt_interval_list,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_args = select_first([m2_extra_args, ""]),
      mem = M2_mem,
      preemptible_tries = preemptible_tries,
      n_cpu = n_cpu
  }

  call Filter as InitialFilter {
    input:
      raw_vcf = CallMt.raw_vcf,
      raw_vcf_index = CallMt.raw_vcf_idx,
      raw_vcf_stats = CallMt.stats,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      run_contamination = false,
      preemptible_tries = preemptible_tries
  }

 
  call SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override
  }

  call GetContamination {
    input:
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker,
      preemptible_tries = preemptible_tries
  }

  call Filter as FilterContamination {
    input:
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_idx,
      raw_vcf_stats = CallMt.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
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

 call MergeVcfs {
    input:
      vcf_no_filter = FilterContamination.filtered_vcf,
      vcf_to_filter = FilterNuc.filtered_vcf
 }

  output {
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_idx
    File joint_vcf = MergeVcfs.merged_vcf
    File joint_vcf_index = MergeVcfs.merged_vcf_index
    File input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    File duplicate_metrics = PreProcessBam.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    File contamination_metrics = GetContamination.contamination_file
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    Float median_coverage = CollectWgsMetrics.median_coverage
    String major_haplogroup = GetContamination.major_hg
    Float? contamination = FilterContamination.contamination
    String hasContamination = GetContamination.hasContamination
    Float contamination_major = GetContamination.major_level
    Float contamination_minor = GetContamination.minor_level
  }
}


task GetContamination {
  input {
    File input_vcf
    # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(input_vcf, "GB")) + 20

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }
  command <<<
  set -e
  PARENT_DIR="$(dirname "~{input_vcf}")"
  java -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

  sed 's/\"//g' output > output-noquotes

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

  grep -v "SampleID" output-noquotes > output-data
  awk -F "\t" '{print $2}' output-data > contamination.txt
  awk -F "\t" '{print $6}' output-data > major_hg.txt
  awk -F "\t" '{print $8}' output-data > minor_hg.txt
  awk -F "\t" '{print $14}' output-data > mean_het_major.txt
  awk -F "\t" '{print $15}' output-data > mean_het_minor.txt
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124"
  }
  output {
    File contamination_file = "output-noquotes"
    String hasContamination = read_string("contamination.txt") 
    String major_hg = read_string("major_hg.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float major_level = read_float("mean_het_major.txt")
    Float minor_level = read_float("mean_het_minor.txt")
  }
}

task PreProcessBam {
  input {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String? read_name_regex

    Int? preemptible_tries
  }

  String basename = basename(input_bam, ".bam")
  String metrics_filename = basename + ".metrics"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") * 4 + ref_size) + 20

  command <<<
    set -e
    
    java -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=~{input_bam} \
      OUTPUT=md.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false

    java -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=md.bam \
      OUTPUT=~{basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      MAX_RECORDS_IN_RAM=300000
  >>>

  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    File output_bam = "~{basename}.bam"
    File output_bam_index = "~{basename}.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task CallNucHaplotypeCaller {
  input {
      File ref_fasta
      File ref_fai
      File ref_dict

      File input_bam
      File input_bai

      Int max_reads_per_alignment_start = 75

      Boolean compress
      File gatk_override = gatk_override
      String gatk_docker_override = gatk_docker_override

      Int? preemptible_tries
      Int? n_cpu
  }

  command <<<
    set -e
    
    
  >>>
}

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    Int? read_length
    Int? coverage_cap
    File? mt_interval_list

    Int? preemptible_tries
  }

  Int read_length_for_optimization = select_first([read_length, 151])
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Collect coverage metrics"
  }
  parameter_meta {
    read_length: "Read length used for optimization only. If this is too small CollectWgsMetrics might fail. Default is 151."
  }

  command <<<
    set -e

    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      ~{"INTERVALS=" + mt_interval_list} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=metrics.txt \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length_for_optimization} \
      ~{"COVERAGE_CAP=" + coverage_cap} \
      INCLUDE_BQ_HISTOGRAM=true \
      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
      write.table(df[,"MEDIAN_COVERAGE"], "median_coverage.txt", quote=F, col.names=F, row.names=F)
    CODE
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
    Float median_coverage = read_float("median_coverage.txt")
  }
}

task MergeVcfs {
  input {
    File vcf_no_filter
    File vcf_to_filter
    String basename = basename(vcf1, ".vcf")

    File? gatk_override
    String? gatk_docker_override

    # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(vcf1, "GB") + size(vcf2, "GB")) + 20

    command<<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      gatk SelectVariants \
      -V ~{vcf_to_filter} \
      --exclude-filtered \
      -O filtered_vcf.vcf

      gatk MergeVcfs \
      I=~{vcf_no_filter} \
      I=filtered_vcf.vcf \
      O=~{basename}.mergedNuc.vcf
    >>>

    runtime {
      disks: "local-disk " + disk_size + " HDD"
      memory: "1200 MB"
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.4.1"])
      preemptible: select_first([preemptible_tries, 5])
    }
    output {
      File merged_vcf = "~{basename}.mergedNuc.vcf"
      File merged_vcf_index = "~{basename}.mergedNuc.vcf.idx"
    }
}

task LiftoverAndCombineVcfs {
  input {
    File shifted_vcf
    File vcf
    String basename = basename(shifted_vcf, ".vcf")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File shift_back_chain

    # runtime
    Int? preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(shifted_vcf, "GB") + ref_size) + 20
  String failed_vars = 'num_failed_vars.tsv'

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -jar /usr/gitc/picard.jar LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{basename}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shift_back_chain} \
      REJECT=~{basename}.rejected.vcf

    java -jar /usr/gitc/picard.jar MergeVcfs \
      I=~{basename}.shifted_back.vcf \
      I=~{vcf} \
      O=~{basename}.merged.vcf

    cat ~{basename}.rejected.vcf | grep ^chrM | wc -l > ~{failed_vars}
    >>>
    runtime {
      disks: "local-disk " + disk_size + " HDD"
      memory: "1200 MB"
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
      preemptible: select_first([preemptible_tries, 5])
    }
    output {
        # rejected_vcf should always be empty
        File rejected_vcf = "~{basename}.rejected.vcf"
        File merged_vcf = "~{basename}.merged.vcf"
        File merged_vcf_index = "~{basename}.merged.vcf.idx"
        Int number_failed = read_int("~{failed_vars}")
    }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    Int max_reads_per_alignment_start = 75
    String? m2_extra_args
    Boolean make_bamout = false
    Boolean nucdna = false
    Boolean disable_filters = false
    Boolean compress
    File? gatk_override

    File? mt_interval_list
    File? force_call_vcf
    File? force_call_vcf_index
    # runtime
    String? gatk_docker_override
    Int mem
    Int? preemptible_tries
    Int? n_cpu
  }

  String output_vcf = "raw" + if compress then ".vcf.gz" else ".vcf"
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
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.4.1"])
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

  # runtime
    Int? preemptible_tries
  }

  String output_vcf = base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20
  
  # hc_contamination will be None if hasContamination is not defined (I think) OR contamination_major not defined OR contamination_minor not defined
  Float? hc_contamination = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float? max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

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

      gatk VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        ~{"--mask-name 'blacklisted_site' --mask " + blacklisted_sites}
        
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.4.1"])
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

task MergeStats {
  input {
    File shifted_stats
    File non_shifted_stats
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  command{
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk MergeMutectStats --stats ~{shifted_stats} --stats ~{non_shifted_stats} -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.4.1"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  }
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

      gatk SelectVariants \
        -V split.vcf \
        -O splitAndPassOnly.vcf \
        --exclude-filtered
  
  }
  output {
    File vcf_for_haplochecker = "splitAndPassOnly.vcf"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.4.1"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}