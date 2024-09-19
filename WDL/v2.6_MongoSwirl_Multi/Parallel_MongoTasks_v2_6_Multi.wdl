version 1.0
# This is a stupid file containing giant tasks for use with pipeline platforms
# that cannot handle a lot of parallel complexity.
# These tasks are serial combinations of more modular tasks that were located in each individual WDL.
# We cannot stream since we need to get coverage estimates.
# This version of MongoTasks contains only tasks that are capable of utilizing multi-core CPUs.
# Commented code storage:
# ~{if force_manual_download then "-I bamfile.cram --read-index bamfile.cram.crai" else "-I ~{d}{this_bam} --read-index ~{d}{this_bai}"} \
# Todo:
# - Run test on large sample to prove that this pipeline still works

task ParallelMongoSubsetBam {
  input {
    Array[File] input_bam
    Array[File] input_bai
    Array[String] sample_name
    
    File? mt_interval_list
    File? nuc_interval_list
    String? contig_name
    String? requester_pays_project
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    File? gatk_override
    String? gatk_docker_override
    String gatk_version
    File JsonTools

    Int batch_size
    Boolean force_manual_download
    Int? mem
    Int? n_cpu
    Int? preemptible_tries
  }

  Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
  Int addl_size = ceil((55 * length(sample_name)) / 1000)
  Int disk_size = ceil(ref_size) + ceil(size(input_bam,'GB')) + 20 + addl_size
  # Int machine_mem = select_first([mem, 4])
  # adjusted so we dont OOM
  # gives ~55.3 GB max worst case, ~9gb per thread

  # Int disk_size = 350
  Int machine_mem = select_first([mem, ceil(batch_size * 1.5) + 4])
  Int command_mem = ceil(1024 * 1.5)
  # overwrite this varaible for now, mem2_ssd1_v2_x16 cpu count
  Int nthreads = select_first([n_cpu, 1])-1
  String requester_pays_prefix = (if defined(requester_pays_project) then "-u " else "") + select_first([requester_pays_project, ""])
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell
  
  meta {
    description: "Subsets a whole genome bam to just Mitochondria reads in parallel"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
  }

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    mkdir out
    touch out/lockfile.lock

    process_sample() {
      local idx=$1

      bam=('~{sep="' '" input_bam}')
      bai=('~{sep="' '" input_bai}')
      sample=('~{sep="' '" sample_name}')

      local this_bam="~{d}{bam[idx]}"
      local this_bai="~{d}{bai[idx]}"
      local this_sample_t="~{d}{sample[idx]}"

      local this_sample="out/~{d}{this_sample_t}"

      echo "curr sample: ~{d}{this_sample_t}"

      ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bam} bamfile_~{d}{idx}.cram" else ""}
      ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bai} bamfile_~{d}{idx}.cram.crai" else ""}
      ~{if force_manual_download then "this_bam=bamfile_~{d}{idx}.cram" else ""}
      ~{if force_manual_download then "this_bai=bamfile_~{d}{idx}.cram.crai" else ""}

      set +e
      SAMERR=$(/usr/bin/samtools-1.9/samtools idxstats "~{d}{this_bam}" --threads ~{nthreads} 2>&1 > "~{d}{this_sample}.stats.tsv")
      thisexit=$?
      set -e

      SAMERRFAIL=$(echo $SAMERR | grep 'samtools idxstats: failed to process \".*.cram\"$' | wc -l | sed 's/^ *//g')
      echo "~{d}{this_sample_t}: samtools exited with status ~{d}{thisexit}. Match status was ~{d}{SAMERRFAIL}."

      if [ $thisexit -eq 0 ] && [ $SAMERRFAIL -eq 0 ]; then
        /usr/bin/samtools-1.9/samtools flagstat "~{d}{this_bam}" --threads ~{nthreads} > "~{d}{this_sample}.flagstat.pre.txt"

        gatk --java-options "-Xmx~{command_mem}m" PrintReads \
          ~{"-R " + ref_fasta} \
          ~{"-L " + mt_interval_list} \
          ~{"-L " + nuc_interval_list} \
          ~{"-L " + contig_name} \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          ~{"--gcs-project-for-requester-pays " + requester_pays_project} \
          -I "~{d}{this_bam}" --read-index "~{d}{this_bai}" \
          -O "~{d}{this_sample}.bam"
        
        # profiling
        cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
        mem_total=$(free | grep Mem | awk '{print $2}')
        mem_used=$(free | grep Mem | awk '{print $3}')
        mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
        echo "top output: CPU Usage: ~{d}cpu_usage%"
        echo "top output: Memory Usage: ~{d}mem_usage%"

        echo "~{d}{this_sample_t}: completed gatk. Writing to json output."
          {
            flock 200
              python ~{JsonTools} \
              --path out/jsonout.json \
              --set-int idx="~{d}{idx}" \
              --set samples="~{d}{this_sample_t}" \
                subset_bam="~{d}{this_sample}.bam" \
                subset_bai="~{d}{this_sample}.bai" \
                idxstats_metrics="~{d}{this_sample}.stats.tsv" \
                flagstat_pre_metrics="~{d}{this_sample}.flagstat.pre.txt"
          } 200>"out/lockfile.lock"
      
      elif [ $thisexit -eq 1 ] && [ $SAMERRFAIL -eq 1 ]; then
        echo "Samtools exited with status ~{d}{thisexit} and ~{d}{SAMERR}."
        echo "Thus, sample ~{d}{this_sample_t} is skipped."
      else
        echo "ERROR: samtools exited with the failure (match ~{d}{SAMERRFAIL}): ~{d}{SAMERR}."
        echo "ERROR: samtools exited with the exit status: ~{d}{thisexit}"
        exit "~{d}{thisexit}"
      fi
    }
    export -f process_sample
    
    # let's overwrite the n cpu by asking bash
    n_cpu_t=$(nproc)
    seq 0 $((~{length(input_bam)}-1)) | xargs -n 1 -P ~{batch_size} -I {} bash -c 'process_sample "$@"' _ {}

    # enforce ordering of json
    python <<EOF
  import json
  with open('out/jsonout.json', 'r') as json_file:
    data = json.load(json_file)
  idx = data['idx']
  reordered_data = {key: [value for _, value in sorted(zip(idx, data[key]))] for key in data.keys()}
  with open('out/jsonout.json', 'w') as json_file:
    json.dump(reordered_data, json_file, indent=4)
  EOF
  
  >>>
  runtime {
    memory: machine_mem + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
    # preemptible: select_first([preemptible_tries, 5])
    cpu: select_first([n_cpu, 1])
    # dx_instance_type: "mem2_ssd1_v2_x16"
  }
  
  output {
    Object obj_out = read_json("out/jsonout.json")
    Array[String] samples = obj_out.samples
    Array[File] subset_bam = obj_out.subset_bam
    Array[File] subset_bai = obj_out.subset_bai
    Array[File] idxstats_metrics = obj_out.idxstats_metrics
    Array[File] flagstat_pre_metrics = obj_out.flagstat_pre_metrics
  }
}

task ParallelMongoProcessBamAndRevert {
  input {
    Array[File] subset_bam
    Array[File] subset_bai
    Array[File] flagstat_pre_metrics
    Array[String] sample_name
    
    File? mt_interval_list
    File? nuc_interval_list
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    Boolean skip_restore_hardclips
    String? read_name_regex
    Int? read_length
    Int? coverage_cap

    File? gatk_override
    String? gatk_docker_override
    String gatk_version
    File JsonTools

    # runtime
    Int batch_size
    Int? n_cpu
    Int? preemptible_tries
  }
  Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
  Int disk_size = ceil(ref_size) + ceil(size(subset_bam,'GB')) + ceil(size(flagstat_pre_metrics,'GB')) + 20
  Int read_length_for_optimization = select_first([read_length, 151])

  Int machine_mem = ceil(select_first([batch_size, 4]) * 1.5) + 4
  Int command_mem = ceil(1024 * 1.5)

  String skip_hardclip_str = if skip_restore_hardclips then "--RESTORE_HARDCLIPS false" else ""
  Int nthreads = select_first([n_cpu,1])-1

  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell
  
  meta {
    description: "Processes a whole genome bam to just Mitochondria reads in parallel"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
  }
  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    mkdir out
    touch out/lockfile.lock

    process_sample_and_revert() {
      local idx=$1

      bam=('~{sep="' '" subset_bam}')
      bai=('~{sep="' '" subset_bai}')
      sample=('~{sep="' '" sample_name}')
      flagstat=('~{sep="' '" flagstat_pre_metrics}')

      local this_bam="~{d}{bam[idx]}"
      local this_bai="~{d}{bai[idx]}"
      local this_flagstat="~{d}{flagstat[idx]}"
      local this_sample_t="~{d}{sample[idx]}"

      local this_sample="out/~{d}{this_sample_t}"

      echo "Starting processBAM ~{d}{this_sample_t} at $(date +"%Y-%m-%d %T.%3N")."

      R --vanilla <<EOF 
        vec <- readLines("~{d}{this_flagstat}")
        titles <- c('total', 'secondary', 'supplementary', 'duplicates', 'mapped', 'paired', 'read1', 'read2', 'properly_paired', 'with_itself_and_mate_mapped', 'singletons', 'mate_diff_chr', 'mate_diff_chr_mapq_5')
        get_ele <- function(x) gregexpr('^[0-9]+',x)[[1]]
        results_vec <- as.numeric(sapply(vec, function(x) substr(x, get_ele(x)[1], get_ele(x)[1] + attr(get_ele(x), 'match.length') - 1)))
        names(results_vec) <- titles
        df <- do.call(data.frame, as.list(results_vec))
        write.table(df, sep ='\t', row.names = F, file = "~{d}{this_sample}.flagstat.txt", quote = F)
  EOF
      
      gatk --java-options "-Xmx~{command_mem}m" CollectQualityYieldMetrics \
      -I "~{d}{this_bam}" \
      ~{"-R " + ref_fasta} \
      -O "~{d}{this_sample}.yield_metrics.tmp.txt"
      cat "~{d}{this_sample}.yield_metrics.tmp.txt" | tail -n 4 | head -n 2 > "~{d}{this_sample}.yield_metrics.tmp2.txt"
      paste -d "\t" "~{d}{this_sample}.yield_metrics.tmp2.txt" "~{d}{this_sample}.flagstat.txt" > "~{d}{this_sample}.yield_metrics.txt"

      echo "Now removing mapping..."
      set +e
      gatk --java-options "-Xmx~{command_mem}m" ValidateSamFile \
        -INPUT "~{d}{this_bam}" \
        -O "~{d}{this_sample_t}.output.txt" \
        -M VERBOSE \
        -IGNORE_WARNINGS true \
        -MAX_OUTPUT 9999999
      cat "~{d}{this_sample_t}.output.txt" | \
        grep 'ERROR.*Mate not found for paired read' | \
        sed -e 's/ERROR::MATE_NOT_FOUND:Read name //g' | \
        sed -e 's/, Mate not found for paired read//g' > "~{d}{this_sample_t}.read_list.txt"
      cat "~{d}{this_sample_t}.read_list.txt" | wc -l | sed 's/^ *//g' > "~{d}{this_sample}.ct_failed.txt"

      if [[ $(tr -d "\r\n" < "~{d}{this_sample_t}.read_list.txt"|wc -c) -eq 0 ]]; then
        cp "~{d}{this_bam}" "~{d}{this_sample_t}.rescued.bam"
      else
        gatk --java-options "-Xmx~{command_mem}m" FilterSamReads \
          -I "~{d}{this_bam}" \
          -O "~{d}{this_sample_t}.rescued.bam" \
          -READ_LIST_FILE "~{d}{this_sample_t}.read_list.txt" \
          -FILTER excludeReadList
      fi
      gatk --java-options "-Xmx~{command_mem}m" RevertSam \
        -INPUT "~{d}{this_sample_t}.rescued.bam" \
        -OUTPUT_BY_READGROUP false \
        -OUTPUT "~{d}{this_sample}.unmap.bam" \
        -VALIDATION_STRINGENCY LENIENT \
        -ATTRIBUTE_TO_CLEAR FT \
        -ATTRIBUTE_TO_CLEAR CO \
        -SORT_ORDER queryname \
        -RESTORE_ORIGINAL_QUALITIES false ~{skip_hardclip_str}
      set -e
      echo "Now getting WGS metrics on the subsetted bam..."
      gatk --java-options "-Xmx~{command_mem}m" CollectWgsMetrics \
        INPUT="~{d}{this_bam}" \
        ~{"INTERVALS=" + mt_interval_list} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE="~{ref_fasta}" \
        OUTPUT="~{d}{this_sample}.wgs_metrics.txt" \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=~{read_length_for_optimization} \
        ~{"COVERAGE_CAP=" + coverage_cap} \
        INCLUDE_BQ_HISTOGRAM=true \
        THEORETICAL_SENSITIVITY_OUTPUT="~{d}{this_sample}.theoretical_sensitivity.txt"

      R --vanilla <<EOF
        df = read.table("~{d}{this_sample}.wgs_metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
        write.table(floor(df[,"MEAN_COVERAGE"]), "~{d}{this_sample}.mean_coverage.txt", quote=F, col.names=F, row.names=F)
        write.table(df[,"MEDIAN_COVERAGE"], "~{d}{this_sample}.median_coverage.txt", quote=F, col.names=F, row.names=F)
  EOF
      echo "Now preprocessing subsetted bam..."
      gatk --java-options "-Xmx~{command_mem}m" MarkDuplicates \
        INPUT="~{d}{this_bam}" \
        OUTPUT="~{d}{this_sample_t}.md.bam" \
        METRICS_FILE="~{d}{this_sample}.duplicate.metrics" \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false

      # profiling
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
      mem_total=$(free | grep Mem | awk '{print $2}')
      mem_used=$(free | grep Mem | awk '{print $3}')
      mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
      echo "top output: CPU Usage: ~{d}cpu_usage%"
      echo "top output: Memory Usage: ~{d}mem_usage%"

      echo "Now sorting md.bam for ~{d}{this_sample_t}"
      gatk --java-options "-Xmx~{command_mem}m" SortSam \
        INPUT="~{d}{this_sample_t}.md.bam" \
        OUTPUT="~{d}{this_sample}.proc.bam" \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000

      echo "~{d}{this_sample_t}: completed gatk. Writing to json output."
      {
        flock 200
        python ~{JsonTools} \
        --path out/jsonout.json \
        --set-int reads_dropped="$(cat ~{d}{this_sample}.ct_failed.txt)" \
                  idx="~{d}{idx}" \
                  mean_coverage="$(cat ~{d}{this_sample}.mean_coverage.txt)" \
        --set samples="~{d}{this_sample_t}" \
          output_bam="~{d}{this_sample}.proc.bam" \
          output_bai="~{d}{this_sample}.proc.bai" \
          unmapped_bam="~{d}{this_sample}.unmap.bam" \
          duplicate_metrics="~{d}{this_sample}.duplicate.metrics" \
          yield_metrics="~{d}{this_sample}.yield_metrics.txt"
      } 200>"out/lockfile.lock"

      echo "~{d}{this_sample_t}: WRITTEN TO JSON."
      exit_code=$?
      if [ $exit_code -ne 0 ]; then
        echo "Command failed with exit code $exit_code"
      fi
    }
    export -f process_sample_and_revert
    
    # let's overwrite the n cpu by asking bash
    n_cpu_t=$(nproc)
    echo "Processing bams across ~{d}{n_cpu_t} CPUs..."
    seq 0 $((~{length(subset_bam)}-1)) | xargs -n 1 -P ~{batch_size} -I {} bash -c 'process_sample_and_revert "$@"' _ {}
    # seq 0 $((~{length(subset_bam)}-1)) | xargs -n 1 -P ~{d}{n_cpu_t} -I {} bash -c 'process_sample_and_revert "$@"' _ {}

    echo "Finished processing BAMs, now producing output from json..."
    # call loop then read and compute mean_coverage stat to return and output for next step. if that fails, this is the place
    # enforce ordering of json
    python <<EOF
  import json
  from math import ceil
  with open('out/jsonout.json', 'r') as json_file:
    data = json.load(json_file)
  idx = data['idx']
  reordered_data = {key: [value for _, value in sorted(zip(idx, data[key]))] for key in data.keys()}
  this_max = ceil(max(reordered_data['mean_coverage']))
  with open('out/jsonout.json', 'w') as json_file:
    json.dump(reordered_data, json_file, indent=4)
  with open('this_max.txt', 'w') as f:
    f.write(str(this_max))
  EOF
  >>>
  runtime {
    memory: machine_mem + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
    cpu: select_first([n_cpu, 1])
    #mem1_ssd1_v2_x2 works well but seems to be susceptible to spotinstance interruptions
    # dx_instance_type: "mem2_ssd1_v2_x16"
  }
  output {
    Object obj_out = read_json("out/jsonout.json")
    Array[String] samples = obj_out.samples
    Array[File] output_bam = obj_out.output_bam
    Array[File] output_bai = obj_out.output_bai
    Array[File] unmapped_bam = obj_out.unmapped_bam
    Array[File] duplicate_metrics = obj_out.duplicate_metrics
    Array[File] yield_metrics = obj_out.yield_metrics
    Array[Int] reads_dropped = obj_out.reads_dropped
    Array[Int] mean_coverage = obj_out.mean_coverage
    Int max_mean_coverage = ceil(read_float("this_max.txt"))
  }
}

task ParallelMongoHC {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict

    Array[File] input_bam
    Array[File] input_bai

    Array[String] sample_name
    String suffix = ""

    Int max_reads_per_alignment_start = 75
    String? hc_extra_args
    Boolean make_bamout = false

    File? nuc_interval_list
    File? force_call_vcf
    File? force_call_vcf_index

    Boolean compress
    String gatk_version
    File? gatk_override
    String? gatk_docker_override
    Float? contamination
    File JsonTools

    Int hc_dp_lower_bound

    Int? mem
    Int batch_size
    Int? preemptible_tries
    Int? n_cpu
  }

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = ceil(batch_size * 1.5) + 4
  Int command_mem = ceil(1024 * 1.5)

  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil((size(input_bam, "GB") * 2) + ref_size) + 22

  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  command <<<
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    mkdir out
    touch out/lockfile.lock

    haplotype_caller_mt() {
      local idx=$1

      sample=('~{sep="' '" sample_name}')
      cram=('~{sep="' '" input_bam}')

      local this_sample_t="~{d}{sample[idx]}"
      local this_sample_suff="~{d}{this_sample_t}~{suffix}"
      local this_cram="~{d}{cram[idx]}"

      local this_sample="out/~{d}{this_sample_t}"
      local this_basename="~{d}{this_sample}""~{suffix}"
      local bamoutfile="~{d}{this_basename}.bamout.bam"
      
      touch "~{d}{bamoutfile}"

      if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{bamoutfile}"; else bamoutstr=""; fi
        
      gatk --java-options "-Xmx~{command_mem}m" HaplotypeCaller \
        -R ~{ref_fasta} \
        -I "~{d}{this_cram}" \
        ~{"-L " + nuc_interval_list} \
        -O "~{d}{this_basename}.raw.vcf" \
        ~{hc_extra_args} \
        -contamination ~{default="0" contamination} \
        ~{"--genotype-filtered-alleles --alleles " + force_call_vcf} \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
        --max-mnp-distance 0 \
        --annotation StrandBiasBySample \
        -G StandardAnnotation -G StandardHCAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 ~{d}{bamoutstr}

      echo "Now applying hard filters..."
      gatk --java-options "-Xmx~{command_mem}m" SelectVariants -V "~{d}{this_basename}.raw.vcf" -select-type SNP -O "~{d}{this_sample_suff}.snps.vcf"
      gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V "~{d}{this_sample_suff}.snps.vcf" \
        -R ~{ref_fasta} \
        -O "~{d}{this_sample_suff}.snps_filtered.vcf" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
        --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
        ~{'--genotype-filter-expression "DP < ' + hc_dp_lower_bound + '" --genotype-filter-name "genoDP' + hc_dp_lower_bound + '"'}

      gatk --java-options "-Xmx~{command_mem}m" SelectVariants -V "~{d}{this_basename}.raw.vcf" -select-type INDEL -O "~{d}{this_sample_suff}.indels.vcf"
      gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V "~{d}{this_sample_suff}.indels.vcf" \
        -R ~{ref_fasta} \
        -O "~{d}{this_sample_suff}.indels_filtered.vcf" \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "SOR > 10.0" --filter-name "SOR10" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
        --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
        ~{'--genotype-filter-expression "DP < ' + hc_dp_lower_bound + '" --genotype-filter-name "genoDP' + hc_dp_lower_bound + '"'}

      gatk --java-options "-Xmx~{command_mem}m" MergeVcfs -I "~{d}{this_sample_suff}.snps_filtered.vcf" -I "~{d}{this_sample_suff}.indels_filtered.vcf" -O "~{d}{this_basename}.vcf"

      # profiling
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
      mem_total=$(free | grep Mem | awk '{print $2}')
      mem_used=$(free | grep Mem | awk '{print $3}')
      mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
      echo "top output: CPU Usage: ~{d}cpu_usage%"
      echo "top output: Memory Usage: ~{d}mem_usage%"

      echo "Now filtering VCF..."
      gatk --java-options "-Xmx~{command_mem}m" SelectVariants \
        -V "~{d}{this_basename}.vcf" \
        --exclude-filtered \
        --set-filtered-gt-to-nocall \
        --exclude-non-variants \
        -O "~{d}{this_basename}.pass.vcf"

      gatk --java-options "-Xmx~{command_mem}m" CountVariants -V "~{d}{this_basename}.pass.vcf" | tail -n1 > "~{d}{this_basename}.passvars.txt"

      echo "Now splitting multi-allelics..."
      gatk --java-options "-Xmx~{command_mem}m" LeftAlignAndTrimVariants \
        -R ~{ref_fasta} \
        -V "~{d}{this_basename}.pass.vcf" \
        -O "~{d}{this_basename}.pass.split.vcf" \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac \
        --create-output-variant-index

      {
        flock 200
        python ~{JsonTools} \
        --path out/jsonout.json \
        --set-int post_filt_vars="$(cat ~{d}{this_basename}.passvars.txt)" \
                  idx="~{d}{idx}" \
        --set samples="~{d}{this_sample_t}" \
          raw_vcf="~{d}{this_basename}.raw.vcf" \
          raw_vcf_idx="~{d}{this_basename}.raw.vcf.idx" \
          output_bamOut="~{d}{bamoutfile}" \
          filtered_vcf="~{d}{this_basename}.vcf" \
          filtered_vcf_idx="~{d}{this_basename}.vcf.idx" \
          full_pass_vcf="~{d}{this_basename}.pass.vcf" \
          full_pass_vcf_index="~{d}{this_basename}.pass.vcf.idx" \
          split_vcf="~{d}{this_basename}.pass.split.vcf" \
          split_vcf_index="~{d}{this_basename}.pass.split.vcf.idx"
      } 200>"out/lockfile.lock"
    }

    export -f haplotype_caller_mt
    # let's overwrite the n cpu by asking bash
    n_cp_t=$(nproc)
    seq 0 $((~{length(input_bam)}-1)) | xargs -n 1 -P ~{batch_size} -I {} bash -c 'haplotype_caller_mt "$@"' _ {}

    # enforce ordering of json
    python <<EOF
  import json
  with open('out/jsonout.json', 'r') as json_file:
    data = json.load(json_file)
  idx = data['idx']
  reordered_data = {key: [value for _, value in sorted(zip(idx, data[key]))] for key in data.keys()}
  with open('out/jsonout.json', 'w') as json_file:
    json.dump(reordered_data, json_file, indent=4)
  EOF
  >>>

  runtime {
    docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
    memory: machine_mem + " GB"
    # disks: "local-disk " + disk_size + " HDD"
    preemptible: select_first([preemptible_tries, 5])
    cpu: select_first([n_cpu,1])
    # dx_instance_type: "mem2_ssd1_v2_x8"
  }
  output {
    Object obj_out = read_json("out/jsonout.json")
    Array[String] samples = obj_out.samples
    Array[File] raw_vcf = obj_out.raw_vcf
    Array[File] raw_vcf_idx = obj_out.raw_vcf_idx
    Array[File] output_bamOut = obj_out.output_bamOut
    Array[File] filtered_vcf = obj_out.filtered_vcf
    Array[File] filtered_vcf_idx = obj_out.filtered_vcf_idx
    Array[File] full_pass_vcf = obj_out.full_pass_vcf
    Array[File] full_pass_vcf_index = obj_out.full_pass_vcf_index
    Array[Int] post_filt_vars = obj_out.post_filt_vars
    Array[File] split_vcf = obj_out.split_vcf
    Array[File] split_vcf_index = obj_out.split_vcf_index
  }
}

task ParallelMongoAlignToMtRegShiftedAndMetrics {
  input {
    Array[File] input_bam
    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -Y $bash_ref_fasta"
    Array[String] sample_base_name
    String suffix

    File selfref_bundle

    Array[String] mt
    Array[String] mt_index
    Array[String] mt_dict

    Array[String] mt_cat
    Array[String] mt_cat_index
    Array[String] mt_cat_dict

    Array[String] mt_shifted
    Array[String] mt_shifted_index
    Array[String] mt_shifted_dict

    Array[String] mt_shifted_cat
    Array[String] mt_shifted_cat_index
    Array[String] mt_shifted_cat_dict

    Array[String] mt_interval_list

    String? read_name_regex
    File JsonTools

    Int? read_length
    Int? coverage_cap
    Int? n_cpu
    Int batch_size

    Int? preemptible_tries
  }

  Int this_cpu = select_first([n_cpu, 2])
  String this_bwa_commandline = bwa_commandline + " -t " + this_cpu
  #Float ref_size = size(mt_cat, "GB") + size(mt_cat_index, "GB")
  #Float shifted_ref_size = size(mt_shifted_cat, "GB") + size(mt_shifted_cat_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") * 4  + ceil(size(selfref_bundle, 'GB')) * 4) + 20

  Int read_length_for_optimization = select_first([read_length, 151])

  Int command_mem = 1024 * 3

  # num gigabytes
  Int machine_mem = batch_size * 3 + 4

  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  meta {
    description: "Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam."
  }
  parameter_meta {
    input_bam: "Unmapped bam"
  }
  command <<<
    export BWAVERSION=$(/usr/gitc/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
    set -o pipefail
    set -e
    tar xf "~{selfref_bundle}"
    mkdir out
    touch out/lockfile.lock

    align_to_mt_reg_shifted_metrics() {
      local idx=$1

      sample=('~{sep="' '" sample_base_name}')
      bam=('~{sep="' '" input_bam}')
      mt_intervals=('~{sep="' '" mt_interval_list}')
      mt_cat_fasta=('~{sep="' '" mt_cat}')
      mt_fasta=('~{sep="' '" mt}')
      mt_shifted_cat_fasta=('~{sep="' '" mt_shifted_cat}')
      mt_shifted_fasta=('~{sep="' '" mt_shifted}')

      local this_sample_t="~{d}{sample[idx]}"
      local this_bam="~{d}{bam[idx]}"
      local this_mt_intervals="~{d}{mt_intervals[idx]}"
      local this_mt_cat_fasta="~{d}{mt_cat_fasta[idx]}"
      local this_mt_fasta="~{d}{mt_fasta[idx]}"
      local this_mt_shifted_cat_fasta="~{d}{mt_shifted_cat_fasta[idx]}"
      local this_mt_shifted_fasta="~{d}{mt_shifted_fasta[idx]}"

      local this_sample_suff="~{d}{this_sample_t}~{suffix}"
      local this_sample=out/"~{d}{this_sample_suff}"

      local this_sample_fastq="~{d}{this_sample_suff}".fastq
      local this_sample_fastq_shifted="~{d}{this_sample_suff}".shifted.fastq
      local this_sample_bam_aligned="~{d}{this_sample_suff}".aligned.bam
      local this_sample_bam_aligned_shifted="~{d}{this_sample_suff}".shifted.aligned.bam

      local this_output_bam_basename=out/"$(basename ~{d}{this_bam} .bam).remap~{suffix}"
      
      echo "tmp files: ~{d}{this_sample_fastq}   ~{d}{this_sample_fastq_shifted}    ~{d}{this_sample_bam_aligned}    ~{d}{this_sample_bam_aligned_shifted}"

      # set the bash variable needed for the command-line
      /usr/gitc/bwa index "~{d}{this_mt_cat_fasta}"
      local bash_ref_fasta="~{d}{this_mt_cat_fasta}"

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT="~{d}{this_bam}" \
        FASTQ="~{d}{this_sample_fastq}" \
        INTERLEAVE=true \
        NON_PF=true

      /usr/gitc/~{this_bwa_commandline} "~{d}{this_sample_fastq}" - 2> >(tee "~{d}{this_output_bam_basename}.bwa.stderr.log" >&2) > "~{d}{this_sample_bam_aligned}"

      ls out/*.aligned.bam
      
      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM="~{d}{this_sample_bam_aligned}" \
        UNMAPPED_BAM="~{d}{this_bam}" \
        OUTPUT="~{d}{this_sample_suff}.mba.bam" \
        REFERENCE_SEQUENCE="~{d}{this_mt_cat_fasta}" \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION=$BWAVERSION \
        PROGRAM_GROUP_COMMAND_LINE="~{this_bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        MarkDuplicates \
        INPUT="~{d}{this_sample_suff}.mba.bam" \
        OUTPUT="~{d}{this_sample_suff}.md.bam" \
        METRICS_FILE="~{d}{this_output_bam_basename}.metrics" \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false
      echo "md BAM files:"
      ls out/*.md.bam

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        SortSam \
        INPUT="~{d}{this_sample_suff}.md.bam" \
        OUTPUT="~{d}{this_output_bam_basename}_pre_mt_filt.bam" \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000
      
      echo "md_filt BAM files:"
      ls out/*_pre_mt_filt.bam
      # now we have to subset to mito and update sequence dictionary
      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        ReorderSam \
        I="~{d}{this_output_bam_basename}_pre_mt_filt.bam" \
        O="~{d}{this_output_bam_basename}.bam" \
        REFERENCE="~{d}{this_mt_fasta}" \
        ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \
        CREATE_INDEX=true

      echo "Now starting on shifted..."
      # set the bash variable needed for the command-line
      /usr/gitc/bwa index "~{d}{this_mt_shifted_cat_fasta}"
      local bash_ref_fasta="~{d}{this_mt_shifted_cat_fasta}"
      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT="~{d}{this_bam}" \
        FASTQ="~{d}{this_sample_fastq_shifted}" \
        INTERLEAVE=true \
        NON_PF=true

      # profiling
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
      mem_total=$(free | grep Mem | awk '{print $2}')
      mem_used=$(free | grep Mem | awk '{print $3}')
      mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
      echo "top output: CPU Usage: ~{d}cpu_usage%"
      echo "top output: Memory Usage: ~{d}mem_usage%"

      /usr/gitc/~{this_bwa_commandline} "~{d}{this_sample_fastq_shifted}" - 2> >(tee "~{d}{this_output_bam_basename}.shifted.bwa.stderr.log" >&2) > "~{d}{this_sample_bam_aligned_shifted}"

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM="~{d}{this_sample_bam_aligned_shifted}" \
        UNMAPPED_BAM="~{d}{this_bam}" \
        OUTPUT="~{d}{this_sample_suff}.mba.shifted.bam" \
        REFERENCE_SEQUENCE="~{d}{this_mt_shifted_cat_fasta}" \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION=$BWAVERSION \
        PROGRAM_GROUP_COMMAND_LINE="~{this_bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        MarkDuplicates \
        INPUT="~{d}{this_sample_suff}.mba.shifted.bam" \
        OUTPUT="~{d}{this_sample_suff}.md.shifted.bam" \
        METRICS_FILE="~{d}{this_output_bam_basename}.shifted.metrics" \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false

      # profiling
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
      mem_total=$(free | grep Mem | awk '{print $2}')
      mem_used=$(free | grep Mem | awk '{print $3}')
      mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
      echo "top output: CPU Usage: ~{d}cpu_usage%"
      echo "top output: Memory Usage: ~{d}mem_usage%"

      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        SortSam \
        INPUT="~{d}{this_sample_suff}.md.shifted.bam" \
        OUTPUT="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bam" \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000

      # now we have to subset to mito and update sequence dictionary
      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        ReorderSam \
        I="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bam" \
        O="~{d}{this_output_bam_basename}.shifted.bam" \
        REFERENCE="~{d}{this_mt_shifted_fasta}" \
        ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \
        CREATE_INDEX=true

      echo "Now collecting wgs metrics..."
      java "-Xms~{command_mem}m" "-Xmx~{command_mem}m" -jar /usr/gitc/picard.jar \
        CollectWgsMetrics \
        INPUT="~{d}{this_output_bam_basename}.bam" \
        INTERVALS="~{d}{this_mt_intervals}" \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE="~{d}{this_mt_fasta}" \
        OUTPUT="~{d}{this_sample}_r2_wgs_metrics.txt" \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=~{read_length_for_optimization} \
        ~{"COVERAGE_CAP=" + coverage_cap} \
        INCLUDE_BQ_HISTOGRAM=true \
        THEORETICAL_SENSITIVITY_OUTPUT="~{d}{this_sample}_r2_wgs_theoretical_sensitivity.txt"

      R --vanilla <<EOF
        df = read.table("~{d}{this_sample}_r2_wgs_metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
        write.table(floor(df[,"MEAN_COVERAGE"]), "~{d}{this_sample}_r2_mean_coverage.txt", quote=F, col.names=F, row.names=F)
        write.table(df[,"MEDIAN_COVERAGE"], "~{d}{this_sample}_r2_median_coverage.txt", quote=F, col.names=F, row.names=F)
  EOF
      {
        flock 200
        python ~{JsonTools} \
        --path out/jsonout.json \
        --set-int mean_coverage="$(cat ~{d}{this_sample}_r2_mean_coverage.txt)" \
                  idx="~{d}{idx}" \
        --set-float median_coverage="$(cat ~{d}{this_sample}_r2_median_coverage.txt)" \
        --set samples="~{d}{this_sample_t}" \
          mt_aligned_bam="~{d}{this_output_bam_basename}.bam" \
          mt_aligned_bai="~{d}{this_output_bam_basename}.bai" \
          nuc_and_mt_aligned_bam="~{d}{this_output_bam_basename}_pre_mt_filt.bam" \
          nuc_and_mt_aligned_bai="~{d}{this_output_bam_basename}_pre_mt_filt.bai" \        
          duplicate_metrics="~{d}{this_output_bam_basename}.metrics" \
          shifted_mt_aligned_bam="~{d}{this_output_bam_basename}.shifted.bam" \
          shifted_mt_aligned_bai="~{d}{this_output_bam_basename}.shifted.bai" \        
          nuc_and_shifted_mt_aligned_bam="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bam" \
          nuc_and_shifted_mt_aligned_bai="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bai" \
          wgs_metrics="~{d}{this_sample}_r2_wgs_metrics.txt" \        
          theoretical_sensitivity="~{d}{this_sample}_r2_wgs_theoretical_sensitivity.txt"
        } 200>"out/lockfile.lock"
    }

    export -f align_to_mt_reg_shifted_metrics
    # let's overwrite the n cpu by asking bash
    n_cp_t=$(nproc)
    seq 0 $((~{length(input_bam)}-1)) | xargs -n 1 -P ~{batch_size} -I {} bash -c 'align_to_mt_reg_shifted_metrics "$@"' _ {}

    # enforce ordering of json and compute maximum mean coverage
    python <<EOF
  import json
  from math import ceil
  with open('out/jsonout.json', 'r') as json_file:
    data = json.load(json_file)
  idx = data['idx']
  reordered_data = {key: [value for _, value in sorted(zip(idx, data[key]))] for key in data.keys()}
  this_max = ceil(max(reordered_data['mean_coverage']))
  with open('out/jsonout.json', 'w') as json_file:
    json.dump(reordered_data, json_file, indent=4)
  with open('this_max_r2.txt', 'w') as f:
    f.write(str(this_max))
  EOF
  echo "Finished processing all BAMs! Delocalizing..."
  sleep 2
  >>>
  runtime {
    # preemptible: select_first([preemptible_tries, 5])
    # memory: "55 GB"
    # cpu: this_cpu
    # disks: "local-disk " + disk_size + " SSD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    cpu : this_cpu
    memory : machine_mem + " GB"
    # dx_instance_type: "mem2_ssd1_v2_x16"
  }
  output {
    Object obj_out = read_json('out/jsonout.json')

    Array[String] samples = obj_out.samples
    Array[File] mt_aligned_bam = obj_out.mt_aligned_bam
    Array[File] mt_aligned_bai = obj_out.mt_aligned_bai
    Array[File] nuc_and_mt_aligned_bam = obj_out.nuc_and_mt_aligned_bam
    Array[File] duplicate_metrics = obj_out.duplicate_metrics

    Array[File] shifted_mt_aligned_bam = obj_out.shifted_mt_aligned_bam
    Array[File] nuc_and_shifted_mt_aligned_bam = obj_out.nuc_and_shifted_mt_aligned_bam

    Array[File] wgs_metrics = obj_out.wgs_metrics
    Array[File] theoretical_sensitivity = obj_out.theoretical_sensitivity
    Array[Int] mean_coverage = obj_out.mean_coverage
    Array[Float] median_coverage = obj_out.median_coverage

    Int max_mean_coverage = ceil(read_float("this_max_r2.txt"))
  }
}

task ParallelMongoCallMtAndShifted {
  input {
    Array[String] sample_base_name
    File selfref_bundle
    Array[String] mt_self
    Array[String] mt_self_index
    Array[String] mt_self_dict
    Array[File] input_bam
    Array[File] input_bai
    Array[String] mt_interval_list
    Array[String] force_call_vcf
    Array[String] force_call_vcf_idx

    Array[String] shifted_mt_self
    Array[String] shifted_mt_self_index
    Array[String] shifted_mt_self_dict
    Array[File] shifted_input_bam
    #Array[File] shifted_input_bai
    Array[String] shifted_mt_interval_list
    Array[String] shifted_force_call_vcf
    Array[String] shifted_force_call_vcf_idx

    String? m2_extra_args
    String? shifted_m2_extra_args

    Int max_reads_per_alignment_start = 75
    Boolean make_bamout = false
    Boolean compress
    String suffix

    # runtime
    File? gatk_override
    String gatk_version
    String? gatk_docker_override
    File JsonTools
    Int batch_size
    Int? preemptible_tries
    Int? n_cpu
  }

  Float ref_size = size(selfref_bundle, 'GB')
  #Float ref_size = size(mt_self, "GB") + size(mt_self_index, "GB") + size(shifted_mt_self, "GB") + size(shifted_mt_self_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + size(shifted_input_bam, "GB") + ref_size) + 20

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = ceil(batch_size * 0.25) + 2
  Int command_mem = ceil(1024 * 0.25)

  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  meta {
    description: "Mutect2 for calling Snps and Indels; runs on both shifted and nonshifted"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
    set -e

    tar xf "~{selfref_bundle}"
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    echo "Extra arguments for mutect2: ""~{m2_extra_args}""$cust_interval"

    mkdir out
    touch out/lockfile.lock

    call_mt_and_shifted() {
      local idx=$1

      sampleNames=('~{sep="' '" sample_base_name}')
      bams=('~{sep="' '" input_bam}')
      intervals=('~{sep="' '" mt_interval_list}')
      fastas=('~{sep="' '" mt_self}')
      force_call_self=('~{sep="' '" force_call_vcf}')
      shifted_bams=('~{sep="' '" shifted_input_bam}')
      shifted_intervals=('~{sep="' '" shifted_mt_interval_list}')
      shifted_fastas=('~{sep="' '" shifted_mt_self}')
      shifted_force_call_self=('~{sep="' '" shifted_force_call_vcf}')

      local this_sample_t="~{d}{sampleNames[idx]}"
      local this_sample=out/"~{d}{this_sample_t}~{suffix}"

      local this_bam="~{d}{bams[idx]}"
      local this_noncontrol="~{d}{intervals[idx]}"
      local this_force_vcf="~{d}{force_call_self[idx]}"
      local this_self_fasta="~{d}{fastas[idx]}"
      local this_shifted_bam="~{d}{shifted_bams[idx]}"
      local this_control="~{d}{shifted_intervals[idx]}"
      local this_shifted_force_vcf="~{d}{shifted_force_call_self[idx]}"
      local this_self_shifted_fasta="~{d}{shifted_fastas[idx]}"

      touch "~{d}{this_sample}.bamout.bam"
      touch "~{d}{this_sample}.shifted.bamout.bam"

      samtools index "~{d}{this_shifted_bam}"

      if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{this_sample}.bamout.bam"; else bamoutstr=""; fi
      if [[ ~{make_bamout} == 'true' ]]; then shiftedbamoutstr="--bam-output ~{d}{this_sample}.shifted.bamout.bam"; else shiftedbamoutstr=""; fi

      echo "Obtaining force calls for specified VCF: ~{d}{this_force_vcf}"

      # Fix for DNANexus weirdness
      gatk IndexFeatureFile -I "~{d}{this_force_vcf}"

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
        -R "~{d}{this_self_fasta}" \
        -I "~{d}{this_bam}" \
        -L "~{d}{this_noncontrol}" \
        -O "~{d}{this_sample}.raw.vcf" \
        --genotype-filtered-alleles \
        --alleles "~{d}{this_force_vcf}" \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
        --max-mnp-distance 0 ~{d}{bamoutstr}

      echo "Obtaining force calls for specified VCF: ~{d}{this_shifted_force_vcf}"

      # Fix for DNANexus weirdness
      gatk IndexFeatureFile -I "~{d}{this_shifted_force_vcf}"

      # profiling
      cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2 + $4}')
      mem_total=$(free | grep Mem | awk '{print $2}')
      mem_used=$(free | grep Mem | awk '{print $3}')
      mem_usage=$(echo "scale=2; ~{d}mem_used/~{d}mem_total*100" | bc)
      echo "top output: CPU Usage: ~{d}cpu_usage%"
      echo "top output: Memory Usage: ~{d}mem_usage%"

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
        -R "~{d}{this_self_shifted_fasta}" \
        -I "~{d}{this_shifted_bam}" \
        -L "~{d}{this_control}" \
        -O "~{d}{this_sample}.shifted.raw.vcf" \
        --genotype-filtered-alleles \
        --alleles "~{d}{this_shifted_force_vcf}" \
        ~{shifted_m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
        --max-mnp-distance 0 ~{d}{shiftedbamoutstr}

      {
        flock 200
        python ~{JsonTools} \
        --path out/jsonout.json \
        --set-int idx="~{d}{idx}" \
        --set samples="~{d}{this_sample_t}" \
          raw_vcf="~{d}{this_sample}.raw.vcf" \
          raw_vcf_idx="~{d}{this_sample}.raw.vcf.idx" \
          stats="~{d}{this_sample}.raw.vcf.stats" \
          output_bamOut="~{d}{this_sample}.bamout.bam" \        
          shifted_raw_vcf="~{d}{this_sample}.shifted.raw.vcf" \
          shifted_raw_vcf_idx="~{d}{this_sample}.shifted.raw.vcf.idx" \
          shifted_stats="~{d}{this_sample}.shifted.raw.vcf.stats" \
          shifted_output_bamOut="~{d}{this_sample}.shifted.bamout.bam"
      } 200>"out/lockfile.lock"
    }

    export -f call_mt_and_shifted
    # n_cpu_t=$(nproc)
    echo "len of input bam: $((~{length(input_bam)}-1))"
    seq 0 $((~{length(input_bam)}-1)) | xargs -n 1 -P ~{batch_size} -I {} bash -c 'call_mt_and_shifted "$@"' _ {}

    # enforce ordering of json
    python <<EOF
  import json
  with open('out/jsonout.json', 'r') as json_file:
    data = json.load(json_file)
  idx = data['idx']
  reordered_data = {key: [value for _, value in sorted(zip(idx, data[key]))] for key in data.keys()}
  with open('out/jsonout.json', 'w') as json_file:
    json.dump(reordered_data, json_file, indent=4)
  EOF
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
      memory: machine_mem + " GB"
      disks: "local-disk " + disk_size + " SSD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: select_first([n_cpu,2])
  }
  output {
    Object obj_out = read_json('out/jsonout.json')
    Array[File] raw_vcf = obj_out.raw_vcf
    Array[File] raw_vcf_idx = obj_out.raw_vcf_idx
    Array[File] stats = obj_out.stats
    Array[File] output_bamOut = obj_out.output_bamOut

    Array[File] shifted_raw_vcf = obj_out.shifted_raw_vcf
    Array[File] shifted_raw_vcf_idx = obj_out.shifted_raw_vcf_idx
    Array[File] shifted_stats = obj_out.shifted_stats
    Array[File] shifted_output_bamOut = obj_out.shifted_output_bamOut
  }
}
