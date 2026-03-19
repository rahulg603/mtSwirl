version 1.0

import "https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/LiftoverTools_v2_6_Multi.wdl" as LiftoverTools_Multi
import "https://raw.githubusercontent.com/rahulg603/mtSwirl/master/WDL/v2.6_MongoSwirl_Multi/MongoTasks_v2_6_Multi.wdl" as MongoTasks_Multi

workflow ProduceSelfReferenceFiles {
  # Running bcftools consensus when an interval file was used upstream to subset the genome
  #  is problematic because the consensus now is a fasta file with blocks per interval,
  #  indexed to 1 as the first base in each interval (so chr1 65000-65050 NUMT -> NUMT 1-51).
  #
  # See documentation for ProduceSelfReference for how we resolve this -- namely,
  #  we manually change the name for each interval to chr#:start-end from the interval file definition.
  #  This results in bcftools consensus producing files that are properly indexed to genomic coordinates,
  #  critical for getting variants in the VCF to apply to the fasta.
  #
  # However this also produces a problem -- the rest of our machinery assumes that nucDNA segments are
  #  1-indexed. This is important because it is not clear that the GATK/Picard tools
  #  can infer positional subsets of genomes based on a FASTA header, and we don't want to carry around
  #  the entire human genome sequence so we want to produce FASTAs with just subsets of the genome. Thus we
  #  revert the renaming within ProduceSelfReference after variants are successfully applied. 
  #
  # We then use CreateSpanIntervalsWithDict to "lift over" the nucDNA and mtDNA intervals files 
  #  by producing the exactly correct intervals that cover each entire region in the consensus FASTA dict (1-indexed).
  #
  # Finally we have to fix the chain files outputted from ProduceSelfReference which
  #  still use genomic coordinates. To do this we use MoveChainToZero, which renames
  #  and shifts each chain file block backwards such that it starts at 0 (correct for 1-indexed regions).
  #
  # In sum here we produce reference files such that variants are applied appropriately
  #  even when we have multiple intervals that do not all start at position 1 of the chromosome.
  #  Outputted files produce "contigs" for each interval, 1-indexed to the start of the interval and not the actual chromosome.
  #  This part is only supported for nucDNA intervals. The machinery is ready for mtDNA intervals but not tested -- we do not support analysis on subsets of the mtDNA.
  meta {
    description: "Produces all relevant self-reference files for version 2.2 of MitochondriaPipeline."
  }

  input {
    Array[String] sample_name
    Array[File] nuc_variants
    Array[File] mtdna_variants
    String suffix

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_interval_list
    File non_control_region_interval_list

    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File nuc_interval_list
    String reference_name = "reference"

    File blacklisted_sites
    File blacklisted_sites_index

    Int n_shift = 8000

    Boolean compute_numt_coverage
    File FaRenamingScript
    File CheckVariantBoundsScript
    File CheckHomOverlapScript
    File JsonTools

    String genomes_cloud_docker
    String ucsc_docker

    #Optional runtime arguments
    Int? preemptible_tries
  }

  call MongoTasks_Multi.MongoProduceSelfReference as ProduceSelfReference {
    input:
      input_mt_vcf = mtdna_variants,
      input_nuc_vcf = nuc_variants,
      sample_name = sample_name,
      suffix = suffix,
      
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      nuc_interval_list = nuc_interval_list,
      mt_ref_fasta = mt_fasta,
      mt_ref_fasta_index = mt_fasta_index,
      mt_interval_list = mt_interval_list,
      fa_renaming_script = FaRenamingScript,
      variant_bounds_script = CheckVariantBoundsScript,
      check_hom_overlap_script = CheckHomOverlapScript,
      non_control_region_interval_list = non_control_region_interval_list,
      n_shift = n_shift,
      intertext = '',
      JsonTools = JsonTools,
      genomes_cloud_docker = genomes_cloud_docker,
      preemptible_tries = preemptible_tries
  }

  call MongoTasks_Multi.MongoChainSwapLiftoverBed as ChainSwapLiftoverBed {
    input:
      selfref_bundle = ProduceSelfReference.produce_self_ref_bundle,
      source_chain = ProduceSelfReference.grch38_to_self_chain,
      input_target_name = ProduceSelfReference.samples,

      input_source_name = reference_name,
      input_bed = blacklisted_sites,
      input_bed_index = blacklisted_sites_index,
      JsonTools = JsonTools,
      ucsc_docker = ucsc_docker,
      preemptible_tries = preemptible_tries
  }

  if (compute_numt_coverage) {
    call CreateSpanIntervalsWithDict as LiftOverNucReference {
      input:
        selfref_bundle = ProduceSelfReference.produce_self_ref_bundle,
        input_intervals = nuc_interval_list,
        target_dict = ProduceSelfReference.self_cat_dict,
        sample_name = ProduceSelfReference.samples,
        intertext = '.nuc',
        JsonTools = JsonTools,
        preemptible_tries = preemptible_tries
    }

    call CreateSpanIntervalsWithDict as LiftOverNucReferenceShifted {
      input:
        selfref_bundle = ProduceSelfReference.produce_self_ref_bundle,
        input_intervals = nuc_interval_list,
        target_dict = ProduceSelfReference.self_shifted_cat_dict,
        sample_name = ProduceSelfReference.samples,
        intertext = '.nuc.shifted',
        JsonTools = JsonTools,
        preemptible_tries = preemptible_tries
    }

    call LiftoverTools_Multi.ChainSwap as SelfToRefNucLiftoverChain {
      input:
        selfref_bundle = ProduceSelfReference.produce_self_ref_bundle,
        source_chain = ProduceSelfReference.grch38_to_self_nuc_chain,
        input_target_name = ProduceSelfReference.samples,
        input_source_name = reference_name,
        target_suffix = "_nuc",
        JsonTools = JsonTools,
        ucsc_docker = ucsc_docker,
        preemptible_tries = preemptible_tries
    }

    call MoveChainToZero {
      input:
        source_chain = SelfToRefNucLiftoverChain.chain,
        sample_name = SelfToRefNucLiftoverChain.samples,
        ref_intervals = nuc_interval_list,
        JsonTools = JsonTools,
        preemptible_tries = preemptible_tries
    }
  }

  output {
    Array[String] samples = ChainSwapLiftoverBed.samples
    File selfref_bundle = ProduceSelfReference.produce_self_ref_bundle
    Array[String] mt_self = ProduceSelfReference.self_fasta
    Array[String] mt_self_index = ProduceSelfReference.self_fasta_index
    Array[String] mt_self_dict = ProduceSelfReference.self_dict

    Array[String] mt_shifted_self = ProduceSelfReference.self_shifted_fasta
    Array[String] mt_shifted_self_index = ProduceSelfReference.self_shifted_fasta_index
    Array[String] mt_shifted_self_dict = ProduceSelfReference.self_shifted_dict

    Array[String] ref_to_self_chain = ProduceSelfReference.grch38_to_self_chain
    Array[File] self_to_ref_chain = ChainSwapLiftoverBed.chain
    Array[String] self_shift_back_chain = ProduceSelfReference.shift_back_chain
    Array[File]? nuc_self_to_ref_chain = MoveChainToZero.chain

    Array[String] mt_andNuc_self = ProduceSelfReference.self_cat_fasta
    Array[String] mt_andNuc_self_index = ProduceSelfReference.self_cat_fasta_index
    Array[String] mt_andNuc_self_dict = ProduceSelfReference.self_cat_dict

    Array[String] mt_andNuc_shifted_self = ProduceSelfReference.self_shifted_cat_fasta
    Array[String] mt_andNuc_shifted_self_index = ProduceSelfReference.self_shifted_cat_fasta_index
    Array[String] mt_andNuc_shifted_self_dict = ProduceSelfReference.self_shifted_cat_dict

    Array[String] ref_homoplasmies_vcf = ProduceSelfReference.filtered_vcf_ref_coord
    Array[String] force_call_vcf = ProduceSelfReference.reversed_hom_vcf
    Array[String] force_call_vcf_idx = ProduceSelfReference.reversed_hom_vcf_idx
    Array[String] force_call_vcf_filters = ProduceSelfReference.reversed_hom_filters_vcf
    Array[String] force_call_vcf_filters_idx = ProduceSelfReference.reversed_hom_filters_vcf_idx
    Array[String] force_call_vcf_shifted = ProduceSelfReference.reversed_hom_vcf_shifted
    Array[String] force_call_vcf_shifted_idx = ProduceSelfReference.reversed_hom_vcf_shifted_idx

    Array[String] mt_interval_list_self = ProduceSelfReference.lifted_mt_intervals
    Array[File]? nuc_interval_list_self = LiftOverNucReference.lifted_intervals
    Array[File]? nuc_interval_list_shifted_self = LiftOverNucReferenceShifted.lifted_intervals
    Array[File] blacklisted_sites_self = ChainSwapLiftoverBed.transformed_bed
    Array[File] blacklisted_sites_index_self = ChainSwapLiftoverBed.transformed_bed_index
    Array[String] non_control_interval_self = ProduceSelfReference.lifted_noncontrol_intervals
    Array[String] control_shifted_self = ProduceSelfReference.lifted_control_intervals

    Array[Int] nuc_variants_dropped = ProduceSelfReference.nuc_variants_dropped
    Array[Int] mtdna_consensus_overlaps = ProduceSelfReference.mtdna_consensus_overlaps
    Array[Int] nuc_consensus_overlaps = ProduceSelfReference.nuc_consensus_overlaps
  }
}

task CreateSpanIntervalsWithDict {
  input {
    File selfref_bundle
    Array[String] sample_name
    File input_intervals
    Array[String] target_dict
    String intertext
    File JsonTools

    Int? preemptible_tries
  }

  Float ref_size = size(selfref_bundle, "GB")
  Int disk_size = ceil(size(input_intervals, "GB")*2 + ref_size)*2 + 10
  String base_intervals = basename(input_intervals, ".interval_list") + intertext
  String suffix_intervals = ".SelfRefLiftover.interval_list"
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  command <<<
    set -e

    tar xf "~{selfref_bundle}"
    
    sampleNames=('~{sep="' '" sample_name}')
    dicts=('~{sep="' '" target_dict}')
    mkdir out

    for i in "~{d}{!sampleNames[@]}"; do

      this_samp_name="~{d}{sampleNames[i]}"
      this_dict="~{d}{dicts[i]}"

      R --vanilla <<CODE
        intervals <- readLines('~{input_intervals}')
        intervals <- intervals[grep('^@', intervals, invert=T)]
        interval_names <- sapply(strsplit(intervals, '\t'),function(x)x[5])
        new_header <- readLines("~{d}{this_dict}")
        lens <- sapply(interval_names, function(x)as.numeric(gsub('LN:','',strsplit(new_header[grep(paste0('^@SQ\tSN:',x), new_header)[1]], '\t')[[1]][3])))
        if(any(is.na(lens))) stop('ERROR: Some NUMT intervals were not found in the mt_andNuc sequence dictionary.')
        new_intervals <- c(new_header, paste(interval_names, 1, lens, '+', interval_names, sep='\t'))
        writeLines(new_intervals, 'out/~{base_intervals}.~{d}{this_samp_name}~{suffix_intervals}')
    CODE

      python ~{JsonTools} \
      --path out/jsonout.json \
      --set samples="~{d}{sampleNames[i]}" \
        lifted_intervals="out/~{base_intervals}.~{d}{this_samp_name}~{suffix_intervals}"

    done
  >>>

  output {
    Object obj_out = read_json("out/jsonout.json")
    Array[String] samples = obj_out.samples
    Array[File] lifted_intervals = obj_out.lifted_intervals
  }

  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "500 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
}

task MoveChainToZero {
  input {
    Array[String] sample_name
    Array[File] source_chain
    File ref_intervals
    File JsonTools

    Int? preemptible_tries
  }

  Int disk_size = ceil(size(source_chain, "GB")*2.5)
  String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

  command <<<
    set -e
    sampleNames=('~{sep="' '" sample_name}')
    chains=('~{sep="' '" source_chain}')
    mkdir out

    for i in "~{d}{!sampleNames[@]}"; do

      this_samp_name="~{d}{sampleNames[i]}"
      this_chain="~{d}{chains[i]}"
      this_output="out/$(basename ~{d}{this_chain} .chain).toZero.chain"

      R --vanilla <<CODE
        chain <- readLines("~{d}{this_chain}")
        intervals <- readLines('~{ref_intervals}')
        intervals <- intervals[grep('^@', intervals, invert=T)]
        intervals_split <- strsplit(intervals, '\t')
        target_names <- sapply(intervals_split, function(x)x[5])
        names(target_names) <- sapply(intervals_split, function(x)paste0(x[1], ':', as.numeric(x[2])-1, '-', x[3]))
        to_edit <- grep('^chain', chain)
        this_sep <- ' '
        new_headers <- sapply(chain[to_edit], function(x) {
          this_header <- strsplit(x, this_sep)[[1]]
          this_search <- paste0(this_header[8], ':', this_header[11], '-', this_header[12])
          if (!this_search %in% names(target_names)) stop('ERROR: chain file has a segment not in intervals list.')
          this_name <- target_names[this_search]
          this_header[c(3, 8)] <- this_name
          this_header[c(4, 6, 7)] <- as.numeric(this_header[c(4, 6, 7)]) - as.numeric(this_header[6])
          this_header[c(9, 11, 12)] <- as.numeric(this_header[c(9, 11, 12)]) - as.numeric(this_header[11])
          return(paste0(this_header, collapse=this_sep))
        })
        new_chain <- chain
        new_chain[to_edit] <- new_headers
        writeLines(new_chain, "~{d}{this_output}")
    CODE

      python ~{JsonTools} \
      --path out/jsonout.json \
      --set samples="~{d}{sampleNames[i]}" \
        chain="~{d}{this_output}"

    done
  >>>

  output {
    Object obj_out = read_json('out/jsonout.json')
    Array[String] samples = obj_out.samples
    Array[File] chain = obj_out.chain
  }

  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "1 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
}