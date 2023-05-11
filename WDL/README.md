# WDL input files

Several reference file and Docker image inputs are required for proper and reproducible functioning of the mtSwirl pipeline.

## Code

Supporting code can be found in `scripts`. Each of these files should be included as an input to the pipeline in the corresponding input field (see `prepopulated_inputs.json`).

## Docker images

We have made the relevant images available on Docker Hub, corresponding to the following inputs:

- `ucsc_docker`: rahulg603/ucsc_docker
- `genomes_cloud_docker`: rahulg603/genomes_cloud_docker
- `haplochecker_docker`: rahulg603/haplochecker_docker
- `gatk_samtools_docker`: rahulg603/gatk_samtools_docker

**Please note that Docker Hub places pull limits making these links unusable for large scale parallelized workflows.** We recommend using `docker pull` locally to obtain the image, then pushing the image to either GCP artifact registry (see [here](https://cloud.google.com/artifact-registry/docs/docker/store-docker-container-images) for a quickstart guide) for Terra/Cromwell-based executions or to DNANexus (see [here](https://documentation.dnanexus.com/developer/apps/dependency-management/using-docker-images)).

## Multi-sample execution on dx

There are some additional notes required for execution on DNANexus. Once you have ensured that you have proper access to the WGS files of interest on the platform, dxCompiler must first be used to upload the workflow:

```
bash run-dxcompiler-docker compile fullMitoPipeline_v2_6_Multi.wdl -project PROJECT_ID -folder DESIRED_DX_PATH -f
```

For DNANexus, we have several recommendations when building input schemas for batch job submission, building off `prepopulated_inputs_dx.json`:

1. Please do not use docker.io paths for large scale workflows, and instead upload images to your dx project (see [here](https://documentation.dnanexus.com/developer/apps/dependency-management/using-docker-images)).
2. We recommend saving all GCP hosted files (as listed in `prepopulated_inputs_dx.json`) into your own project directory on the dx platform, and suggest you change paths accordingly to the corresponding [DNANexus Link](https://documentation.dnanexus.com/developer/api/running-analyses/job-input-and-output).
3. To run this pipeline, three additional items must be added to the `prepopulated_inputs_dx.json` input used for each workflow run:
   1. `stage-common.sample_name`: an array of sample names
   2. `stage-common.wgs_aligned_input_bam_or_cram`: an array of DNANexus links corresponding to WGS CRAM files
   3. `stage-common.wgs_aligned_input_bam_or_cram_index`: an array of DNANexus links corresponding to WGS CRAI files