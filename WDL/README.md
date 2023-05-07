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