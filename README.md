# ![nf-core/fastqtoconsensus](docs/images/nf-core-fastqtoconsensus_logo_light.png#gh-light-mode-only) ![nf-core/fastqtoconsensus](docs/images/nf-core-fastqtoconsensus_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/fastqtoconsensus/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/fastqtoconsensus/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/fastqtoconsensus/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/fastqtoconsensus/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/fastqtoconsensus/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/fastqtoconsensus)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23fastqtoconsensus-4A154B?logo=slack)](https://nfcore.slack.com/channels/fastqtoconsensus)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/fastqtoconsensus** is a bioinformatics best-practice analysis pipeline for Fastq -> Filtered Consensus Pipeline.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/fastqtoconsensus/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3. Create an tagged unaligned bam ([`fgbio FastqToBam`](http://fulcrumgenomics.github.io/fgbio/))
4. Align fastq and merge with tagged bam ([`fgbio Zipperbams`](http://fulcrumgenomics.github.io/fgbio/) ; [`samtools fastq`](http://www.htslib.org/) ; [`bwa mem`](https://github.com/lh3/bwa))
5. Group reads by umi ([`fgbio GroupReadsByUmi`](http://fulcrumgenomics.github.io/fgbio/))
6. Call consensus ([`fgbio CallMollecularConsensus`](http://fulcrumgenomics.github.io/fgbio/))
7. Filter consensus ([`fgbio FilterConsensusReads`](http://fulcrumgenomics.github.io/fgbio/))
## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Clone the pipeline

3. Clone the pipeline in your scratch folder

    ```console
    git clone https://github.com/chelauk/nf-core-demultiplex-methylation.git
    ```


4. Edit your .bashrc file to set the following variables:

   <pre><lang ="bash"><code>
   # Set all the Singularity cache dirs to Scratch
   export SINGULARITY_CACHEDIR=<b>/your/selected/scratch/folder/singularity_imgs</b>
   export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
   export SINGULARITY_LOCALCACHEDIR=$SINGULARITY_CACHEDIR/localcache
   export SINGULARITY_PULLFOLDER=$SINGULARITY_CACHEDIR/pull
   # match the NXF_SINGULARITY_CACHEDIR
   export NXF_SINGULARITY_CACHEDIR=<b>/your/selected/scratch/folder/singularity_imgs</b>
   </code></pre>

5. Start running your own analysis
   edit a sbatch script <code>runNextflow.sh</code>

    <pre><lang ="bash"><code>
    #!/bin/bash -l
    #SBATCH --job-name=demultiplex
    #SBATCH --output=nextflow_out.txt
    #SBATCH --partition=master-worker
    #SBATCH --ntasks=1
    #SBATCH --time=120:00:00

    nextflow run <b>/location/of/your/nextflow_pipelines/nf-core-fastqtoconsensus</b> \
    --input input.csv \
    -profile slurm,singularity \
    -c local.config \
    -resume
    </code></pre>

6. Start your sbatch job:

   ```console
   sbatch runNextflow.sh
   ````

7. local.config

    Adjust your local config file to match requirements.
    parameters can be set for individual processes or processes can be grouped with labels

    <pre><lang ="bash"><code>
    process {
      executor = 'slurm'
      errorStrategy = {task.exitStatus in [143,137,104,134,139,255] ? 'retry' : 'finish'}
      maxErrors = '-1'
      maxRetries = 5<br>
      withLabel:process_high {
        memory = 64.GB
        cpus   = 24
        time   = 48.h
      }
      withLabel:process_low {
        cpus   = 1
        memory = 8.GB
        time   = 2.h
      }
      withLabel:process_long {
        memory = 16.GB
        cpus = 1
        time = 72.h
      }
      withLabel:process_medium {
        memory = 16.GB
        time = 8.h
      }
      withName:ALIGN {
        cpus = 48
        memory = 384.GB
        time = 48.h
      }
    }
</pre></code>


   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run nf-core/fastqtoconsensus --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/fastqtoconsensus pipeline comes with documentation about the pipeline [usage](https://nf-co.re/fastqtoconsensus/usage), [parameters](https://nf-co.re/fastqtoconsensus/parameters) and [output](https://nf-co.re/fastqtoconsensus/output).

## Credits

nf-core/fastqtoconsensus was originally written by Chela James.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#fastqtoconsensus` channel](https://nfcore.slack.com/channels/fastqtoconsensus) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/fastqtoconsensus for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
