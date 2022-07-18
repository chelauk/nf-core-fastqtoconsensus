# nf-core/fastqtoconsensus: Documentation
The nf-core/fastqtoconsensus documentation is split into the following pages:
- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.

- [Output](output.md)

  - An overview of the different results produced by the pipeline and how to interpret them.
The pipeline processes UMI fastqs through to bam
```mermaid
flowchart TD
step1(umi-fastqs)
step1 -->|QC|step2(fastqc)
step1 -->step3(fastq -> ubam fgbio:FastqToBam)
step3 -->|tagged unaligned bam|step4(ubam,ref -> aligned bam samtools:fastq,bwa:mem,fgbio:ZipperBams)
step4 -->|aligned and tagged bam|step5(tagged bam -> grouped bam fgio:GroupReadsByUmi)
step5 -->|grouped umis|step6(grouped bam -> consensus bam fgbio:CallMolecularConsensus)
step6 -->|consensus bam|step7(filter bam fgbio:FilterConsensusReads)
step7-->step8(filtered bam file)
```
