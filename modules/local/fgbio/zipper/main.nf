process FGBIO_ZIPPER {
    tag "$meta.id"
    label 'process_high'

    container 'zipper.sif'

    input:
    tuple val(meta), path(bam)
    path index
	path fasta
    path fasta_fai
	path dict

    output:
    tuple val(meta), path("*mapped.bam"), emit: zipperbam
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    avail_mem = 3
    if (!task.memory) {
        log.info '[FGBIO ZipperBams] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')
    samtools sort --no-PG -n $bam -o sorted.bam
    samtools fastq $bam \\
    | bwa mem -t $task.cpus -p -Y -K 150000000 \$INDEX - \\
    | samtools sort --no-PG -n \\
    | fgbio -Xmx${avail_mem}g --compression 1 ZipperBams \\
    --unmapped sorted.bam \\
    --ref $fasta \\
    --output ${prefix}.mapped.bam \\
    --tags-to-reverse Consensus \\
    --tags-to-revcomp Consensus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    avail_mem = 3
    if (!task.memory) {
        log.info '[FGBIO ZipperBams] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e  "samtools fastq $bam | \n
    bwa mem -t $task.cpus() -p -K 15000000 -Y $fasta - | \n
    fgbio -Xmx${avail_mem}g --compression 1 --async-io ZipperBams \n
    --unmapped $bam \n
    --ref fasta \n
    --output ${prefix}.mapped.bam \n
    --tags-to-reverse Consensus \n
    --tags-to-revcomp Consensus "
    touch ${prefix}.mapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 2.0.1
    END_VERSIONS
    """
}
