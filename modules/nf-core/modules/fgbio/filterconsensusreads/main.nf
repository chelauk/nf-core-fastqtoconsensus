process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium'

    container 'zipper.sif'
    
    input:
    tuple val(meta), path(bam)
    path fasta
	path fasta_fai
	path fasta_dict
	val (min_reads)
    val (min_base_quality)
    val (max_error_rate)

    output:
    tuple val(meta), path("*filtered.bam"), emit: filteredbam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    avail_mem = 3
    if (!task.memory) {
        log.info '[FGBIO FilterConsensusReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgbio -Xmx${avail_mem}g --compression 0 \\
        FilterConsensusReads \\
        --input $bam \\
        --ref $fasta \\
        --output /dev/stdout \\
        --min-reads $min_reads \\
        --min-base-quality $min_base_quality \\
        --max-base-error-rate $max_error_rate \\
       | samtools sort --threads $task.cpus -o ${meta.id}.cons.filtered.bam --write-index \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    avail_mem = 3
    if (!task.memory) {
        log.info '[FGBIO FilterConsensusReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "fgbio -Xmx${avail_mem}g --compression 0 FilterConsensusReads \\
        --input $bam \\
        --output /dev/stdout \\
        --min-reads $min_reads \\
        --min-base-quality $min_base_quality \\
        --max-base-error-rate $max_error_rate \\
        | samtools sort --threads $task.cpus -o ${meta.id}.cons.filtered.bam --write-index"
    
    touch  ${meta.id}.cons.filtered.bam
    touch  ${meta.id}.cons.filtered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 2.0.1 
    END_VERSIONS
    """
}
