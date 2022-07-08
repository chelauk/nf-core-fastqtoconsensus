process FGBIO_CALLMOLECULARCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fgbio=2.0.2--hdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    val(min_reads)
    val(min_base_qual)

    output:
    tuple val(meta), path("*.consensus.bam"), emit: consensusbam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgbio \\
        CallMolecularConsensusReads \\
        --input $bam \\
        $args \\
        --min-reads $min_reads \\
        --min-input-base-quality $min_base_qual \\
        --output ${prefix}.consensus.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "fgbio \\
        CallMolecularConsensusReads \\
        --input $bam \\
        $args \\
        --min-reads $min_reads \\
        --min-input-base-quality $min_base_qual \\
        --output ${prefix}.consensus.bam"

    touch  ${prefix}.consensus.bam 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 2.0.1 
    END_VERSIONS
    """
}
