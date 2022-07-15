process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fgbio=fgbio:2.0.2--hdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val read_structure

    output:
    tuple val(meta), path("*_umi_converted.bam"), emit: umibam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir tmp

        fgbio \\
        --tmp-dir=${PWD}/tmp \\
        --compression 1 --async-io \\
        FastqToBam \\
        --input $reads \\
        --output "${prefix}_umi_converted.bam" \\
        --read-structures $read_structure \\
        --sample $meta.id \\
        --library $meta.id \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix =  task.ext.prefix ?: "${meta.id}"
    """
    echo "fgbio \\
        --tmp-dir=${PWD}/tmp \\
        --compression 1 --async-io \\
        FastqToBam \\
        -i $reads \\
        -o "${prefix}_umi_converted.bam" \\
        --read-structures $read_structure \\
        --sample $meta.id \\
        --library $meta.id \\
        $args "

    touch ${prefix}_umi_converted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 2.0.1 
    END_VERSIONS
    """
}
