//
// Initialize options with default values.
//
def initParams(Map params) {
    params.args = params.args ?: ''
    params.publishDir = params.publishDir ?: ''
    params.publishDirMode = params.publishDirMode ?: ''
    params.publishDirOverwrite = params.publishDirMode ?: false
    return params
}

params = initParams(params)

process sambamba_markdup {
  label "process_high"
  tag "${sampleName}"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(bam)
    path index

  output:
    path outputBam
    path outputIndex

  script:
    outputBam = "${bam.simpleName}_dedup.bam"
    outputIndex = "${outputBam}.bai"
    """
    sambamba markdup --tempdir tmp -t ${task.cpus} ${params.join(' ')} ${bam} ${outputBam}
    """
}
