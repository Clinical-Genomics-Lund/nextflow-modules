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

process chewbbaca_allelecall {
  tag "${workflow.runName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    val sampleName
    path batchInput
    path schemaDir
    path trainingFile

  output:
    val sampleName, emit: sampleName
    path 'output_dir/results_alleles.tsv', emit: calls
    //tuple val(sampleName), path("${missingLoci}"), emit: missing

  script:
    def args = task.ext.args ?: ''
    missingLoci = "chewbbaca.missingloci"
    trainingFile = trainingFile ? "--ptf ${trainingFile}" : "" 

    """
    chewie AlleleCall \\
    -i ${batchInput} \\
    ${args} \\
    --cpu ${task.cpus} \\
    --output-directory output_dir \\
    ${trainingFile} \\
    --schema-directory ${schemaDir}
    #bash parse_missing_loci.sh batch_input.list 'output_dir/*/results_alleles.tsv' ${missingLoci}
    """
}

process chewbbaca_create_batch_list {
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path maskedAssembly

  output:
    path "batch_input.list"

  script:
    output = "batch_input.list"
    """
    realpath $maskedAssembly > $output
    """
}

process chewbbaca_split_results {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    each sampleName
    path input

  output:
    tuple val(sampleName), path("${output}")

  script:
    output = "${sampleName}.chewbbaca"
    """
    head -1 ${input} > ${output}
    grep ${sampleName} ${input} >> ${output}
    """
}

process chewbbaca_split_missing_loci {
  tag "${assembly.simpleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path input

  output:
    path("${output}")

  script:
    id = "${input.simpleName}"
    output = "${id}.chewbbaca"
    """
    grep ${id} ${input} > ${output}
    """
}
