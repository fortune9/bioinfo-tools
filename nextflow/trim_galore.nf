/*
 * This pipeline does the following things:
 * 1. trim adapter sequences
 */

// nextflow.enable.dsl=2

// default output dir
params.inFiles=null
params.outdir="trimmed"
params.keepUnpaired=false
params.lengthCutoff=20
params.rrbs=false
params.extraParams=""
params.paired=true
params.directional=true

log.info """\
         Trim fastq files with TrimGalore - N F   P I P E L I N E
         ===================================
         inFiles : ${params.inFiles}
         outdir : ${params.outdir}
         awsQueue: '${params.awsQueue}'
         awsContainer: ${params.container}
         extraParams: ${params.extraParams}
         keepUnpaired: ${params.keepUnpaired}
         lengthCutoff: ${params.lengthCutoff}
         rrbs: ${params.rrbs}
         paired: ${params.paired}
         dirctional: ${params.directional}
         """
         .stripIndent()

if(! params.inFiles) {
    exit 1, "Parameter --inFiles is needed"
}

if( (! params.directional) & (! params.rrbs) ) {
    exit 1, "set --rrbs to true, required by --directional true"
}

if(workflow.profile == 'batch') {
    if(!params.awsQueue) {
        exit 1, "For profile 'batch', parameter 'awsQueue' is required"
    }
    if(!params.container) {
        exit 1, "For profile 'batch', 'process.container' need be set"
    }
}

Channel
    .fromFilePairs( params.inFiles, checkIfExists:true) {f -> "${f.name}".replaceAll(/_S\d+_.*$/, "")}
    .set {fastqFiles}

process trim_fastq {
    tag "$sampleId" 
    echo true
    cpus 2
    memory '8 GB'
    publishDir params.outdir, mode: 'copy',
        saveAs: { filename ->
            if(filename.endsWith("_fastqc.html") ) "qc/$filename"
            else if ( filename.endsWith("_fastqc.zip") ) "qc/$filename"
            else filename
        }

    input:
    tuple val(sampleId), path(reads) from fastqFiles

    output:
    tuple val(sampleId), path("*{_val_1,_val_2,_trimmed}.fq.gz") into ch_trimmed_fastq
    path "*_trimming_report.txt" into ch_trimmed_report
    path "*_fastqc.{zip,html}" into ch_trimmed_fastqc

    script:
    rrbs = params.rrbs? "--rrbs" : ""
    keepUnpaired=params.keepUnpaired? "--retain_unpaired -r1 50 -r2 50" : ""
    lengthCutoff="--length ${params.lengthCutoff}"
    paired=params.paired? "--paired" : ""
    nonDirectional=params.directional? "" : "--non_directional"
    """
    # renaming input files
    [ ! -f  ${sampleId}_R1.fastq.gz ] && ln -s ${reads[0]} ${sampleId}_R1.fastq.gz
    [ ! -f  ${sampleId}_R2.fastq.gz ] && ln -s ${reads[1]} ${sampleId}_R2.fastq.gz
    echo Trimming fastq files of $sampleId
    trim_galore $paired ${sampleId}_R1.fastq.gz \
    ${sampleId}_R2.fastq.gz --gzip --fastqc $rrbs \
    --cores 1 $lengthCutoff $keepUnpaired \
    $nonDirectional ${params.extraParams}
    """
}

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Results are in --> $params.outdir\n" : "Oops .. something went wrong" )
}
