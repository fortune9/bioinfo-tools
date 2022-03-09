/*
 * This pipeline does the following things:
 * 1. align reads to genome
 */

// nextflow.enable.dsl=2

prog="bismark.nf"

// default output dir
params.indexFolder=null
params.inFiles=null
params.outdir="bam"
params.directional=false
params.extra=""


log.info """\
         Call Methylation from Bam files - N F   P I P E L I N E
         ===================================
         indexFolder : ${params.indexFolder}
         inFiles : ${params.inFiles}
         outdir : ${params.outdir}
         directional: ${params.directional}
         params.extra: ${params.extra}
         awsQueue: '${params.awsQueue}'
         awsContainer: ${params.container}

         e.g.: $prog --indexFolder s3://to/bismark/index \
            --inFiles  s3://to/*_R{1,2}.fq.gz \
            --outdir s3://path/to/save
         """
         .stripIndent()

if(params.indexFolder==null | params.inFiles==null ) {
    exit 1, "parameter --indexFolder and --inFiles are both required"
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
    .fromFilePairs( params.inFiles, checkIfExists:true)
    .set {fastqFiles}

/*
bismark uses many resources because it calls other programs
under the hood.
(1) for each fastq file (pair), it calls 2 (directional) or 4
(nondirectional) threads to align to genome.
(2) it also calls gzip, samtools to process fastq and output file.
(3) bismark itself processing in perl.
As you can see, one sample may use 4 or 6 cores, depending on whether
the input is directional or not, assuming (2) and (3) just use 2
cores.
Also bismark has options --parallel (aka --multicore) and -p to add
parallelization: the former creates more independent bismark run while
the latter change the number of threads used by each alignment
process. So increasing the latter by 1 will increase the total number
of cores by 2 or 4, depending on directional input or not.
*/

process bismark {
    echo true
    cpus 8
    memory "31 GB"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(reads) from fastqFiles
    path indexFolder from params.indexFolder

    output:
    path "*.bam" into ch_bam
    path "*_report.txt" into ch_bismark_report

    script:
    // calculate the threads used by bowtie2 aligner
    cpusForAligner=(task.cpus as int) - 2
    cpusPerAligner=(cpusForAligner/(params.directional ? 2 : 4)) as int
    alignerThreads=cpusPerAligner > 1 ? "-p $cpusPerAligner" : ""
    nonDirectional=params.directional ? "" : "--non_directional"
    fqFiles="-1 ${reads[0]} -2 ${reads[1]}"
    """
    bismark $nonDirectional --genome_folder $indexFolder \
        $alignerThreads $fqFiles -B $sampleId
    """
}

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Results are in --> $params.outdir\n" : "Oops .. something went wrong" )
}
