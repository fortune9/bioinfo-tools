#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def usage() {

log.info """\
    Combine fastq files with the sample sample id
    into one file, .e.g, those from different runs

    Options:

    Provide the following 2 options to get paired fastq files
    from AWS s3 path:

    --s3path <path> the folder where fastq files are located

    --prj-id <str> the project id, which will be matched against
        the fastq filenames to find the project files. E.g., zr1234
    
    Alternatively, one can use the following option to directly 
    provide a tab-delimited file with the following 3 fields:
    id  read_1  read_2

    --fqPairFile <path> file containing the paths to read 1 and read 2
        fastq files.
    """
    .stripIndent()
}

params.verbose = false
params.s3path = ''
params.prjId = ''
params.outdir = 'combined_fqs'
params.r1tag = "_R1_"
params.r2tag = "_R2_"
params.fqPairFile = ''

if(!params.s3path || !params.prjId) {
    if(!params.fqPairFile) {
        usage()
        exit 1, "Please provide required parameters"
    }
}

include {pair_fastq; combine_fastq as combine_r1; combine_fastq as combine_r2} from \
    "../process/proc_s3_fastq.nf" params (
        r1tag: params.r1tag,
        r2tag: params.r2tag,
        verbose: params.verbose,
        outdir: "${params.outdir}"
    )

workflow {
    // if fqPairFile is provided, use this file directly
    if(params.fqPairFile) {
        ch_paired_fq=Channel.fromPath(params.fqPairFile)
    } else {
        if(workflow.profile != 'standard') {
            error 2, "'pair_fastq' can only be run in profile 'standard'"
        }
        pair_fastq(params.prjId, params.s3path)
        ch_paired_fq=pair_fastq.out.ch_paired_fq
    }
    ch_paired_fq
        .splitCsv(header: true, sep: "\t")
        .map{ row -> [ row.id, row.read_1, row.read_2] }
        .groupTuple(by: 0)
       // .take(2)
        .multiMap{ row -> 
            R1: ["${row[0]}.R1.fastq.gz", row[1]]
            R2: ["${row[0]}.R2.fastq.gz", row[2]] }
        .set {ch_fqs_by_id}
    //ch_fqs_by_id.R1.view()
    //ch_fqs_by_id.R2.view()
    combine_r1(ch_fqs_by_id.R1)
    combine_r2(ch_fqs_by_id.R2)
}

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Results are in --> $params.outdir\n" : "Oops .. something went wrong" )
}

