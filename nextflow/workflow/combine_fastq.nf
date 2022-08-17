#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def usage() {

log.info """\
    Combine fastq files with the sample sample id
    into one file, .e.g, those from different runs

    Options:

    --s3path <path> the folder where fastq files are located

    --prj-id <str> the project id, which will be matched against
        the fastq filenames to find the project files. E.g., zr1234
        
    """
    .stripIndent()
}

params.s3path = ''
params.prjId = ''
params.outdir = 'combined_fqs'
params.r1tag = "_R1_"
params.r2tag = "_R2_"

if(!params.s3path || !params.prjId) {
    usage()
    exit 1, "Please provide required parameters"
}

include {pair_fastq; combine_fastq} from \
    "../process/proc_s3_fastq.nf" params (
        r1tag: params.r1tag,
        r2tag: params.r2tag,
        verbose: params.verbose,
        outdir: "${params.outdir}"
    )

workflow {
    pair_fastq(params.prjId, params.s3path)
    ch_paired_fq=pair_fastq.out.ch_paired_fq
    ch_paired_fq
        .splitCsv(header: true, sep: "\t")
        .map{ row -> [ row.id, row.read_1, row.read_2] }
        .groupTuple(by: 0)
        .take(3)
        .multiMap{ row -> 
            R1: ["${row[0]}.R1.fastq.gz", row[1]]
            R2: ["${row[0]}.R2.fastq.gz", row[2]] }
        .set {ch_fqs_by_id}
    //ch_fqs_by_id.R1.view()
    //ch_fqs_by_id.R2.view()
    combine_fastq(ch_fqs_by_id.R1)
}

/*
Channel
    .fromPath(params.infile)
    .splitCsv(header: false, skip:1, sep: "\t")
    .map{ row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]), file(row[4]) ) }
    .set {fqPairs}

fqPairs
    .view() { it[0] }

return

process merge_fq {
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(id), file(read11), file(read12), file(read21), file(read22) from fqPairs

    output:
    path "*.fastq.gz" optional true

    script:
    """
    echo "hello" >tmp.txt
    #merge_fastq.sh $id $read11 $read12 $read21 $read22
    """
}

*/

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Results are in --> $params.outdir\n" : "Oops .. something went wrong" )
}

