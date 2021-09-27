/*
 * This pipeline does the following things:
 * 1. call methylation values using Methyldackel
 */

// nextflow.enable.dsl=2

// default output dir
params.outdir = "bedgraph"
params.taskCpus = 1

log.info """\
         Call Methylation from Bam files - N F   P I P E L I N E
         ===================================
         genome : ${params.genome}
         bamFiles : ${params.bamfiles}
         outdir : ${params.outdir}
         awsQueue: '${params.awsQueue}'
         awsContainer: ${params.container}
         """
         .stripIndent()

if(workflow.profile == 'batch') {
    if(!params.awsQueue) {
        exit 1, "For profile 'batch', parameter 'awsQueue' is required"
    }
    if(!params.container) {
        exit 1, "For profile 'batch', 'process.container' need be set"
    }
}

Channel
    .fromPath( params.bamfiles, checkIfExists:true)
    .set {bamFiles}

/*
println "This is a string"
println bamFiles
bamFiles.view()
// exit 2, "aborted". // uncomment this line will lead to skipping of
// previous line's output
*/

// index genome file
process index {
    input:
    path genome from params.genome // generates value channel

    output:
    tuple path(genome), path(faIndex) into fa_index // also value channel

    script:
    faIndex=genome + '.fai'
    """
    samtools faidx $genome
    """
}

process call_meth {
    
    publishDir params.outdir, mode: 'copy'
    cpus params.taskCpus
    echo true

    input:
    path f from bamFiles
    tuple path(genome), path(faIndex) from fa_index

    output:
    path "*.bedGraph.gz", emit: outfile

    script:
    """
    #echo \$PATH
    echo Processing $f
    out=\$(basename $f)
    out=\${out%_R1_val_*}
    MethylDackel extract -@ ${task.cpus} -o \$out --CHH --CHG \
        $genome $f && gzip *.bedGraph
    """
}

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Results are in --> $params.outdir\n" : "Oops .. something went wrong" )
}
