/*
 * This pipeline does the following things:
 * 1. merge the cytosines from 
 */

// nextflow.enable.dsl=2

// default output dir
params.outdir = "mergeContext"
//params.taskCpus = 1

log.info """\
         Call Methylation from Bam files - N F   P I P E L I N E
         ===================================
         genome : ${params.genome}
         bgFiles : ${params.bgFiles}
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
    .fromPath( params.bgFiles, checkIfExists:true)
    .set {bgFiles}

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

process merge_context {
    
    publishDir params.outdir, mode: 'copy'
    cpus 1
    echo true

    input:
    path f from bgFiles
    tuple path(genome), path(faIndex) from fa_index

    output:
    path "*.bedGraph.gz", emit: outfile

    script:
    """
    #echo \$PATH
    echo Processing $f
    out=\$(basename $f)
    out=\${out%_R1_val_*}
    MethylDackel mergeContext $genome $f | \
        gzip -c >\${out}_CpG.bedGraph.gz
    """
}

workflow.onComplete {
        log.info ( 
            workflow.success ? 
                "\nDone! Results are in --> $params.outdir\n" : 
                "Oops .. something went wrong" 
                )
}
