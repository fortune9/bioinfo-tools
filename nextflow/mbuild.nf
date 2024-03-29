/*
 * This pipeline builds methylation matrices using mbuilder
 * with bedgraphs as input
 */

// nextflow.enable.dsl=2

// customizable parameters
params.outdir = "mbuilder_output"
params.taskCpus = 1
params.taskMemory = "2 GB"
params.useRefFile = false // use chromosome name reference file?
params.cytoContexts = 'CpG'
params.bgFiles = null
params.fileList = null
params.pathPrefix = null
//params.mbuilderContainer = \$mbuilderContainer

log.info """\
         Build Methylation Matrix - N F   P I P E L I N E
         ===================================
         bgFiles      : '${params.bgFiles}'
         outdir       : ${params.outdir}
         cytoContexts : ${params.cytoContexts}
         useRefFile   : ${params.useRefFile}
         awsQueue     : ${params.awsQueue}
         awsContainer : ${params.container}

         Example use:
         nextflow run mbuild.nf --bgFiles 's3://path/*.bedgraph.gz' \
            --useRefFile -with-docker <mbuilder-docker-img:tag>
         
         Alternatively, one can put file paths in a text file, one
         path per line, and feed this file to the parameter
         --fileList, e.g.:

         nextflow run mbuild.nf --fileList 'files.txt' \
            --pathPrefix  '/folder/to/files' \
            --useRefFile -with-docker <mbuilder-docker-img:tag>

         ===================================
         """
         .stripIndent()

/*
if(!params.mbuilderContainer) {
    exit 1, "Please set environment variable 'mbuilderContainer'"
}
*/

if(workflow.profile == 'batch') {
    if(!params.awsQueue) {
        exit 1, "For profile 'batch', parameter 'awsQueue' is required"
    }
    if(!params.container) {
        exit 1, "For profile 'batch', 'process.container' need be set"
    }
}

if(params.bgFiles != null) {
    log.info "Reading files from command line"
    Channel
    //.fromPath(params.bgFiles, checkIfExists:true)
        .fromPath(params.bgFiles.split(/\s+/).flatten(), checkIfExists:true)
        .set {bedgraphO}
} else if(params.fileList != null) {
    log.info "Reading files from a file ${params.fileList}"
    paths=file(params.fileList).readLines()
    if(params.pathPrefix != null) {
        paths=paths.collect() { "${params.prefix}/$it" }
    }
    // create channel now
    Channel
        .fromPath(paths, checkIfExists:true)
        .set {bedgraphO}
        
} else {
    exit 1, "Both --fileList and --bgFiles parameters are empty"
}

process gzip_bg {
    input:
    path f from bedgraphO

    output:
    // use includeInputs to include input files
    path "*.gz" includeInputs true into bedgraph

    script:
    if(! "$f".endsWith(".gz")) // "$f" is to convert an object to a string
        """
        gzip -f $f
        """
    else
        """
        echo $f is already gzipped
        """
}

process build_chr_ref {
    cpus 2
    memory '2 GB'

    input:
    path f from bedgraph

    output:
    path f into bedgraph_ch
    path 'refFile.txt' into ref_file_ch

    script:
    if(params.useRefFile) {
        """
        zcat $f | grep -v '^track' | \
        cut -f 1 | uniq | sort | uniq >refFile.txt
        """
    } else {
        """
        echo "No chr ref file" >refFile.txt
        """
    }
}

// merge all refFile.txt into one
process merge_chr_ref {
    cpus 2
    memory "2 GB"
    echo true

    input:
    path 'chr_ref*' from ref_file_ch.toList()

    output:
    path 'chrs.ref.csv' into chr_ref_ch

    script:
    """
    #echo `ls chr_ref*`
    cat chr_ref* | sort | uniq >chrs.ref.csv
    """
}

process mbuilder {
//    container "${params.mbuilderContainer}"
    publishDir params.outdir, mode: 'copy'
    cpus params.taskCpus
    memory "${params.taskMemory}"
    echo true

    input:
    path chr_ref from chr_ref_ch
    path bedgraph from bedgraph_ch.collect()

    output:
    path "mbuilder_output/*.txt", emit: logs
    path "mbuilder_output/*.csv", emit: csv
    path "mbuilder_output/{meth_cov,perc,total_cov}_matrix_*.txt.gz", emit: matrices

    script:
    cContexts = params.cytoContexts? "-x ${params.cytoContexts}" : ""
    refOpt = params.useRefFile? "-u -R" : ""
    """
    # echo $bedgraph
    echo curdir: `pwd`
    mbuilder bg -i . -d 0 -p ${task.cpus} $cContexts $refOpt
    """
}

workflow.onComplete {
        log.info ( workflow.success ? 
            "\nDone! Results are in --> $params.outdir\n" : 
            "Oops .. something went wrong" )
}

