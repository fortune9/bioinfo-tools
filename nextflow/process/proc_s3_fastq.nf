params.r1tag="_R1_"
params.r2tag="_R2_"

process pair_fastq {
    debug params.verbose
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val prjId
    val s3path

    output:
    path 'fq_list.tsv', emit: ch_paired_fq

    script:
    """
    pair_fastq_files --r1-tag '${params.r1tag}' \
        --r2-tag '${params.r2tag}' $prjId $s3path \
        >fq_list.tsv
    """
}

// combine the fastq files into one, .e.g, from different
// runs
process combine_fastq {
    debug params.verbose
    label 'processMedium'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(outfile), path(fqs)

    output:
    path outfile, emit: ch_fastq_file

    script:
    """
    # if fqs is size 1, just need rename input
    for f in $fqs
    do
        mycat=cat
        if [[ "\$f" =~ \\.(gz|GZ)\$ ]]; then
            mycat=zcat
        fi
        if [[ $outfile =~ \\.(gz|GZ)\$ ]]; then
            \$mycat \$f | gzip -c >>$outfile
        else
            \$mycat \$f >>$outfile
        fi
    done
    """
}
