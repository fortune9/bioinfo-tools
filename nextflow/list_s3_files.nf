
Channel
    .fromPath(params.s3path, checkIfExists:true)
    .set {file_ch}

file_ch.println()
//file_ch.view() // The same

/* section for checking channel property
println file_ch.getClass()
//println file_ch.collect().getClass()
println file_ch.toList().getClass() // the same as above
*/
