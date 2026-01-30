// a process that downloads a FASTA file given an accession number

process download_refseq {
    conda "bioconda::entrez-direct=24.0"

    input:
        val accession

    output:
        path "reference.fasta" 

    script:

    """
    esearch -db nucleotide -query "$accession" \
    | efetch -format fasta > "reference.fasta" 

    """

}


process compile_seq {

   
    input:
        path  samples

    output:
        path "combined_seq.fasta"

    script:

    """
    cat ${samples} > combined_seq.fasta

    """


}

process align_seq {

    conda "bioconda::mafft=7.525"

    input:
        path combined_fasta
        

    output:
        path "sequences_aligned.fasta", emit: aligned

    script:

    """
    mafft --auto $combined_fasta > sequences_aligned.fasta

    """
}




// parameters defined here
params.accession = "M21012"
params.sample_seq = null



workflow {

 if (params.sample_seq == null) {
        println("no sample sequences provided")
        exit 1
    }

ch_input = channel.fromPath("${params.sample_seq}/*.fasta").collect()


download_refseq(params.accession)

compile_seq(ch_input)

align_seq()


}