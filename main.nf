// parameters

params.accession = "M21012"          // default if no NCBI accession number is provided
params.in = null               // Folder with sample FASTA files
params.outdir  = "./results"        // Results folder

process download_refseq {
    conda "bioconda::entrez-direct=24.0"

    tag "$accession"

    input:
        val accession

    output:
        path "reference.fasta"

    publishDir params.outdir,  mode: 'copy'

    script:
    """

    esearch -db nucleotide -query "${accession}" | efetch -format fasta > reference.fasta

    """
}

process combine_fastas {

    tag "combine_fastas"

    input:
    path fastas

    output:
    path "combined.fasta"

    publishDir params.outdir, mode: 'copy'

    script:
    """
    echo "Combining sample FASTA files"
    cat ${fastas.join(' ')} > combined.fasta
    """
}

process mafft_align {
    conda "bioconda::mafft=7.525"

    tag "mafft_alignment"

    input:
    path ref
    path combined

    output:
    path "aligned.fasta"

    publishDir params.outdir, mode: 'copy'
    script:
    """
    echo "Running MAFFT alignment"
    cat ${ref} ${combined} > alignment_input.fasta
    mafft --auto alignment_input.fasta > aligned.fasta
    """
}

process trimal_clean {
    conda "bioconda::trimal=1.5.0"

    tag "trimal_clean"

    input:
    path aligned

    output:
    path "aligned_trimmed.fasta"
    path "trimal_report.html"

    publishDir params.outdir, mode: 'copy'
    script:
    """
    echo "Cleaning alignment with TrimAl"
    trimal \\
        -in ${aligned} \\
        -out aligned_trimmed.fasta \\
        -htmlout trimal_report.html \\
        -automated1
    """
}



workflow {
    if ( !params.in ) {
        error "ERROR: You must provide a directory with FASTA files using --in"
    }

    Channel
        .fromPath("${params.in}/*.fasta")
        .ifEmpty { error "ERROR: No FASTA files found in ${params.in}" }
        .set { sample_fastas }

    ref_ch = download_refseq(params.accession)

    // Step 2: Combine sample FASTAs
    combined_ch = combine_fastas(sample_fastas)

    // Step 3: Align reference + samples
    aligned_ch = mafft_align(ref_ch, combined_ch)

    // Step 4: Clean alignment with TrimAl
    trimal_clean(aligned_ch)
}