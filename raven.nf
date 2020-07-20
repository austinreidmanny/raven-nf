#!/usr/bin/env nextflow

/*
   ------------------------------------------------------------------------------------------------
   WELCOME TO THE raven-nf PIPELINE
   ------------------------------------------------------------------------------------------------
   Usage:

   nextflow run raven.nf --sra "SRR123456"

   optional parameters:
   --diamondDB "path/to/db.dmnd" [default = "/n/data1/hms/mbin/nibert/austin/diamond/nr.dmnd"]
   --nucl "dna" or "rna" [default = "rna"]
   --output "output/dir/" [default="current-directory/output/"]
   --threads # (number of CPUs) [default=4]
   --memory # (amount of RAM) [in GB, default=8]
   --tempdir "tmp/dir" [default="/tmp"]
   ------------------------------------------------------------------------------------------------
   Contact:

   Austin R. Manny
   Nibert Lab @ Harvard Medical School
   austinmanny@g.harvard.edu
   github.com/austinreidmanny
   ------------------------------------------------------------------------------------------------
   Complete background, project objectives, and information available @
   github.com/austinreidmanny/raven-nf/README.md
   ------------------------------------------------------------------------------------------------
*/

// ---------------------------------------------------------------------------------------------- //
// Welcome the user
// ---------------------------------------------------------------------------------------------- //

println "\n====================================================================================\n" +
        "Welcome to the 'raven-nf' pipeline! To run it, just type the following: \n\n" +
        "   nextflow run raven.nf --sra 'SRA00001,SRA00002,SRA00003' --diamondDB 'path/to/db.dmnd' --nucl 'dna' or 'rna' \n\n" +
        "For detailed information on this pipeline, refer to the 'README.md' file or visit \n" +
        "www.github.com/austinreidmanny/raven-nf/ \n" +
        "====================================================================================\n"

// ---------------------------------------------------------------------------------------------- //

// ---------------------------------------------------------------------------------------------- //
// Main code
// ---------------------------------------------------------------------------------------------- //

// Create a 'run_name' parameter for naming files;
//     if single SRA accession given, name it that; if multple SRAs: "SRA_first-SRA_last"
if (params.sra.split(",").size() > 1) {

    run_name = params.sra.split(",")[0] +
               "-" +
               params.sra.split(",")[-1]
} else {
    run_name = params.sra
}

process log_inputs {
    // Save all inputs used for this run; important because files are named according to 'run_name'
    // which only indicates the first and last sample

    publishDir "${params.output}/00_analysis_info/", mode: "copy"

    input:
    val samples from params.sra

    output:
    file "${run_name}.readme.txt"

    """
    echo "Pipeline began at: \$(date)" > \
         "${run_name}.readme.txt"

    echo "Input samples: $samples" >> \
         "${run_name}.readme.txt"

    echo "Reference database to map to: $params.diamondDB" >> \
         "${run_name}.readme.txt"
    """

}

process parse_sra_ids {
    // Process the SRA accessions provided by the user

    input:
    val sra_id from params.sra.split(",")

    output:
    val sra_id into sra_accessions

    """
    echo "Preparing run for SRA ID: $sra_id ..."
    """

}

process download_sra_files {
    // Take in each SRA accession number, download the files, and send them to the mapping process

    input:
    val sra_id from sra_accessions

    output:
    file "${sra_id}.fq.gz" into sra_fastqs

    """
    fasterq-dump \
      -o "${sra_id}.fq" \
      -O ./ \
      -b 100MB \
      -c 500MB \
      --mem "${task.memory.toGiga()}G" \
      --temp $params.tempdir \
      --threads ${task.cpus} \
      --progress \
      --split-spot \
      --skip-technical \
      --rowid-as-name \
      $sra_id

    # Compress the output
    gzip "${sra_id}.fq"
    """

}

process combine_reads {

    // Combine SRA FASTQs into one file

    input:
    file reads from sra_fastqs.collect()

    output:
    file "merged_reads.fq.gz" into merged_fastq, merged_fastq_for_taxonomy_analysis, reads_for_refinement

    """
    zcat $reads | gzip > "merged_reads.fq.gz"
    """

}

process de_novo_assembly {

    // Assemble reads into fuller-length transcripts

    publishDir path: "${params.output}/01_de_novo_assembly",
               pattern: "${run_name}.transcripts.fasta",
               mode: "copy"

    input:
    file reads from merged_fastq

    output:
    tuple file("${run_name}.transcripts.fasta"), file(reads) into contigs_and_reads
    file "${run_name}.transcripts.fasta" into contigs_for_identifying_viruses

    """
    # Build contigs with rnaSPAdes [note: need to add feature to switch to metaSPAdes for DNA]
    rnaspades.py \
      -s $reads \
      -o unfiltered_assemblies/ \
      --memory "${task.memory.toGiga()}" \
      --threads $task.cpus \
      --tmp-dir $params.tempdir

    # Filter out (drop) any short contigs <300 nt
    seqtk seq -L 300 unfiltered_assemblies/transcripts.fasta > "${run_name}.transcripts.fasta"
    """

}

process coverage {

    // Map the reads to the contigs to determine per-contig coverage

    publishDir path: "${params.output}/02_coverage/",
               pattern: "${run_name}.contigs_coverage.txt",
               mode: "copy"

    input:
    tuple file(contigs), file(reads) from contigs_and_reads

    output:
    tuple file(contigs), file("${run_name}.contigs_coverage.txt") into contigs_with_coverage

    """
    # Index contigs for BWA
    bwa index -p "${run_name}_index" $contigs

    # Map reads to contigs with BWA-mem
    bwa mem -t $task.cpus "${run_name}_index" $reads | \
    samtools sort --threads $task.cpus -o "${run_name}.mapped.bam"

    # Calculate the mean-depth (i.e., coverage) per contig; keep each contig's name & coverage; throw away header; sort by coverage
    samtools coverage "${run_name}.mapped.bam" | \
    cut -f 1,7 > "${run_name}.contigs_coverage.txt"
    """
}

process classify_contigs {

    // Use DIAMOND to taxonomically classify each assembly

    publishDir path: "${params.output}/03_contigs_classification/",
               pattern: "${run_name}.contigs_classification.txt",
               mode: "copy"

    input:
    tuple file(contigs), file(coverage) from contigs_with_coverage

    output:
    tuple file("${run_name}.contigs_classification.txt"), file(coverage) into classified_contigs

    """
    # Run diamond
    diamond \
    blastx \
    --verbose \
    --more-sensitive \
    --db $params.diamondDB \
    --query $contigs \
    --out "${run_name}.contigs_classification.txt" \
    --outfmt 102 \
    --max-hsps 1 \
    --top 1 \
    --block-size $params.blocksize \
    --index-chunks 2 \
    --threads $task.cpus \
    --tmpdir $params.tempdir
    """

}

process taxonomy {

    // Save translated classification files containing the full taxonomic lineages

    publishDir path: "${params.output}/04_contigs_taxonomy/",
               pattern: "${run_name}.classification.taxonomy.txt",
               mode: "copy"

    publishDir path: "${params.output}/04_contigs_taxonomy/",
              pattern: "${run_name}.final_table.txt",
              mode: "copy"

    input:
    tuple file(classifications), file(coverage) from classified_contigs

    output:
    file "${run_name}.final_table.txt" into table_with_coverage_and_taxonomy

    """
    # Translate the DIAMOND results to full lineages
    diamond_to_taxonomy.py $classifications

    # Join the coverage values and the taxonomy results
    join \
        -j 1 \
        -t \$'\t' \
        --check-order \
        <(sort -k1,1 $coverage) \
        <(grep -v "^#" "${run_name}.contigs_classification.taxonomy.txt" | sort -k1,1) | \
    sort -rgk2,2 > \
    "${run_name}.contigs_coverage_taxonomy.txt"

    # Make a header for a final results table
    echo -e \
        "#Contig\t" \
        "#Coverage\t" \
        "#TaxonID\t" \
        "#e-value\t" \
        "#Domain\t" \
        "#Kingdom\t" \
        "#Phylum\t" \
        "#Class\t" \
        "#Order\t" \
        "#Family\t" \
        "#Genus_species" \
    > "${run_name}.final_table.txt"

    # Add the data to the final with just the header
    cat "${run_name}.contigs_coverage_taxonomy.txt" >> \
        "${run_name}.final_table.txt"
    """

}

process classify_reads {

    // Now that the contigs are assembled and classified, I would like to also do a metatranscriptomic
    // census of just the unassembled reads

    publishDir path: "${params.output}/05_unassembled_reads_taxonomy/",
               pattern: "${run_name}.taxonomy-of-reads.output.txt",
               mode: "copy"

    publishDir path: "${params.output}/05_unassembled_reads_taxonomy/",
               pattern: "${run_name}.taxonomy-of-reads.report.txt",
               mode: "copy"

    input:
    file reads from merged_fastq_for_taxonomy_analysis

    output:
    file "${run_name}.taxonomy-of-reads.report.txt" into classified_reads
    file "${run_name}.taxonomy-of-reads.output.txt"

    """
    kraken2 \
    --db $params.krakenDB \
    --gzip-compressed --memory-mapping \
    --threads $task.cpus \
    --output "${run_name}.taxonomy-of-reads.output.txt" \
    --report "${run_name}.taxonomy-of-reads.report.txt" \
    $reads
    """

}

process visualize_reads {

    // Visualize the classification of the reads from the 'classify_reads' process

    publishDir path: "${params.output}/05_unassembled_reads_taxonomy/",
               mode: "copy"

    input:
    file classifications from classified_reads

    output:
    file "${run_name}.taxonomy-of-reads.visualization.html"

    """
    ImportTaxonomy.pl \
    -m 3 -t 5 \
    $classifications \
    -o "${run_name}.taxonomy-of-reads.visualization.html"
    """

}

process identify_viral_assemblies {

    // Identify the viruses for refinement

    publishDir path: "${params.output}/06_viruses/",
               pattern: "${run_name}.viruses.txt",
               mode: "copy"

    input:
    file table from table_with_coverage_and_taxonomy
    file contigs from contigs_for_identifying_viruses

    output:
    file "${run_name}.viruses.txt" into viruses_table
    file "${run_name}.viruses.fasta" into viral_assemblies

    """
    awk '\$5 == "Viruses" {print}' $table > "${run_name}.viruses.txt"
    seqtk subseq $contigs <(cut -f 1 "${run_name}.viruses.txt") $contigs > "${run_name}.viruses.fasta"
    """

}

process refine_viral_assemblies {

    /*
       ------------------------------------------------------------------------
       Map the raw reads to the denovo-assembled-contigs with BWA;
       with the reads mapped to their denovo-assemblies, refine the assembly by
       finding any mismatches where rnaSPAdes called something different than
       what is shown by a pileup of the reads themselves
       ------------------------------------------------------------------------
    */

    // Save the refined viral contigs
    publishDir path: "${params.output}/06_viruses",
               pattern: "${run_name}.viruses.fasta.gz",
               mode: "copy"

    // Save the variants-called bcf file
    publishDir path: "${params.output}/06_viruses/mapping_files/",
               pattern: "${run_name}.variants_called.bcf",
               mode: "copy"

    // Save the BAM file with the name of the TVV species
    publishDir path: "${params.output}/06_viruses/mapping_files/",
               pattern: "${run_name}.reads_mapped_to_contigs.sorted.bam",
               mode: "copy"

    // Save the mapping-statistics file with the name of the TVV species
    publishDir path: "${params.output}/06_viruses/mapping_files/",
              pattern: "${run_name}.reads_mapped_to_contigs.sorted.stats",
              mode: "copy"

    input:
    file viral_assemblies
    file reads from reads_for_refinement

    output:
    file "${run_name}.refined_contigs.fasta.gz"
    val "true" into pipeline_complete

    """
    # Run the refinement script
    bash refine_contigs.sh \
    -s "${run_name}" -r $reads -c $viral_assemblies -o "./" -t $task.cpus
    
    # Rename the refined viruses file with a more descriptive name
    mv "${run_name}.refined_contigs.fasta.gz" "${run_name}.viruses.fasta.gz"
    """
}

process print_results {

    // Create a short summary with the number of virus assemblies, number of unique virus families,
    // and info on the longest viral contig

    input:
    file viruses_table
    val pipeline_complete

    output:
    stdout final_results

    """
    # Count number of virus contigs
    echo "Number of viral sequence assemblies in ${run_name}: \$(wc -l $viruses_table)"

    # Print mapped reads per virus family
    echo "Mapped reads per each identified virus family:"
    awk '{a[\$10] += \$2} END{for (i in a) print a[i], i}' < $viruses_table

    # Print the longest virus assembly constructed
    echo "Longest viral sequence assembled:"
    head -n 1 $viruses_table
    """

}

final_results.view{ it }

// ---------------------------------------------------------------------------------------------- //
