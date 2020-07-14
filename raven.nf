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
    file "${sra_id}*fastq" into sra_fastqs

    """
    download_sra.sh -s $sra_id -t $params.tempdir -m $task.memory -n $task.cpus -o ./
    """

}

process combine_fastqs {
    
    // Combine SRA FASTQs into one file

    input:
    file reads from sra_fastqs.collect()
    
    output:
    file "merged_reads.fq.gz" into merged_fastq
    //file "merged_reads.fq" into merged_fastq
    
    """
    #cat $reads | gzip > "merged_reads.fq.gz"
    cat $reads > "merged_reads.fq.gz"
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
    file "${run_name}.transcripts.fasta"


    """
    # Build contigs with rnaSPAdes & drop any short contigs <300 nt

    rnaspades.py -s $reads -o unfiltered_assemblies/ --memory $params.memory --threads $params.threads
    seqtk seq -L 300 unfiltered_assemblies/transcripts.fasta > "${run_name}.transcripts.fasta"
    """

}

