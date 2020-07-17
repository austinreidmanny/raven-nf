# raven-nf
Virus discovery pipeline built upon the Nextflow pipeline framework

Like the bird, `raven-nf` is smart and good at finding things (like viruses). Unlike the bird, this pipeline has not been cited by Edgar Allen Poe.

## Contact
**Austin Reid Manny**

Nibert Lab, Harvard Medical School

`austinmanny@g.harvard.edu`

## Objective
`raven-nf` is a virus discovery pipeline for finding viruses in publicly available transcriptome data.

Extension of the original `raven` pipeline (writte in Bash, available at github.com/austinreidmanny/raven)

This pipeline simply takes in one or more SRA accession number (for NGS datasets in the NCBI Sequence Read Archive)
and returns a list of viruses and their sequences. All software dependencies are handled by the `conda` package manager.

The steps executed by the pipeline are are:
- Log inputs
- Download reads [`NCBI sra-tools`]
- Adapter trimming [`TrimGalore!`]
- De novo contig assembly [`rnaSPAdes`]
- Assembly refinement to map reads to newly assembled contigs [`bwa-mem`]
- Taxonomic classification of assemblies [`DIAMOND`]
- Translate taxonomy [custom scripts using `JGI-ISF taxonomy server`]
- Mapping reads to assemblies to determine coverage values [`bwa-mem`]
- Binning [custom code to parse taxonomy results > `seqtk`]
- Taxonomic assignment of unassembled reads [`kraken2`]
- Visualization of reads-based taxonomic distribution [`krona`]
- Print results to user

## Requirements
Download the latest version of `raven-nf` at https://github.com/austinreidmanny/raven-nf

Dependencies & requirements are
 * Operating system: Linux/macOS
 * Nextflow (https://www.nextflow.io)
 * Conda (https://docs.conda.io/en/latest/miniconda.html)
 * DIAMOND database
 * Kraken database
 * Krona set up (Conda will install it, but might not put it in your Bash path or initiate the taxonomy database for this)

## Usage

```
nextflow run raven.nf --sra "SRR123456"

   optional parameters:
   --nucl "dna" or "rna" [default = "rna"]
   --output "output/dir/" [default="results/"]
   --threads # (number of CPUs) [default=4]
   --memory # (amount of RAM) [in GB, default=8]
   --tempdir "tmp/dir" [default="/tmp/"]
   --diamondDB "path/to/db.dmnd" [default = "/n/data1/hms/mbin/nibert/austin/diamond/nr.dmnd"]
   --krakenDB "path/to/kraken_db" [default = "/n/data1/hms/mbin/nibert/austin/tools/kraken2/kraken_viruses"]
 ```

