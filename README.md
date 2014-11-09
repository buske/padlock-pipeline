padlock-pipeline
================

Pipelines for analyzing padlock-captured DNA sequencing data



### 1. Initial configuration:


#### Install necessary scripts (bowtie, bismark, trimmomatic):

    make install

#### Reference preparation:

* Download fasta files for reference genome into subdirectory of `sequences` folder (or symlink them), e.g.:

    mkdir sequences/hg19
    cd sequences/hg19
    ln -s /path/to/hg19.fa

* Prepare bismark genomic indexes:

    lib/bismark_v0.12.3/bismark_genome_preparation sequences/phiX
    lib/bismark_v0.12.3/bismark_genome_preparation sequences/hg19


### 2. Run pipeline on a sample directory

Sample directory should contain paired-end files of the form: `<SAMPLE>_L001_R1_001.fastq.gz`, `<SAMPLE>_L001_R2_001.fastq.gz`)

Pipeline can then be run with SGE-compatible command:

    src/bis_seq/run_sample.sh hg19 /path/to/sample/
