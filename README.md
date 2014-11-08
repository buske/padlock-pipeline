padlock-pipeline
================

Pipelines for analyzing padlock-captured DNA sequencing data



1) Initial configuration:


a) Install necessary scripts (bowtie, bismark, trimmomatic):

  `$ make install`

b) Reference preparation:

- Download fasta files for reference genome into subdirectory of `sequences` folder (or symlink them)

  e.g. `$ mkdir sequences/hg19; cd sequences/hg19; ln -s /path/to/hg19.fa`

- Prepare bismark genomic indexes:

  `$ lib/bismark_v0.12.3/bismark_genome_preparation sequences/phiX`
  `$ lib/bismark_v0.12.3/bismark_genome_preparation sequences/hg19`


2) Running a sample:

  `$ src/bis_seq/run_sample.sh hg19 /path/to/sample/`
