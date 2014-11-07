padlock-pipeline
================

Pipelines for analyzing padlock-captured DNA sequencing data



Initial configuration:

1) Download fasta files for reference genome into subdirectory of `sequences` folder (or symlink them)

  e.g. `$ mkdir sequences/hg19; cd sequences/hg19; wget ....fq.gz`

2) Prepare bismark genomic indexes:

  e.g. `$ lib/bismark_v0.12.3/bismark_genome_preparation refs/hg19`

3) Install necessary scripts:

  `$ make install`


Running a sample:

