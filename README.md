# mirkwater
Workflow that uses SARS-CoV-2 'constellations' from [https://github.com/cov-lineages/constellations](https://github.com/cov-lineages/constellations) to identify the best match lineages of consensus sequence clusters obtained from PacBio HiFi Amplicon sequencing of wastewater metagenomic samples.

Requirements:
* [Nextflow](https://www.nextflow.io/)
* Container Engine - (Docker/Singularity)

Usage:
`./mirkwater.nf --pacbio <path to dir with ccs bam files>`


*Name note: Mirkwood is a dark fictional forest in novels written by Sir Walter Scott, William Morris, and J. R. R. Tolkien. Since this workflow looks SARS-CoV-2 genomes in wastewater the name seemed fitting.*
