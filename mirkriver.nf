#!/usr/bin/env nextflow

//Description: Workflow for examining SARS-CoV-2 genomes from wastewater
//Author: Kelsey Florek
//email: kelsey.florek@slh.wisc.edu

//setup channel to read in sequencing reads
if (params.pacbio) {
  Channel
    .fromPath("${params.pacbio}/*.bam")
    .set { pacbio_reads }
}
else {
  println("No input...exiting.")
  System.exit(1)
}

//get reference geneome
Channel
  .fromPath(params.reference_genome, checkIfExists: true, type: 'file')
  .into { reference_genome; reference_genome_cluster }

//get reference genome annotation
Channel
  .fromPath(params.reference_gff, checkIfExists: true, type: 'file')
  .set { reference_genome_gff }


//combine reads with reference
pacbio_reads
  .combine(reference_genome)
  .combine(reference_genome_gff)
  .set { pacbio_reads_reference }

//Step1: extract fastq from bamfile
process preProcessPacBioBAM {
  tag "${reads.baseName}"

  input:
  set file(reads), file(reference), file(reference_gff) from pacbio_reads_reference

  when:
  params.pacbio

  output:
  tuple val("${reads.baseName}"), file("*.fastq"), file(reference), file(reference_gff) into fastq_reads_mapping, fastq_reads_pbaa
  file("*fasta.fai") into reference_genome_index_pbaa
  file("*fastq.fai") into read_index

  script:
  """
    samtools fastq ${reads} > ${reads.baseName}.fastq
    samtools faidx ${reference}
    samtools fqidx ${reads.baseName}.fastq
  """
}

//Step2.a: get pbaa clusters
process pbaa {
  tag "$name"
  publishDir "${params.outdir}/pbaa_clustering", pattern: '*passed*', mode: 'copy', saveAs: {"${name}.clusters.passing.${file(it).getExtension()}"}
  publishDir "${params.outdir}/pbaa_clustering",pattern: '*failed*', mode: 'copy', saveAs: {"${name}.clusters.failed.${file(it).getExtension()}"}

  input:
  set val(name), file(reads), file(reference), file(reference_gff) from fastq_reads_pbaa
  file(index) from reference_genome_index_pbaa
  file(rd_index) from read_index

  output:
  file("*passed_cluster_sequences.fasta") into passed_clusters
  file("*failed_cluster_sequences.fasta")

  script:
  """

  pbaa cluster ${reference} ${reads} ${name}.sam --filter ${params.minreaddepth} --min-var-frequency ${params.minfreq}
  """
}

//Step3.a: align all clusters
process cluster_seq_alignment {
  publishDir "${params.outdir}/cluster_alignment", mode: 'copy'

  input:
  file(clusters) from passed_clusters.collect()
  file(reference) from  reference_genome_cluster

  output:
  file("cluster_alignment.fasta") into cluster_alignment

  script:
  """
  cat *.fasta > all.fa
  mafft all.fa > cluster_alignment.fasta
  """
}

//Step4.a: get variants
process cluster_variant_sites {
  publishDir "${params.outdir}/cluster_variants", mode: 'copy'

  input:
  file(alignment) from cluster_alignment

  output:
  file("cluster_variants.vcf") into cluster_variants

  script:
  """
  snp-sites -v -o cluster_variants.vcf ${alignment}
  """
}

//Step2.b: map fastq to reference using minimap2
process mapping_reads {
  tag "$name"
  publishDir "${params.outdir}/mapped", mode: 'copy'

  input:
  set val(name), file(reads), file(reference), file(reference_gff) from fastq_reads_mapping

  output:
  tuple name, file("*.sam"), file(reference), file(reference_gff) into mapped_reads

  script:
  """
  minimap2 -ax splice:hq ${reference} ${reads} > ${name}.sam
  """
}

//Step3.b: convert sam to sorted indexed bam
process ivar_variant_calling {
  tag "$name"
  publishDir "${params.outdir}/variant_calls", mode: 'copy', saveAs: {"${name}.${file(it).getExtension()}"}

  memory {2.GB * task.attempt}
  errorStrategy {'retry'}
  maxRetries 2

  input:
  set val(name), file(reads), file(reference), file(annotation) from mapped_reads

  output:
  file("*.tsv")

  script:
  """
    samtools view -b ${reads} | samtools sort > ${name}.bam
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference ${reference} ${name}.bam | ivar variants -p ${name} -q ${params.minquality} -t ${params.minfreq} -m ${params.minreaddepth} -r ${reference} -g ${annotation}
  """
}
