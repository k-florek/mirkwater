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

//get snpEff annotation database
Channel
  .fromPath(params.snpEff_database, checkIfExists: true, type: 'file')
  .set { snpeff_database_input }


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
  cat ${clusters} > all.fa
  mafft --6merpair --adjustdirection --addfragments all.fa ${reference} > cluster_alignment.fasta
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
  snp-sites -cb -v -o cluster_variants.vcf ${alignment}
  """
}

//Step5.a: annotate variants
process prep_variant_calls {
  publishDir "${params.outdir}/cluster_variants", mode: 'copy'

  input:
  file(vcf) from cluster_variants

  output:
  file("*sample*.vcf") into prep_variant_calls

  script:
  """
  bcftools annotate --rename-chrs <(printf '1\tNC_045512.2') ${vcf} > cluster_variants_chrom.vcf
  for sample in \$(bcftools query -l cluster_variants_chrom.vcf); do
    bcftools view -c1 -Ov -s \$sample -o \$sample.vcf cluster_variants_chrom.vcf
  done
  mv cluster_variants_chrom.vcf ${vcf}
  """
}

//combine vcfs with snpEff database
prep_variant_calls
  .flatten()
  .combine(snpeff_database_input)
  .set { annotation_input }

//Step6.a: annotate variants
process annotate_variants {
  publishDir "${params.outdir}/annotated_variants", mode: 'copy'

  input:
  tuple file(vcf), file(compressed_database) from annotation_input

  output:
  file(vcf) into annotated_variants

  script:
  """
  tar -xzf ${compressed_database}
  mv ${vcf} input.vcf
  snpEff ann sc2 input.vcf > ${vcf}
  """
}

//Step7.a find lineage hits
process find_hits {
  publishDir "${params.outdir}/lineage_matches", mode: 'copy'

  input:
  file(annotation) from annotated_variants

  output:
  file("*.results.csv") into lineage_hits

  script:
  """
  python /usr/src/app/partial_matching.py /usr/src/app/definitions ${annotation}
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

//Final Step
process summary {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(hits) from lineage_hits.collect()

  output:
  file("summary.csv") into results

  script:
  """
  echo 'PBAA Cluster,Best Lineage Match,Matching Mutations,All Identified Mutations' > summary.csv
  cat *.results.csv >> summary.csv
  """
}
