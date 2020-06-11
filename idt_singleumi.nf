#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*  Check if all necessary input parameters are present */
if (!params.fastq_path){
  exit 1, "Please provide a fastq_path!"
}

if (!params.out_dir){
  exit 1, "No 'out_dir' parameter found in config file!"
}
if (!params.genome){
  exit 1, "No 'genome' parameter found in in <analysis_name>.config file or on the commandline (add -- in front of the parameter)."
}
if ( !params.genomes[params.genome] || !params.genomes[params.genome].fasta ){
  exit 1, "'genome' parameter ${params.genome} not found in list of genomes in resources.config!"
}

params.genome_fasta = params.genomes[params.genome].fasta



include './NextflowModules/Utils/fastq.nf' params(params)
include MergeFastqs from './NextflowModules/bash/MergeFastQs.nf' params(params)
include MakeUMIBam from './NextflowModules/python/MakeUmiBam.nf' params()
include SortBam from './NextflowModules/fgbio/1.1.0/SortBam.nf' params(optional: '--sort-order TemplateCoordinate',mem: "${params.sortbam.mem}")
include CallMolecularConsensusReads from './NextflowModules/fgbio/1.1.0/CallMolecularConsensusReads.nf' params(optional: '--error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=30 --min-reads 1 ',mem: "${params.callmolecularconsensusreads.mem}")
include FilterConsensusReads from './NextflowModules/fgbio/1.1.0/FilterConsensusReads.nf' params(optional: '--reverse-per-base-tags=true -M 1 -N 40 -e 0.1 -n 0.1',mem: "${params.filterconsensusreads.mem}")
include CountUMIFamilies from './NextflowModules/python/CountUmiFamilies.nf' params(optional: '')
include SamToFastq from './NextflowModules/GATK/4.1.3.0/SamToFastq.nf' params(optional: '',mem: "${params.samtofastq.mem}", genome_fasta: "${params.genome_fasta}")

workflow{
  main:
    def input_fastqs



    // Gather input FastQ's
    if (params.fastq_path){
      input_fastqs = extractAllFastqFromDir(params.fastq_path)
    }
    // input_fastqs.view()
    if (params.singleEnd){
      input_fastqs
        .map{
          sample_id, rg_id, machine ,run_nr,fastqs ->
          def flowcell = rg_id.split('_')[1]
          [sample_id, flowcell, fastqs[0], fastqs[1]]
        }
        .groupTuple()
        .multiMap{
          sample_id, flowcell, r1, umi ->
          R1_channel: [sample_id,flowcell[0], r1.toSorted()]
          UMI_channel: [sample_id,flowcell[0], umi.toSorted()]
        }
        .set{fastq_channels}
      MergeFastqs(
        fastq_channels.R1_channel
          .mix(fastq_channels.UMI_channel)
      )
    }else{
      input_fastqs
        .map{
          sample_id, rg_id, machine, run_nr, fastqs ->
          def flowcell = rg_id.split('_')[1]
          [sample_id, flowcell, fastqs[0],fastqs[1], fastqs[2]]
        }
        .groupTuple()
        .multiMap{
          sample_id, flowcell, r1, umi, r2 ->
          R1_channel: [sample_id,flowcell[0],r1.toSorted()]
          R2_channel: [sample_id,flowcell[0],r2.toSorted()]
          UMI_channel: [sample_id,flowcell[0],umi.toSorted()]
        }
        .set{fastq_channels}
      MergeFastqs(
          fastq_channels.R1_channel
            .mix(fastq_channels.R2_channel)
            .mix(fastq_channels.UMI_channel)
      )
    }
    to_annotate = MergeFastqs.out
      .groupTuple()
      .map{
        sample_id, flowcells, tags, fastqs ->
        def (flowcell, lane, machine, run_nr) = flowcellLaneFromFastq(fastqs[0])
        [sample_id, flowcells[0],machine, run_nr,fastqs.toSorted{a, b -> a.baseName <=> b.baseName} ]
      }
    MakeUMIBam(to_annotate)


    SortBam(MakeUMIBam.out).view()
  //
    to_fastq = CallMolecularConsensusReads(SortBam.out)
    if (params.filter){
      to_fastq = FilterConsensusReads(CallMolecularConsensusReads.out)
    }
    SamToFastq(to_fastq)
    CountUMIFamilies(to_fastq)
  emit:
    SamToFastq.out



}
