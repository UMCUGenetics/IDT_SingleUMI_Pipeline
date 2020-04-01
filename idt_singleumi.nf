#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf' params(params)
include MergeFastQs from './NextflowModules/bash/MergeFastQs.nf' params(params)
include DemuxFastqs from './NextflowModules/fgbio/1.1.0/DemuxFastqs.nf' params(params)
include MakeUmiBam from './NextflowModules/python/MakeUmiBam.nf' params(params)
include SortBam from './NextflowModules/fgbio/1.1.0/SortBam.nf' params(params)
include CallMolecularConsensusReads from './NextflowModules/fgbio/1.1.0/CallMolecularConsensusReads.nf' params(params)
include FilterConsensusReads from './NextflowModules/fgbio/1.1.0/FilterConsensusReads.nf' params(params)
include CountUmiFamilies from './NextflowModules/python/CountUmiFamilies.nf' params(params)
include SamToFastq from './NextflowModules/GATK/4.1.3.0/SamToFastq.nf' params(params)

workflow{
  main:
    def input_fastqs

    /*  Check if all necessary input parameters are present */
    if (!params.fastq_path){
      exit 1, "Please provide a fastq_path!"
    }

    if (!params.out_dir){
      exit 1, "No 'out_dir' parameter found in config file!"
    }

    // Gather input FastQ's
    if (params.fastq_path){
      input_fastqs = extractAllFastqFromDir(params.fastq_path)
    }

    if (params.singleEnd){
      R1_channel = input_fastqs
        .map{
          sample_id, rg_id, machine ,run_nr,fastqs ->
          def flowcell = rg_id.split('_')[1]
          // fastqs[0].split("/")[-1].split
          [sample_id, flowcell, fastqs[0], fastqs[1], fastqs[2]]
        }
        .groupTuple()
        .multiMap{
            sample_id, flowcell, r1, i1, i2 ->
            R1_channel: [sample_id,flowcell[0], r1.toSorted()]
            I1_channel: [sample_id,flowcell[0], i1.toSorted()]
            I2_channel: [sample_id,flowcell[0], i2.toSorted()]
        }
        .set{fastq_channels}
      // MergeFastQs(R1_channel).view()
      MergeFastQs(
          fastq_channels.R1_channel
            .mix(fastq_channels.I1_channel)
            .mix(fastq_channels.I2_channel)
      )
    }else{
      input_fastqs
        .map{
          sample_id, rg_id, machine, run_nr, fastqs ->
          def flowcell = rg_id.split('_')[1]
          [sample_id, flowcell, fastqs[0],fastqs[1], fastqs[2],fastqs[3]]
        }
        .groupTuple()
        .multiMap{
          sample_id, flowcell, r1, r2, i1, i2 ->
          R1_channel: [sample_id,flowcell[0],r1.toSorted()]
          R2_channel: [sample_id,flowcell[0],r2.toSorted()]
          I1_channel: [sample_id,flowcell[0],i1.toSorted()]
          I2_channel: [sample_id,flowcell[0],i2.toSorted()]
        }
        .set{fastq_channels}
      MergeFastQs(
        fastq_channels.R1_channel
          .mix(fastq_channels.R2_channel)
          .mix(fastq_channels.I1_channel)
          .mix(fastq_channels.I2_channel)
      )
    }
    to_demux = MergeFastQs.out
      .groupTuple()
      .map{
        sample_id, flowcell, tags, fastqs ->
        [ params.sample_sheet, flowcell, fastqs.sort{it.name}, params.read_structures ]
      }
    to_annotate = DemuxFastqs(to_demux)
      .map{metrics,flowcell, fastqs ->
         fastqs.findAll{ ! (it =~ /.*unmatched.*/)}
       }
       .flatten()
       .map{ fastq ->
         def name = fastq.baseName.split('_')[0]
         def (flowcell, lane, machine, run_nr) = flowcellLaneFromFastq(fastq)
         [name, flowcell, machine, run_nr, fastq]
       }
       .groupTuple()
       .map{name, flowcells, machines, run_nrs, fastqs -> [name, flowcells[0], machines[0], run_nrs[0], fastqs]}

    MakeUmiBam(to_annotate)
    SortBam(MakeUmiBam.out)

    to_fastq = CallMolecularConsensusReads(SortBam.out)
    if (params.filter){
      to_fastq = FilterConsensusReads(CallMolecularConsensusReads.out)
    }
    SamToFastq(to_fastq)
    CountUmiFamilies(to_fastq)
  emit:
    SamToFastq.out



}
