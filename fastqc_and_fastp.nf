#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.input_file  = ""
params.output_dir = ""

params.min_read_length = 200
params.max_read_length = 800

include { FASTQC as fastqc_pre_fastp } from './processes'
include { FASTQC as fastqc_post_fastp } from './processes'

include { FASTP as fastp } from './processes'




                         
workflow {
    def input_file_basename = params.input_file.split('/')[-1].split("\\.")[0]
    
    fastqc_pre_fastp(params.input_file)
    fastp(params.input_file, input_file_basename, params.min_read_length, params.max_read_length)
    fastqc_post_fastp(fastp.out.trimmed_file)

}




