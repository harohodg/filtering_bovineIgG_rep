process FASTQC {
    module 'StdEnv/2020:fastqc/0.11.9'
    publishDir "${params.output_dir}/fastqc/${task.process}", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.html"
    publishDir "${params.output_dir}/fastqc/${task.process}", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.zip"
    
    input:
        path input_file
    output:
        path "*.html"
        path "*.zip"
    
    shell:
        '''
        fastqc --threads $(nproc) --outdir . !{input_file}
        '''   
}



process FASTP {
    module 'StdEnv/2020:fastp/0.23.1'
    publishDir "${params.output_dir}/fastp", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/fastp", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.json"
    publishDir "${params.output_dir}/fastp", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.html"
    
    input:
        path input_file
        val output_prefix
        val min_length
        val max_length
        
    output:
        path "${output_prefix}-trimmed.fastq.gz", emit: trimmed_file
        path "${output_prefix}-fastp_failed.fastq.gz", emit: failed_file
       
    shell:
        '''
            fastp --thread $(nproc) --dedup --disable_adapter_trimming --length_required !{min_length} --length_limit !{max_length} --json !{output_prefix}-fastp_report.json --html !{output_prefix}-fastp-report.html  --in1 !{input_file} --failed_out !{output_prefix}-fastp_failed.fastq.gz -o !{output_prefix}-trimmed.fastq.gz
        '''   
       
}


process RELABEL_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        val  header_match
        val  header_replacement
        val  output_file
    output:
        path "${output_file}"
    
    script:
        """
        seqkit replace --pattern '${header_match}' --replacement '${header_replacement}' --line-width 0 ${input_file} --out-file ${output_file}
        """
}


process FILTER_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.faa.gz"
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        
        val search_pattern
        val num_mismatches
        val invert_match
        val output_file
        
    output:
        path "${output_file}"
    
    script:
        def invert_flag = invert_match != false ? "--invert-match" : ''
        """
        seqkit grep --by-seq --only-positive-strand  --immediate-output --line-width 0 --pattern '${search_pattern}' --max-mismatch ${num_mismatches} ${invert_flag} ${input_file} --out-file ${output_file}
        """
}


process REVERSE_COMPLEMENT_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'

    output:
        path 'reverse_complimented_file.fastq.gz'
        
    script:
        """
        seqkit seq --reverse --complement --seq-type DNA --line-width 0 ${input_file} --out-file reverse_complimented_file.fastq.gz
        """
    
}




process FIND_COMMON_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    
    input:
        path input_file1, stageAs: 'input_file1.fastq.gz'
        path input_file2, stageAs: 'input_file2.fastq.gz'
        val by_sequence
        
    output:
        path 'common_sequences.fastq.gz'
    
    script:
        def by_seq_flag = by_sequence != false ? "--only-positive-strand --by-seq" : '' 
        """
        seqkit common ${input_file1} ${input_file2} ${by_seq_flag} --line-width 0  -o common_sequences.fastq.gz
        """
}

process CONCAT_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    input:
        val input_files
        
    output:
        path 'merged_sequences.fastq.gz'
        
    script:
        """
            seqkit seq  ${input_files.join(' ')} --out-file merged_sequences.fastq.gz
        """
}

process REMOVE_DUPLICATE_SEQUENCES {
    module 'StdEnv/2020:seqkit/2.3.1'
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        val by_sequence
        val output_file
        
    output:
        path "${output_file}"
    
    script:
        def by_seq_flag = by_sequence != false ? "--only-positive-strand --by-seq" : ''
        """
        seqkit rmdup ${by_seq_flag} ${input_file} --out-file ${output_file}
        """
}


process LOCATE_REGEX_MATCHES {
    module 'StdEnv/2020:seqkit/2.3.1'
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        val patterns
    output:
        path 'matches.bed'
    
    shell:
        '''
        seqkit locate --bed --only-positive-strand --use-regexp --pattern '!{patterns.join('|')}' !{input_file} > matches.bed
        '''
    
}

process EXTRACT_MATCHES {
    module 'StdEnv/2020:seqkit/2.3.1'
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.fastq.gz"
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        path bed_file
        val output_file
        
    output:
        path "${output_file}"
        
    script:
        """
        seqkit subseq --bed ${bed_file} ${input_file} --out-file ${output_file}
        """


}

process GET_MATCH_LENGTHS {
    module 'StdEnv/2020:bioawk/1.0'
    publishDir "${params.output_dir}/filtering", enabled: params.output_dir as boolean, mode: 'copy', overwrite: false, pattern: "*.tsv"
    
    input:
        path input_file
        val offset
        val output_file
    output:
        path "${output_file}"
    
    shell:
        '''
        bioawk -c fastx 'BEGIN {OFS="\t";print "read","match_length","match_length-!{offset}","sequence"} {print $name,length($seq), length($seq)-!{offset},$seq }' !{input_file} > !{output_file}
        '''
}

process TRANSLATE_TO_AA_SEQUENCE {
    module 'StdEnv/2020:seqkit/2.3.1'
    
    input:
        path input_file, stageAs: 'input_file.fastq.gz'
        val frames
    output:
        path 'translated_sequences.faa.gz'    
    
    script:
        """
        seqkit translate --append-frame --frame ${frames.join(',')} ${input_file} --out-file 'translated_sequences.faa.gz'
        """    
}




