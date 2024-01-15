#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.input_file = ""
params.output_dir = ""


params.forward_primer                = 'AGATGAACCCACTGTGGACC'
params.reverse_primer                = 'GCTGTGGTGGAGGCTGAG'
params.forward_primer_num_mismatches = 4
params.reverse_primer_num_mismatches = 3

params.IgG_motif = 'CTCAGCCTCCACCACAGC'
params.IgG_motif_num_mismatches = 3

params.pre_CDR3_motifs = ['GA[AG]GA[TC][ATG][GC][ATGC]GC[GCTA][ATC][CGT][ATCG].*',
                          '[ATG]C[GA]GCC[AG][CT][AG][TC]A[CT].*'
                         ]

params.post_CDR3_motifs = ['[ATG][GC][ATGC]GC.*?TGGGG[GC]C[GCA][AG]',
                           '[ATG][GC][ATGC]GC.*?TGGGCCAA',
                           '[ATG][GC][ATGC]GC.*?CCGGGCCAA',
                           '[ATG][GC][ATGC]GC.*?TGCGGCCGA',
                           '[ATG][GC][ATGC]GC.*?TGGGGTC[AG]G',
                           '[ATG][GC][ATGC]GC.*?TGTGGCCAG',
                           '[ATG][GC][ATGC]GC.*?TGGGGCTCA',
                           '[ATG][GC][ATGC]GC.*?AGGGGCCAA',
                           '[ATG][GC][ATGC]GC.*?TGGAGCCAG'
                          ] 
                         

params.CDR3_lengths_offset = 27
params.productive_sequences_offset = 9
params.translation_frame = [1]


include { RELABEL_SEQUENCES as simplify_read_labels } from './processes'
include { RELABEL_SEQUENCES as label_forward_primer_reads } from './processes'
include { RELABEL_SEQUENCES as label_reverse_primer_reads } from './processes'
include { RELABEL_SEQUENCES as label_reverse_complimented_reverse_primer_reads } from './processes'

include { RELABEL_SEQUENCES as label_IgG_motif_reads } from './processes'

include { RELABEL_SEQUENCES as clean_up_extracted_pre_CDR3_motif_headers } from './processes'
include { RELABEL_SEQUENCES as clean_up_extracted_post_CDR3_motif_headers } from './processes'

include { RELABEL_SEQUENCES as label_pre_CDR3_sequences } from './processes'
include { RELABEL_SEQUENCES as label_post_CDR3_sequences } from './processes'

include { FILTER_SEQUENCES as get_forward_primer_reads } from './processes'
include { FILTER_SEQUENCES as get_reverse_primer_reads } from './processes'

include { FILTER_SEQUENCES as get_IgG_motif_reads } from './processes'

include { FILTER_SEQUENCES as get_productive_sequences } from './processes'

include { REVERSE_COMPLEMENT_SEQUENCES as reverse_compliment_reverse_primer_reads } from './processes'

include { CONCAT_SEQUENCES as merge_forward_reverse_primer_reads } from './processes'

include { REMOVE_DUPLICATE_SEQUENCES as remove_duplicate_merged_forward_reverse_reads } from './processes'
include { REMOVE_DUPLICATE_SEQUENCES as remove_duplicate_pre_CDR3_sequences } from './processes'
include { REMOVE_DUPLICATE_SEQUENCES as remove_duplicate_post_CDR3_sequences } from './processes'

include { LOCATE_REGEX_MATCHES as get_pre_CDR3_locations } from './processes'
include { EXTRACT_MATCHES as extract_pre_CDR3_sequences } from './processes'

include { LOCATE_REGEX_MATCHES as get_post_CDR3_locations } from './processes'
include { EXTRACT_MATCHES as extract_post_CDR3_sequences } from './processes'


include { GET_MATCH_LENGTHS as pre_and_post_CDR3_match_lengths } from './processes'
include { GET_MATCH_LENGTHS as productive_sequence_lengths } from './processes'

include { TRANSLATE_TO_AA_SEQUENCE as translate_pre_and_post_CDR3_sequences } from './processes'

                         
workflow {
    def input_file_basename = params.input_file.split('/')[-1].split("\\.")[0]
    
    
    relabeled_sequences = simplify_read_labels(params.input_file, '(.*?)\s.*', '$1', 'dont_save.gz')

    reads_with_forward_primer = get_forward_primer_reads(relabeled_sequences, params.forward_primer, params.forward_primer_num_mismatches, false, 'dont_save.gz')
    reads_with_reverse_primer = get_reverse_primer_reads(relabeled_sequences, params.reverse_primer, params.reverse_primer_num_mismatches, false, 'dont_save.gz')
     
    labeled_reads_with_forward_primer = label_forward_primer_reads(reads_with_forward_primer, '(.*)', '$1\thas_forward_primer', "${input_file_basename}-reads_with_forward_primer.fastq.gz")
    labeled_reads_with_reverse_primer = label_reverse_primer_reads(reads_with_reverse_primer, '(.*)', '$1\thas_reverse_primer', "${input_file_basename}-reads_with_reverse_primer.fastq.gz")
    
    reverse_complimented_reverse_primer_reads              = reverse_compliment_reverse_primer_reads(reads_with_reverse_primer)
    labeled_reverse_complimented_reads_with_reverse_primer = label_reverse_complimented_reverse_primer_reads(reverse_complimented_reverse_primer_reads, '(.*)', '$1\thas_reverse_primer:reverse_complimented', 'dont_save.gz')
    
    reads_with_IgG_motif = get_IgG_motif_reads(labeled_reads_with_forward_primer, params.IgG_motif, params.IgG_motif_num_mismatches, false, 'dont_save.gz')
    
    labeled_reads_with_IgG_motif = label_IgG_motif_reads(reads_with_IgG_motif, '(.*)', '$1\thas_IgG_motif', 'dont_save.gz')
   
  
    merged_forward_reverse_primer_reads = merge_forward_reverse_primer_reads( labeled_reads_with_IgG_motif.collect() + labeled_reverse_complimented_reads_with_reverse_primer.collect() )
    merged_forward_reverse_primer_reads_duplicates_removed = remove_duplicate_merged_forward_reverse_reads(merged_forward_reverse_primer_reads, false, "${input_file_basename}-reads_with_an_IgG_motif_or_reverse_complimented_reverse_primer.fastq.gz")
    
    pre_CDR3_locations                           = get_pre_CDR3_locations(merged_forward_reverse_primer_reads_duplicates_removed, params.pre_CDR3_motifs)
    pre_CDR3_sequences                           = extract_pre_CDR3_sequences(merged_forward_reverse_primer_reads, pre_CDR3_locations, 'dont_save.gz')
    pre_CDR3_sequences_with_cleaned_header       = clean_up_extracted_pre_CDR3_motif_headers(pre_CDR3_sequences, '(.*?)_.*', '$1', 'dont_save.gz')
    
    pre_CDR3_sequnces_duplicates_removed         = remove_duplicate_pre_CDR3_sequences( pre_CDR3_sequences_with_cleaned_header, false, 'dont_save.gz' )
    labeled_pre_CDR3_sequnces_duplicates_removed = label_pre_CDR3_sequences(pre_CDR3_sequnces_duplicates_removed, '(.*?):.*', '$1\textracted_pre_CDR3_motif', "${input_file_basename}-pre_CDR3_motif_sequences.fastq.gz")


    post_CDR3_locations                           = get_post_CDR3_locations(pre_CDR3_sequnces_duplicates_removed, params.post_CDR3_motifs)
    post_CDR3_sequences                           = extract_post_CDR3_sequences(pre_CDR3_sequnces_duplicates_removed, post_CDR3_locations, 'dont_save.gz')
    post_CDR3_sequences_with_cleaned_header       = clean_up_extracted_post_CDR3_motif_headers(post_CDR3_sequences, '(.*?)_.*', '$1', 'dont_save.gz')
    post_CDR3_sequences_duplicates_removed        = remove_duplicate_post_CDR3_sequences( post_CDR3_sequences_with_cleaned_header, false, 'dont_save.gz' )
    labeled_post_CDR3_sequnces_duplicates_removed = label_post_CDR3_sequences(post_CDR3_sequences_duplicates_removed, '(.*?):.*', '$1\thas_pre_CDR3_motif\textracted_pre_and_post_CDR3_motif', "${input_file_basename}-pre_and_post_CDR3_motifs_sequences.fastq.gz")
    
    pre_and_post_CDR3_match_lengths( post_CDR3_sequences_duplicates_removed, params.CDR3_lengths_offset, "${input_file_basename}-CDR3_lengths.tsv" )
    
    
    translated_sequences = translate_pre_and_post_CDR3_sequences(post_CDR3_sequences_duplicates_removed, params.translation_frame)
    productive_sequences = get_productive_sequences(translated_sequences, '*', 0, true, "${input_file_basename}-productive_sequences.faa.gz")
    productive_sequence_lengths(productive_sequences, params.productive_sequences_offset, "${input_file_basename}-productive_sequences_lengths.tsv")
}




