include { concatenate_phased } from "../modules/postprocessing"
include { fasta_index } from "../modules/preprocessing"

workflow postprocessing_ref{
    take:
    phased_vcf
    output_prefix
    fasta
    main:
    fasta_index(fasta)
    concatenate_phased(phased_vcf, output_prefix, fasta_index.out.index)
    emit:
    phased_variants = concatenate_phased.out.phased_variants
}