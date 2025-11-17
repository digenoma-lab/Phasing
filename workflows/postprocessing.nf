include { concatenate_phased } from "../modules/postprocessing"

workflow postprocessing_ref{
    take:
    phased_vcf
    output_prefix
    main:
    concatenate_phased(phased_vcf, output_prefix)
    emit:
    phased_variants = concatenate_phased.out.phased_variants
}