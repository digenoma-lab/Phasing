include { concatenate_phased } from "../modules/postprocessing"

workflow postprocessing_ref{
    take:
    phased_vcf
    main:
    concatenate_phased(phased_vcf)
    emit:
    phased_variants = concatenate_phased.out.phased_variants
}