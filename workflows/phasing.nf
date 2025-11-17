include { beagle } from '../modules/phasing'

workflow phasing_beagle {
    take:
    vcf_unphased
    refcsv

    main:
    beagle(vcf_unphased, refcsv)
    emit:
    vcf = beagle.out.vcf
    vcf_single = beagle.out.vcf_single
}
