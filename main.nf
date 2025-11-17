include { preprocessing } from "./workflows/preprocessing"
include { phasing_beagle } from "./workflows/phasing"
include { postprocessing_ref } from "./workflows/postprocessing"
workflow {

    ch_reference_vcf = Channel.value(file (params.refcsv)) \
    | splitCsv(header:true) \
    | map { row-> tuple (row.chr, file(row.ref_vcf), file(row.ref_vcf_index), file(row.gmap))}
    ch_unphased_vcf = Channel.value(tuple file(params.vcf_unphased), file("${params.vcf_unphased}.csi"))

    preprocessing(ch_unphased_vcf, ch_reference_vcf)
    phasing_beagle(preprocessing.out.unphased, preprocessing.out.ref_1000G)
    postprocessing_ref(phasing_beagle.out.vcf_single.collect())
    postprocessing_ref.out.view()
}