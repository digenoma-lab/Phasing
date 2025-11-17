include { fix_ploidy; fill_AC; remove_duplicates; remove_missing; create_index;
    create_index_1000G; convert_to_vcf; clean; remove_same_ref_alt; remove_multiallelic_snps } from "../modules/preprocessing"

workflow preprocessing_unphased{
    take:
    unphased
    main:
    create_index(unphased)
    fill_AC(create_index.out.bcf)
    remove_duplicates(fill_AC.out.bcf)
    remove_missing(remove_duplicates.out.bcf)
    emit:
    unphased = remove_missing.out.bcf
}
workflow preprocessing_1000G{
    take:
    ref
    main:
    create_index_1000G(ref)
    emit:
    ref = create_index_1000G.out.bcf
}

workflow preprocessing{
    take:
    unphased
    ref_1000G
    main:
    preprocessing_unphased(unphased)
    preprocessing_1000G(ref_1000G)
    emit:
    unphased = preprocessing_unphased.out.unphased
    ref_1000G = preprocessing_1000G.out.ref
}