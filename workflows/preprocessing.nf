include {fill_AC; remove_duplicates; remove_missing; create_index;
    create_index_1000G; fasta_index; reheader; get_contigs; fix_ploidy} from "../modules/preprocessing"

workflow preprocessing_unphased{
    take:
    unphased
    main:
    fix_ploidy(unphased)
    create_index(fix_ploidy.out.vcf)
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