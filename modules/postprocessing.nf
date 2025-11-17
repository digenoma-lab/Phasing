process concatenate_phased {
    publishDir "$params.outdir/phased_variants/", mode: 'copy'
    input:
    path(vcf_files)
    val(output_prefix)
    output:
    tuple path("${output_prefix}.bcf"), path("${output_prefix}.bcf.csi"), emit: phased_variants
    script:
    """
    bcftools concat -Ob -o merged_variants.bcf ${vcf_files.collect{it}.join(' ')} -W --threads ${task.cpus}
    bcftools sort -Ob -o ${output_prefix}.bcf merged_variants.bcf -W
    """
    stub:
    """
    touch ${output_prefix}.bcf
    touch ${output_prefix}.bcf.csi
    """                 
}