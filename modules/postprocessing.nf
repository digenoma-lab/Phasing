process concatenate_phased {
    publishDir "$params.outdir/phased_variants/", mode: 'copy'
    input:
    path(vcf_files)
    val(output_prefix)
    tuple path(fasta), path(fasta_index)
    output:
    tuple path("${output_prefix}.vcf.gz"), path("${output_prefix}.vcf.gz.csi"), emit: phased_variants
    script:
    """
    bcftools concat --naive \\
        -Oz -o merged_variants.vcf.gz ${vcf_files.collect{it}.join(' ')} --threads ${task.cpus}
    bcftools reheader merged_variants.vcf.gz -f ${fasta_index} -o temp.vcf.gz
    bcftools sort -Oz -o ${output_prefix}.vcf.gz temp.vcf.gz -W
    """
    stub:
    """
    touch ${output_prefix}.vcf.gz
    touch ${output_prefix}.vcf.gz.csi
    """
}