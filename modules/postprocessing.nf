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
process publish {
    publishDir "$params.outdir/${output_folder}/", mode:"copy"
    input:
    tuple val(chr), path(input_file), path(input_file_index), path(refvcf), path(refvcfindex), path(map)
    val(output_folder)
    output:
    path(input_file)
    script:
    """
    echo "Publishing"
    """
    stub:
    """
    echo "Publishing"
    """
}