process concatenate_phased {
    publishDir "$params.outdir/phased_variants/", mode: 'copy'
    input:
    path(vcf_files)
    output:
    tuple path("phased_variants.bcf"), path("phased_variants.bcf.csi"), emit: phased_variants
    script:
    """
    bcftools concat -Ob -o merged_variants.bcf ${vcf_files.collect{it}.join(' ')} -W --threads ${task.cpus}
    bcftools sort -Ob -o phased_variants.bcf merged_variants.bcf -W
    """
    stub:
    """
    touch phased_variants.bcf
    touch phased_variants.bcf.csi
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