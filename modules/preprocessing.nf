
process fill_AC{
    input:
    tuple path(panel), path(panel_index)
    output:
    tuple path("AC_${panel.simpleName}.bcf"), path("AC_${panel.simpleName}.bcf.csi"), emit: bcf
    script:
    """
    bcftools +fill-AN-AC ${panel} -Ob -o AC_${panel.simpleName}.bcf -W --threads ${task.cpus}
    """
    stub:
    """
    echo "bcftools +fill-AN-AC ${panel}" > AC_${panel.simpleName}.bcf
    touch AC_${panel.simpleName}.bcf.csi
    """
}

process remove_duplicates{
    input:
    tuple path(panel), path(panel_index)
    output:
    tuple path("rm_dup_${panel.simpleName}.bcf"), path("rm_dup_${panel.simpleName}.bcf.csi"), emit: bcf
    script:
    """
    bcftools norm ${panel} -D -Ob -o rm_dup_${panel.simpleName}.bcf -W --threads ${task.cpus}
    """
    stub:
    """
    echo "bcftools norm ${panel} -D -Ob -o rm_dup_${panel.simpleName}.bcf -W --threads ${task.cpus}" > rm_dup_${panel.simpleName}.bcf
    touch rm_dup_${panel.simpleName}.bcf.csi
    """
}

process remove_missing { 
    input:
    tuple path(panel), path(panel_index)
    output:
    tuple path("rm_miss_${panel.simpleName}.bcf"), path("rm_miss_${panel.simpleName}.bcf.csi"), emit: bcf
    script:
    """
    bcftools view -e 'GT[*]="mis"' ${panel} -Ob -o rm_miss_${panel.simpleName}.bcf -W --threads ${task.cpus}
    """
    stub:
    """
    echo "bcftools view -e 'GT[*]=\"mis\"' ${panel} -Ob -o rm_miss_${panel.simpleName}.bcf -W --threads ${task.cpus}" > rm_miss_${panel.simpleName}.bcf
    touch rm_miss_${panel.simpleName}.bcf.csi
    """
}

process create_index {
    input:
    tuple path(panel), path(panel_index)
    output:
    tuple path("id_${panel.simpleName}.vcf.gz"), path("id_${panel.simpleName}.vcf.gz.csi"), emit: bcf
    script:
    """
    bcftools annotate -I "%CHROM:%POS:%REF:%FIRST_ALT" ${panel} -Oz -o id_${panel.simpleName}.vcf.gz -W
    """
    stub:
    """
    touch id_${panel.simpleName}.vcf.gz
    touch id_${panel.simpleName}.vcf.gz.csi
    """
}

process create_index_1000G {
    input:
    tuple val(chr), path(panel), path(panel_index), path(gmap)
    output:
    tuple val(chr), path("id_${panel.simpleName}.vcf.gz"), path("id_${panel.simpleName}.vcf.gz.csi"), path(gmap), emit: bcf
    script:
    """
    bcftools annotate -I "%CHROM:%POS:%REF:%FIRST_ALT" ${panel} -Oz -o id_${panel.simpleName}.vcf.gz -W
    """
    stub:
    """
    touch id_${panel.simpleName}.vcf.gz
    touch id_${panel.simpleName}.vcf.gz.csi
    """
}

