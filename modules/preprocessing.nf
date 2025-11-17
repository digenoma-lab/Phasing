
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
    tuple path("rm_miss_${panel.simpleName}.vcf.gz"), path("rm_miss_${panel.simpleName}.vcf.gz.csi"), emit: bcf
    script:
    """
    bcftools view -e 'GT[*]="mis"' ${panel} -Oz -o rm_miss_${panel.simpleName}.vcf.gz -W --threads ${task.cpus}
    """
    stub:
    """
    echo "bcftools view -e 'GT[*]=\"mis\"' ${panel} -Ob -o rm_miss_${panel.simpleName}.vcf.gz -W --threads ${task.cpus}" > rm_miss_${panel.simpleName}.vcf.gz
    touch rm_miss_${panel.simpleName}.vcf.gz.csi
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

process fasta_index {
    input:
    path(fasta)
    output:
    tuple path(fasta), path("${fasta}.fai"), emit: index
    script:
    """
    samtools faidx ${fasta}
    """
    stub:
    """
    touch ${fasta}.fai
    """
}

process get_contigs {
    input:
    tuple path(fasta), path(fasta_index)
    output:
    path ("contigs.txt"), emit: contig
    script:
    """
    awk '{print "##contig=<ID="\$1",length="\$2">"}' ${fasta_index} > contigs.txt
    """
    stub:
    """
    touch contigs.txt
    """
}

process reheader {
    input:
    tuple path(fasta), path(fasta_index)
    tuple path(vcf), path(vcf_index)
    output:
    tuple path("header_${vcf}"), path("header_${vcf}.tbi"), emit: vcf
    script:
    """
    bcftools reheader ${vcf} -f ${fasta_index} -o header_${vcf}
    tabix -p vcf header_${vcf}
    """
    stub:
    """
    touch header_${vcf}
    touch header_${vcf}.tbi
    """
}