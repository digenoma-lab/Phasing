process fix_ploidy {
    input:
    tuple path(array), path(array_index)
    output:
    tuple path("ploidy_${array.simpleName}.vcf.gz"), path("ploidy_${array.simpleName}.vcf.gz.csi"), emit: diploid_array
    script:
    """
    echo "chrX 1 156040895 M 2" > ploidy.txt
    bcftools +fixploidy ${array}  --threads ${task.cpus} -Ov -- -p ploidy.txt | bcftools view -i 'ALT !="-" & ALT !="."' --threads ${task.cpus} | bcftools sort -Oz -o ploidy_${array.simpleName}.vcf.gz -W 
    """
    stub:
    """
    echo "chrX 1 156040895 M 2" > ploidy.txt
    echo "bcftools +fixploidy ${array} --threads ${task.cpus} -Ov -- -p ploidy.txt| bcftools view -i 'ALT !=\"-\" & ALT !=\".\"' --threads ${task.cpus} | bcftools sort -Ob -o GWAS_diploid.bcf -W"  > ploidy_${array.simpleName}.vcf.gz
    touch ploidy_${array.simpleName}.vcf.gz.csi
    """
}

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

process bcftools_index{
    input:
    path(bcf)
    output:
    tuple path(bcf), path("${bcf}.csi"), emit: bcf_index
    script:
    """
    bcftools index ${bcf} -o ${bcf}.csi --threads ${task.cpus} -f
    """
    stub:
    """
    echo bcftools index ${bcf} -O z -o ${bcf}.csi --threads ${task.cpus} > ${bcf}.csi
    """
}

process convert_to_vcf{
    input:
    tuple val(chr), path(ref_vcf), path(ref_vcf_index), path(gmap)
    output:
    tuple val(chr), path("${ref_vcf.simpleName}.vcf.gz"), path(ref_vcf_index), path(gmap), emit: vcf
    script:
    """
    bcftools query -l ${ref_vcf} | shuf -n 100 > sample_ids.txt
    bcftools query -f '%ID\n' ${ref_vcf} | shuf -n 1000000 > snps.txt
    bcftools view ${ref_vcf} -S sample_ids.txt --include ID=@snps.txt -Oz -o sample_${ref_vcf.simpleName}.vcf.gz
    """
    stub:
    """
    echo " bcftools convert ${ref_vcf} -Ov -o ${ref_vcf.simpleName}.vcf --threads ${task.cpus} " > sample_${ref_vcf.simpleName}.vcf
    """
}

process clean{
    input:
    tuple path(panel), path(panel_index)
    output:
    tuple path("clean_${panel.simpleName}.vcf.gz"), path("clean_${panel.simpleName}.vcf.gz.csi"), emit: bcf
    script:
    """
    bcftools view -i 'ALT !="-" & ALT !="."' ${panel} | bcftools sort -Oz -o clean_${panel.simpleName}.vcf.gz -W
    """
    stub:
    """
    touch clean_${panel.simpleName}.vcf.gz
    touch clean_${panel.simpleName}.vcf.gz.csi
    """
}

process remove_same_ref_alt {
    input:
    tuple path(bcf), path(bcf_index)
    output:
    tuple path("filtered_${bcf.simpleName}.vcf.gz"), path("filtered_${bcf.simpleName}.vcf.gz.csi"), emit: bcf

    script:
    """
    bcftools view -e 'REF=ALT' ${bcf} -Oz -o filtered_${bcf.simpleName}.vcf.gz
    bcftools index filtered_${bcf.simpleName}.vcf.gz
    """
    stub:
    """
    touch filtered_${bcf.simpleName}.vcf.gz
    touch filtered_${bcf.simpleName}.vcf.gz.csi
    """
}

process remove_multiallelic_snps {
    publishDir "$params.outdir/array/", mode: "copy"
    input:
    tuple path(bcf), path(bcf_index)
    output:
    tuple path("biallelic_${bcf.simpleName}.vcf.gz"), path("biallelic_${bcf.simpleName}.vcf.gz.csi"), emit: bcf

    script:
    """
    bcftools view -m2 -M2 -v snps ${bcf} -Oz -o biallelic_${bcf.simpleName}.vcf.gz
    bcftools index biallelic_${bcf.simpleName}.vcf.gz
    """
    stub:
    """
    touch biallelic_${bcf.simpleName}.vcf.gz
    touch biallelic_${bcf.simpleName}.vcf.gz.csi
    """
}

