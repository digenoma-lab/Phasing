process beagle{
    tag "$chr"
    input:
    tuple path(unphased), path(unphasedindex)
    tuple val(chr), path(refvcf), path(refvcfindex), path(map)
    output:
    tuple val(chr), path ("phased_${chr}.vcf.gz"), path("phased_${chr}.vcf.gz.csi"), path(refvcf), path(refvcfindex), path(map), emit: vcf
    path("phased_${chr}.vcf.gz"), emit: vcf_single
    script:
    def memory = task.memory.getGiga()
    """
    bcftools view ${unphased} --regions ${chr} -Ov -o ${unphased.simpleName}.vcf --threads ${task.cpus}
    bcftools view ${refvcf} -Ov -o ${refvcf.simpleName}.vcf --threads ${task.cpus}
    awk '{print "chr"\$0}' ${map} > ${map}.chr
    beagle -Xmx${memory}g gt=${unphased.simpleName}.vcf \
        out=phased_${unphased.simpleName} impute=false ref=${refvcf.simpleName}.vcf \
        nthreads=${task.cpus} chrom=${chr} map=${map}.chr
    bcftools index phased_${unphased.simpleName}.vcf.gz --threads ${task.cpus}
    """
    stub:
    """
    touch phased_${chr}.vcf.gz
    touch phased_${chr}.vcf.gz.csi
    """
}
