rule call_variants_via_haplotype_caller:
    input:
        dupmarked_sorted_filtered_reads = "{main_dir}/{SRR}/filter_out_poor_reads/aligned.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        variants = "{main_dir}/{SRR}/call_variants_via_haplotype_caller/variants.vcf.gz"
    params:
        min_allele_fraction = "0.2",
        sample_ploidy = "1"
    log:
        stderr = "{main_dir}/{SRR}/call_variants_via_haplotype_caller/stderr",
        stdout = "{main_dir}/{SRR}/call_variants_via_haplotype_caller/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-version "-Xms3000m" HaplotypeCaller \
        --reference {input.ref} \
        --sample_ploidy {params.sample_ploidy} \
        --input {input.dupmarked_sorted_filtered_reads} \
        --output {output.variants} 2> {log.stderr} > {log.stdout}
        """
