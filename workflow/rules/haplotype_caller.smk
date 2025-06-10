rule lifted_over_and_combined_vcfs:
    input:
        shifted_vcf = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/variants.vcf",
        vcf = "{main_dir}/{SRR}/gen_haplotype_caller_variants/variants.vcf",
        shiftback_chain = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shiftback.chain",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        ref_fa = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        rejected_vcf = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.rejected.vcf",
        final_vcf = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.final.vcf",
        final_vcf_index = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.final.vcf.idx"
    log:
        stderr_lift_over = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/lift_over_stderr",
        stdout_lift_over = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/lift_over_stdout",
        stderr_merge = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/merge_stderr",
        stdout_merge = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/merge_stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk LiftoverVcf \
        -I {input.shifted_vcf} \
        -O {main_dir}/{wildcards.SRR}/lifted_over_and_combined_vcfs/{wildcards.SRR}.shifted_back.vcf \
        -R {input.ref_fa} \
        --CHAIN {input.shiftback_chain} \
        --REJECT {output.rejected_vcf} 2> {log.stderr_lift_over} > {log.stdout_lift_over}; \

        gatk MergeVcfs \
        -I {main_dir}/{wildcards.SRR}/lifted_over_and_combined_vcfs/{wildcards.SRR}.shifted_back.vcf \
        -I {input.vcf} \
        -O {output.final_vcf} 2> {log.stderr_merge} > {log.stdout_merge}
        """

rule gen_haplotype_caller_variants:
    input:
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam.bai"
    output:
        variants = "{main_dir}/{SRR}/gen_haplotype_caller_variants/variants.vcf",
        variants_index = "{main_dir}/{SRR}/gen_haplotype_caller_variants/variants.vcf.idx",
        assembled_regions = "{main_dir}/{SRR}/gen_haplotype_caller_variants/assembly_regions.tsv",
        assembled_reads = "{main_dir}/{SRR}/gen_haplotype_caller_variants/assembled_reads.bam",
        assembled_reads_index = "{main_dir}/{SRR}/gen_haplotype_caller_variants/assembled_reads.bai"
    params:
        min_allele_fraction = "0.2",
        sample_ploidy = "1"
    log:
        stderr = "{main_dir}/{SRR}/gen_haplotype_caller_variants/stderr",
        stdout = "{main_dir}/{SRR}/gen_haplotype_caller_variants/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms3000m" HaplotypeCaller \
        --input {input.reads} \
        --output {output.variants} \
        --reference {input.ref} \
        --sample-ploidy {params.sample_ploidy} \
        --read-index {input.reads_index} \
        --assembly-region-out {output.assembled_regions} \
        --bam-output {output.assembled_reads} 2> {log.stderr} > {log.stdout}
        """

use rule gen_haplotype_caller_variants as gen_haplotype_caller_variants_shifted with:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.dict",
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/reads.bam",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/reads.bam.bai"
    output:
        variants = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/variants.vcf",
        variants_index = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/variants.vcf.idx",
        assembled_regions = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/assembly_regions.tsv",
        assembled_reads = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/assembled_reads.bam",
        assembled_reads_index = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/assembled_reads.bai"
    params:
        min_allele_fraction = "0.2",
        sample_ploidy = "1"
    log:
        stderr = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_haplotype_caller_variants_shifted/stdout"
    container:
        "docker://broadinstitute/gatk"
