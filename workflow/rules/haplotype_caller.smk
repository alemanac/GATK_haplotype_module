rule call_variants_via_haplotype_caller:
    input:
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam.bai"
    output:
        variants = "{main_dir}/{SRR}/call_variants_via_haplotype_caller/variants.vcf"
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
        gatk --java-options "-Xms3000m" HaplotypeCaller \
        --input {input.reads} \
        --output {output.variants} \
        --reference {input.ref} \
        --sample-ploidy {params.sample_ploidy} \
        --read-index {input.reads_index} 2> {log.stderr} > {log.stdout}
        """

rule gen_mutect2_vcfs: 
    input:
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam.bai"
    output:
        variants = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf.stats",
        local_assemblies = "{main_dir}/{SRR}/gen_mutect2_vcfs/locally_assemblies.bam",
        assembly_region_out = "{main_dir}/{SRR}/gen_mutect2_vcfs/assembly_region_out.tsv",
        graph_out = "{main_dir}/{SRR}/gen_mutect2_vcfs/graph_out"
    threads: 8
    log:
        stderr = "{main_dir}/{SRR}/gen_mutect2_vcfs/stderr",
        stdout = "{main_dir}/{SRR}/gen_mutect2_vcfs/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell: 
        """
        gatk --java-options "-Xmx3000m" Mutect2 \
        -R {input.ref} \
        -I {input.reads} \
        -O {output.variants} \
        --bam-output {output.local_assemblies} \
        --annotation StrandBiasBySample \
        --num-matching-bases-in-dangling-end-to-recover 1 \
        --max-reads-per-alignment-start 75 \
        --native-pair-hmm-threads 8 \
        --assembly-region-out {output.assembly_region_out} \
        --graph-output {output.graph_out} 2> {log.stderr} > {log.stdout}
        """

rule gen_filtered_vcfs:
    input:
        variants = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf.stats",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        filtered_variants = "{main_dir}/{SRR}/gen_filtered_vcfs/filtered.vcf"
    params:
        allele_fraction_thres = "0.0",
        microbial_mode = "--microbial-mode"
    log: 
        stderr = "{main_dir}/{SRR}/gen_filtered_vcfs/stderr",
        stdout = "{main_dir}/{SRR}/gen_filtered_vcfs/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms2500m" FilterMutectCalls \
        -V {input.variants} \
        -R {input.ref} \
        -O {output.filtered_variants} \
        --stats {input.stats} \
        {params.microbial_mode}\
        --min-allele-fraction {params.allele_fraction_thres} 2> {log.stderr} > {log.stdout}
        """	
