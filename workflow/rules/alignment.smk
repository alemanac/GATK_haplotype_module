reference_fasta = Path(config['ref_fna']).name
ref_base = Path(reference_fasta).with_suffix('')

#TODO: trim with cutadapt as specified in the paper

rule gen_bam_index:
    input:
        reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam"
    output:
        bam_index = "{main_dir}/{SRR}/gen_bam_index/reads.bai"
    log:
        stderr = "{main_dir}/{SRR}/gen_bam_index/stderr",
        stdout = "{main_dir}/{SRR}/gen_bam_index/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        samtools index -b {input.reads} {output.bam_index} 2> {log.stderr} > {log.stdout}
        """

rule gen_mq_filtered_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam"
    output:
        mq_filtered_reads = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam.bai"
    params:
        minimum_mq = "10"
    log:
        stderr = "{main_dir}/{SRR}/gen_mq_filtered_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_mq_filtered_reads/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        samtools view \
        --min-MQ 10 \
        --bam \
        --output {output.mq_filtered_reads} \
        {input.reads} 2> {log.stderr} > {log.stdout}; \

        samtools index {output.mq_filtered_reads} 2>> {log.stderr} > {log.stdout}
        """

#TODO: consider playing with the optical_duplicate_pixel_distance
rule gen_duplicates_marked_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_read_group_added_reads/reads.bam"
    output:
        duplicates_marked_reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam",
        metric_file = "{main_dir}/{SRR}/gen_duplicates_marked_reads/metrics.txt"
    log: 
        stderr = "{main_dir}/{SRR}/gen_duplicates_marked_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_duplicates_marked_reads/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms3000m" MarkDuplicates \
        INPUT={input.reads} \
        OUTPUT={output.duplicates_marked_reads} \
        METRICS_FILE={output.metric_file} \
        CLEAR_DT="false" 2> {log.stderr} > {log.stdout}
        """

#TODO: the A in @RG may be wrong
rule gen_read_group_added_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_sorted_reads/reads.bam"
    output:
        read_group_added_reads = "{main_dir}/{SRR}/gen_read_group_added_reads/reads.bam"
    log:
        stderr = "{main_dir}/{SRR}/gen_read_group_added_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_read_group_added_reads/stdout"
    container: "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms3000m" AddOrReplaceReadGroups \
        --INPUT {input.reads} \
        --OUTPUT {output.read_group_added_reads} \
        --RGLB {config[read_group_library]} \
        --RGPL {config[read_group_platform]} \
        --RGPU {config[read_group_platform_unit]} \
        --RGSM {wildcards.SRR} 2> {log.stderr} > {log.stdout}
        """


#TODO: it shouldn't matter, but consider running this before mark duplicates in the pipeline
#TODO: test indexing with samtools, although it shouldn't matter
"""
I'm unsure what this is sorting, but I think it's the indexes for the reads
"""
rule gen_sorted_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_aligned_reads/{SRR}.bam"
    output:
        sorted_reads = "{main_dir}/{SRR}/gen_sorted_reads/reads.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_sorted_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_reads/stdout"
    container: "docker://broadinstitute/gatk"
    shell:
        """
        samtools sort {input.reads} \
        -o {output.sorted_reads} 2> {log.stderr} > {log.stdout}
        """

#TODO: removed the -k flag. perform the same test but with the k flag
#TODO": removed the -y flag. perform a test with it
rule gen_aligned_reads:
    input:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.amb",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.bwt",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.ann",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.sa",
    output:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads/{SRR}.bam"
    params:
        ref_processing_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = "{main_dir}/{SRR}/gen_aligned_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_aligned_reads/stdout"
    container: 
        "file:///vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: replace this an online-hosted image
    shell:
        """
        cd {params.ref_processing_dir}; \

        bwa mem \
        -v 3 \
        -t 2 \
        ../../../../{input.ref} ../../../../{input.fq_1} ../../../../{input.fq_2} \
        -o ../../../../{output.aligned_reads} 2> ../../../../{log.stderr} > ../../../../{log.stdout}
        """

#TODO: turn this into a wrapper :3
rule fetch_reads:
    output:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq"
    log: 
        stderr = "{main_dir}/{SRR}/fetch_reads/{SRR}_fetch_reads_stderr",
        stdout = "{main_dir}/{SRR}/fetch_reads/{SRR}_fetch_reads_stdout"
    params:
        rule_dir = "{main_dir}/{SRR}/fetch_reads"
    container: "/vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: consider changing to public image
    shell:
        """
        fasterq-dump {wildcards.SRR} \
        -O {params.rule_dir} 2> {wildcards.SRR}_log_stderr > {wildcards.SRR}_log_stdout; \

        mv {wildcards.SRR}_log_stderr {log.stderr}; mv {wildcards.SRR}_log_stdout {log.stdout}
        """