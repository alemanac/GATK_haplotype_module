reference_fasta = Path(config['ref_fna']).name
ref_base = Path(reference_fasta).with_suffix('')

#TODO: trim with cutadapt as specified in the paper

rule filter_out_poor_reads:
    input:
        dupmarked_sorted_reads = "{main_dir}/{SRR}/gen_sorted_marked_duplicates_bam/{SRR}_dupmarked_sorted.bam"
    output:
        dupmarked_sorted_filtered_reads = "{main_dir}/{SRR}/filter_out_poor_reads/aligned.bam"
    params:
        minimum_mq = "10"
    log:
        stderr = "{main_dir}/{SRR}/filter_out_poor_reads/stderr",
        stdout = "{main_dir}/{SRR}/filter_out_poor_reads/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        samtools view \
        --min-MQ 10 \
        --bam \
        --output {output.dupmarked_sorted_filtered_reads} \
        {input.dupmarked_sorted_reads} >2 {log.stderr} > {log.stdout}
        """


#TODO: it shouldn't matter, but consider running this before mark duplicates in the pipeline
#TODO: test indexing with samtools, although it shouldn't matter
"""
I'm unsure what this is sorting, but I think it's the indexes for the reads
"""
rule gen_sorted_marked_duplicates_bam:
    input:
        dupmarked_reads = "{main_dir}/{SRR}/gen_marked_duplicates_bam/{SRR}_dupmarked.bam"
    output:
        dupmarked_sorted_reads = "{main_dir}/{SRR}/gen_sorted_marked_duplicates_bam/{SRR}_dupmarked_sorted.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_sorted_marked_duplicates_bam/logs/{SRR}_generate_sorted_marked_duplicates_bam_stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_marked_duplicates_bam/logs/{SRR}_generate_sorted_marked_duplicates_bam_stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms3000m" SortSam \
        INPUT={input.dupmarked_reads} \
        OUTPUT={output.dupmarked_sorted_reads} \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000 2> {log.stderr} > {log.stdout}
        """
        "docker://broadinstitute/gatk"

#TODO: consider playing with the optical_duplicate_pixel_distance
rule gen_marked_duplicates_bam:
    input:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_bam/{SRR}.bam"
    output:
        dupmarked_reads = "{main_dir}/{SRR}/gen_marked_duplicates_bam/{SRR}_dupmarked.bam",
        metric_file = "{main_dir}/{SRR}/gen_marked_duplicates_bam/metrics.txt"
    log: 
        stderr = "{main_dir}/{SRR}/gen_marked_duplicates_bam/{SRR}_generate_marked_duplicates_bam.stderr",
        stdout = "{main_dir}/{SRR}/gen_marked_duplicates_bam/{SRR}_generate_marked_duplicates_bam.stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xms3000m" MarkDuplicates \
        INPUT={input.aligned_reads} \
        OUTPUT={output.dupmarked_reads} \
        METRICS_FILE={output.metric_file} \
        CLEAR_DT="false" 2> {log.stderr} > {log.stdout}
        """

"""
generate index (.bai) for unmapped reads (.bam, below)
"""
rule create_unmapped_bai:
    input:
        bam = "{main_dir}/{SRR}/gen_unmapped_bam/{SRR}.unmapped.bam"
    output:
        bai = "{main_dir}/{SRR}/create_umapped_bai/{SRR}.unmapped.bai"
    log:
        stderr = "{main_dir}/{SRR}/create_umapped_bai/logs/{SRR}_create_unmapped_bai_stderr",
        stdout = "{main_dir}/{SRR}/create_umapped_bai/logs/{SRR}_create_unmapped_bai_stdout"
    container:
        "docker://broadinstitute/gatk" 
    shell:
        """
        samtools index {input.bam} {output.bai} 2> {log.stderr} > {log.stdout}
        """

"""
generates index (.bai) for mapped sorted reads (.bam. below)
"""
use rule create_unmapped_bai as create_mapped_sorted_bai with:
    input:
        bam = "{main_dir}/{SRR}/gen_sorted_mapped_bam/{SRR}.mapped.sorted.bam"
    output:
        bai = "{main_dir}/{SRR}/create_mapped_sorted_bai/{SRR}.mapped.sorted.bai"
    log:
        stderr = "{main_dir}/{SRR}/create_mapped_sorted_bai/logs/{SRR}_create_mapped_bai_stderr",
        stdout = "{main_dir}/{SRR}/create_mapped_sorted_bai/logs/{SRR}_create_mapped_bai_stdout"
    container:
        "docker://broadinstitute/gatk" 

"""
sorts mapped reads (.bam)
"""
rule gen_sorted_mapped_bam:
    input:
        bam = "{main_dir}/{SRR}/gen_aligned_bam/{SRR}.mapped.bam"
    output:
        sorted_bam =  "{main_dir}/{SRR}/gen_sorted_mapped_bam/{SRR}.mapped.sorted.bam"
    log:
        stderr = "{main_dir}/{SRR}/gen_sorted_mapped_bam/logs/{SRR}_generate_sorted_mapped_bam_stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_mapped_bam/logs/{SRR}_generate_sorted_mapped_bam_stdout"
    container:
        "docker://broadinstitute/gatk" 
    shell:
        """
        samtools sort {input.bam} -o {output.sorted_bam} 2> {log.stderr} > {log.stdout}
        """

#TODO: removed the -k flag. perform the same test but with the k flag
#TODO": removed the -y flag. perform a test with it
rule gen_aligned_bam:
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
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_bam/{SRR}.bam"
    params:
        ref_processing_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = "{main_dir}/{SRR}/gen_aligned_bam/logs/stderr",
        stdout = "{main_dir}/{SRR}/gen_aligned_bam/logs/stdout"
    container: 
        "file:///vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: replace this an online-hosted image
    shell:
        """
        cd {params.ref_processing_dir}; \

        bwa mem \
        -v 3 \
        -t 2 \
        ../../../{input.ref} ../../../{input.fq_1} ../../../{input.fq_2} \
        -o ../../../{output.aligned_reads} 2> ../../../{log.stderr} > ../../../{log.stdout}
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