reference_fasta = Path(config['ref_fna']).name
ref_base = Path(reference_fasta).with_suffix('')

use rule gen_mq_filtered_reads as gen_mq_filtered_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/reads.bam",
    output:
        mq_filtered_reads = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/reads.bam",
        reads_index = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/reads.bam.bai"
    params:
        minimum_mq = "10"
    log:
        stderr = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/stdout"
    container:
        "docker://broadinstitute/gatk"

use rule gen_duplicates_marked_reads as gen_duplicates_marked_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_read_group_added_reads_shifted/reads.bam"
    output:
        duplicates_marked_reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/reads.bam",
        metric_file = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/metrics.txt"
    log: 
        stderr = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/stdout"
    container:
        "docker://broadinstitute/gatk"

use rule gen_read_group_added_reads as gen_read_group_added_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_sorted_reads_shifted/reads.bam"
    output:
        read_group_added_reads = "{main_dir}/{SRR}/gen_read_group_added_reads_shifted/reads.bam"
    log:
        stderr = "{main_dir}/{SRR}/gen_read_group_added_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_read_group_added_reads_shifted/stdout"
    container: "docker://broadinstitute/gatk"

use rule gen_sorted_reads as gen_sorted_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_aligned_reads_shifted/{SRR}.bam"
    output:
        sorted_reads = "{main_dir}/{SRR}/gen_sorted_reads_shifted/reads.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_sorted_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_reads_shifted/stdout"
    container: "docker://broadinstitute/gatk"

use rule gen_aligned_reads as gen_aligned_reads_shifted with:
    input:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq",
        ref = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna",
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.amb",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.ann",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.bwt",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.sa"
    output:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads_shifted/{SRR}.bam"
    params:
        ref_processing_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = "{main_dir}/{SRR}/gen_aligned_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_aligned_reads_shifted/stdout"
    container: 
        "file:///vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: replace this an online-hosted image
