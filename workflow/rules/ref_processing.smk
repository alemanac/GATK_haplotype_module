#generates indexes for ref

rule gen_ref_bwa_indexes:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.amb",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.bwt",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.ann",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.sa"
    params:
        ref_fa = f"{reference_fasta}",
        ref_base = f"{ref_base}",
        rule_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_stdout"
    container:
        "file:///vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: replace this with an online-hosted image
    shell:
        """
        cd {params.rule_dir}; \
        
        bwa index {params.ref_fa} 2> ../../../../{log.stderr} > ../../../../{log.stdout}
        """

"""
creates BWA index for shifted fasta
"""
use rule gen_ref_bwa_indexes as gen_ref_bwa_indexes_shifted with: 
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna"
    params:
        ref_fa = f"{ref_base}.shifted.fna",
        ref_base = f"{ref_base}.shifted",
        rule_dir = f"{main_dir}/{ref_base}/ref_processing"
    output:
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.amb",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.ann",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.bwt",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.sa"
    log:
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_shifted_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_shifted_stdout"
    container:
        "file:///vol/patric3/production/containers/ubuntu-045-12.sif" #TODO: replace this an online-hosted image

"""
creates shifted fasta file, index, and other supplemental files for such shift
"""
rule shift_fasta:
    input:
        ref_index = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        shifted_ref_fa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna",
        shifted_ref_fa_index = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.fai",
        shiftback_chain = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shiftback.chain",
        shifted_intervals = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.intervals",
        unshifted_intervals = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.intervals",
        shifted_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.dict"
    params:
        rule_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_shift_fasta_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_shift_fasta_stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk ShiftFasta \
        -R {input.ref} \
        --read-index {input.ref_index} \
        -O {output.shifted_ref_fa} \
        --interval-file-name {params.rule_dir}/{ref_base} \
        --shift-back-output {output.shiftback_chain} 2> {log.stderr} > {log.stdout}
        """

rule gen_ref_faidx:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai"
    log:
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_faidx_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_faidx_stdout",
    container: "file:///vol/patric3/production/containers/ubuntu-045-12.sif"
    shell:
        """
        samtools faidx {input.ref} \
        -o {output.fai} 2> {log.stderr} > {log.stdout}
        """

rule gen_ref_dict:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict"
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_dict_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_dict_stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk CreateSequenceDictionary \
        -R {config[ref_fna]} \
        -O {output.ref_dict} 2> {log.stderr} > {log.stdout}; \
        """    

rule move_ref:
    output:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    shell:
        """
        cp {config[ref_fna]} {output.ref}
        """
rule move_ref_python_pipeline:
    output:
        ref = f"{main_dir}/run_gatk_python_pipeline/{reference_fasta}"