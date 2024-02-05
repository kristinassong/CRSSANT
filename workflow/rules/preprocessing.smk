###############################################################
# PARIS & PARIS2
###############################################################

rule trim3_PARIS:
    input:
        fq_gz = "resources/fastq/PARIS/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/PARIS/{accession}.trim3.fastq"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS/{accession}_trim3.log"
    message:
        "3'-end adapter trimming for PARIS {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "ILLUMINACLIP:resources/adapters.fa:3:20:10 SLIDINGWINDOW:4:20 MINLEN:18"


rule read_collapse_PARIS:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output,
        fq = rules.trim3_PARIS.output
    output:
        "resources/preprocessed_fastq/PARIS/{accession}_trim3_nodup.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from PARIS {wildcards.accession}."
    shell:
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {input.fq} -o {output}"


rule trim5_PARIS:
    input:
        rules.read_collapse_PARIS.output
    output:
        "resources/preprocessed_fastq/PARIS/{accession}_preprocessed.fastq"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS/{accession}_trim5.log"
    message:
        "5'-end adapter trimming for PARIS {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "HEADCROP:17 MINLEN:20"


rule trim3_PARIS2:
    input:
        fq_gz = "resources/fastq/PARIS2/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}.trim3.fastq"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS2/{accession}_trim3.log"
    message:
        "3'-end adapter trimming for PARIS2 {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "ILLUMINACLIP:resources/adapters.fa:3:20:10 SLIDINGWINDOW:4:20 MINLEN:18"


rule read_collapse_PARIS2:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output,
        fq = rules.trim3_PARIS2.output
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}_trim3_nodup.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from PARIS2 {wildcards.accession}."
    shell:
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {input.fq} -o {output}"


rule trim5_PARIS2:
    input:
        rules.read_collapse_PARIS2.output
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}_preprocessed.fastq"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS2/{accession}_trim5.log"
    message:
        "5'-end adapter trimming for PARIS2 {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "HEADCROP:17 MINLEN:20"


###############################################################
# LIGR-seq -- Already preprocessed
###############################################################

rule unzip_LIGR_seq:
    input:
        fq_gz = "resources/fastq/LIGR_seq/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/LIGR_seq/{accession}_preprocessed.fastq"
    message:
        "Unzip LIGR-seq {wildcards.accession} FASTQ file. Already preprocessed."
    shell:
        "gunzip -f -c {input} > {output}"


###############################################################
# SPLASH -- Already preprocessed
###############################################################

rule unzip_SPLASH:
    input:
        fq_gz = "resources/fastq/SPLASH/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/SPLASH/{accession}_preprocessed.fastq"
    message:
        "Unzip SPLASH {wildcards.accession} FASTQ file. Already preprocessed."
    shell:
        "gunzip -f -c {input} > {output}"


###############################################################
# FASTQC
###############################################################

rule fastqc_raw:
    input:
        "resources/fastq/{experiment}/{accession}.fastq.gz"
    output:
        "results/fastqc/raw/{experiment}/{accession}_fastqc.html"
    params:
        "results/fastqc/raw/{experiment}"
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/fastqc/raw/{experiment}/{accession}.log"
    message:
        "Quality control check on the raw sequencing data of {wildcards.experiment} {wildcards.accession}."
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input} "
        "&> {log}"


rule preprocessing_status:
    input:
        expand(rules.unzip_SPLASH.output,accession=config['SPLASH']),
        expand(rules.unzip_LIGR_seq.output,accession=config['LIGR_seq']),
        expand(rules.trim5_PARIS.output,accession=config['PARIS']),
        expand(rules.trim5_PARIS2.output,accession=config['PARIS2'])
    output:
        "resources/preprocessed_fastq/preprocessing_status.txt"
    shell:
        "echo date > {output} && "
        "echo \'All preprocessing steps complete. Run FASTQC.\' >> {output}"


rule fastqc_preprocessed:
    input:
        status = rules.preprocessing_status.output,
        fq = "resources/preprocessed_fastq/{experiment}/{accession}_preprocessed.fastq"
    output:
        "results/fastqc/preprocessed/{experiment}/{accession}_fastqc.html"
    params:
        "results/fastqc/preprocessed/{experiment}"
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/fastqc/preprocessed/{experiment}/{accession}.log"
    message:
        "Quality control check on the preprocessed sequencing data of {wildcards.experiment} {wildcards.accession}."
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "{input.fq} "
        "&> {log}"