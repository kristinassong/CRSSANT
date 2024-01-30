######################
# Preprocess PARIS
######################

rule read_collapse_PARIS:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output.dir,
        fq_gz = "resources/fastq/PARIS/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/PARIS/{accession}_read_collapse.fastq"
    params:
        fq = "resources/fastq/PARIS/{accession}.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from PARIS {wildcards.accession}."
    shell:
        "gunzip -f -c {input.fq_gz} > {params.fq} && "
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {params.fq} -o {output} "
        "&& rm {params.fq}"


rule simple_trim_PARIS:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output.dir,
        fq = rules.read_collapse_PARIS.output
    output:
        "resources/preprocessed_fastq/PARIS/{accession}_simple_trim.fastq"
    params:
        head_crop = 13
    log:
        "results/logs/simple_trim/PARIS/{accession}.log"
    message:
        "Adapter trimming for PARIS {wildcards.accession} - Step 1."
    shell:
        "{input.icSHAPE_dir}/bin/simpleTrim -U {input.fq} -o {output} -d {params.head_crop} "
        "&> {log}"


rule trimmomatic_PARIS:
    input:
        rules.simple_trim_PARIS.output
    output:
        "resources/preprocessed_fastq/PARIS/{accession}.fastq"
    params:
        adapters_dir = "resources/icSHAPE/data/adapter/"
    threads:
        8
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS/{accession}.log"
    message:
        "Adapter trimming for PARIS {wildcards.accession} - Step 2 (via Trimmomatic)."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "ILLUMINACLIP:{params.adapters_dir}/TruSeq2-SE.fa:2:30:4 TRAILING:20 HEADCROP:6 MINLEN:25"


######################
# Preprocess PARIS2
######################

rule read_collapse_PARIS2:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output.dir,
        fq_gz = "resources/fastq/PARIS2/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}_read_collapse.fastq"
    params:
        fq = "resources/fastq/PARIS2/{accession}.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from PARIS2 {wildcards.accession}."
    shell:
        "gunzip -f -c {input.fq_gz} > {params.fq} && "
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {params.fq} -o {output} "
        "&& rm {params.fq}"


rule simple_trim_PARIS2:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output.dir,
        fq = rules.read_collapse_PARIS2.output
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}_simple_trim.fastq"
    params:
        head_crop = 13
    log:
        "results/logs/simple_trim/PARIS2/{accession}.log"
    message:
        "Adapter trimming for PARIS2 {wildcards.accession} - Step 1."
    shell:
        "{input.icSHAPE_dir}/bin/simpleTrim -U {input.fq} -o {output} -d {params.head_crop} "
        "&> {log}"


rule trimmomatic_PARIS2:
    input:
        rules.simple_trim_PARIS2.output
    output:
        "resources/preprocessed_fastq/PARIS2/{accession}.fastq"
    params:
        adapters_dir = "resources/icSHAPE/data/adapter/"
    threads:
        8
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/PARIS2/{accession}.log"
    message:
        "Adapter trimming for PARIS2 {wildcards.accession} - Step 2 (via Trimmomatic)."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "ILLUMINACLIP:{params.adapters_dir}/TruSeq2-SE.fa:2:30:4 TRAILING:20 HEADCROP:6 MINLEN:25"


######################
# Preprocess LIGR-seq
######################

rule read_collapse_LIGR:
    input:
        icSHAPE_dir = rules.download_icSHAPE_git.output.dir,
        fq_gz = "resources/fastq/LIGR/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/LIGR/{accession}_read_collapse.fastq"
    params:
        fq = "resources/fastq/LIGR/{accession}.fastq"
    conda:
        "../envs/git.yaml"
    message:
        "Remove PCR duplicates from LIGR-seq {wildcards.accession}."
    shell:
        "gunzip -f -c {input.fq_gz} > {params.fq} && "
        "perl {input.icSHAPE_dir}/scripts/readCollapse.pl -U {params.fq} -o {output} "
        "&& rm {params.fq}"


rule trimmomatic_LIGR:
    input:
        rules.read_collapse_LIGR.output
    output:
        "resources/preprocessed_fastq/LIGR/{accession}.fastq"
    params:
        adapters_dir = "resources/icSHAPE/data/adapter/"
    threads:
        8
    conda:
        "../envs/trimmomatic.yaml"
    log:
        "results/logs/trimmomatic/LIGR/{accession}.log"
    message:
        "Adapters and barcodes trimming for LIGR-seq {wildcards.accession} via Trimmomatic."
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "-trimlog {log} "
        "{input} {output} "
        "HEADCROP:5 "
        "ILLUMINACLIP:{params.adapters_dir}/adapters/TruSeq3-SE.fa:2:30:4 "
        "TRAILING:20 "
        "MINLEN:25"


######################
# Preprocess SPLASH
######################

rule unzip_SPLASH:
    input:
        "resources/fastq/SPLASH/{accession}.fastq.gz"
    output:
        "resources/preprocessed_fastq/SPLASH/{accession}.fastq"
    shell:
        "gunzip -f -c {input} > {output}"