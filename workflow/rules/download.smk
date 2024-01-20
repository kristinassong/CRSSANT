rule prefetch_PARIS:
    output:
        "resources/sra/PARIS/{accession}/{accession}.sra"
    params:
        max_size = "30g",
        out_dir = "resources/sra/PARIS"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Prefetch all necessary files for PARIS {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} --max-size {params.max_size} -O {params.out_dir}"


rule fastq_dump_PARIS:
    input:
        rules.prefetch_PARIS.output
    output:
        "resources/fastq/PARIS/{accession}.fastq.gz"
    params:
        working_dir = "resources/sra/PARIS",
        out_dir = "../../fastq/PARIS"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Download PARIS {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "parallel-fastq-dump --sra-id {wildcards.accession} --threads {threads} "
        "--skip-technical --outdir {params.out_dir} --split-3 --gzip"


rule prefetch_PARIS2:
    output:
        "resources/sra/PARIS2/{accession}/{accession}.sra"
    params:
        out_dir = "resources/sra/PARIS2"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Prefetch all necessary files for PARIS2 {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} -O {params.out_dir}"


rule fastq_dump_PARIS2:
    input:
        rules.prefetch_PARIS2.output
    output:
        "resources/fastq/PARIS2/{accession}.fastq.gz"
    params:
        working_dir = "resources/sra/PARIS2",
        out_dir = "../../fastq/PARIS2"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Download PARIS2 {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "parallel-fastq-dump --sra-id {wildcards.accession} --threads {threads} "
        "--skip-technical --outdir {params.out_dir} --split-3 --gzip"


rule prefetch_LIGR:
    output:
        "resources/sra/LIGR/{accession}/{accession}.sra"
    params:
        out_dir = "resources/sra/LIGR"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Prefetch all necessary files for LIGR {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} -O {params.out_dir}"


rule fastq_dump_LIGR:
    input:
        rules.prefetch_LIGR.output
    output:
        "resources/fastq/LIGR/{accession}.fastq.gz"
    params:
        working_dir = "resources/sra/LIGR",
        out_dir = "../../fastq/LIGR"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Download LIGR {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "parallel-fastq-dump --sra-id {wildcards.accession} --threads {threads} "
        "--skip-technical --outdir {params.out_dir} --split-3 --gzip"


rule prefetch_SPLASH:
    output:
        "resources/sra/SPLASH/{accession}/{accession}.sra"
    params:
        out_dir = "resources/sra/SPLASH"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Prefetch all necessary files for SPLASH {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} -O {params.out_dir}"


rule fastq_dump_SPLASH:
    input:
        rules.prefetch_SPLASH.output
    output:
        "resources/fastq/SPLASH/{accession}.fastq.gz"
    params:
        working_dir = "resources/sra/SPLASH",
        out_dir = "../../fastq/SPLASH"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Download SPLASH {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "parallel-fastq-dump --sra-id {wildcards.accession} --threads {threads} "
        "--skip-technical --outdir {params.out_dir} --split-3 --gzip"