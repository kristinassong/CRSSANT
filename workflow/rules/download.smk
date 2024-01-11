rule prefetch_PARIS:
    output:
        directory("resources/sra/PARIS/{accession}")
    params:
        max_size = "30g",
        out_dir = "resources/sra/PARIS"
    conda:
        "../envs/sra_toolkit.yaml"
    log:
        "results/logs/prefetch/PARIS/{accession}.log"
    message:
        "Prefetch all necessary files for PARIS {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} --max-size {params.max_size} -O {params.out_dir} &> {log}"


rule fasterq_dump_PARIS:
    input:
        rules.prefetch_PARIS.output
    output:
        "resources/fastq/PARIS/{accession}.fastq"
    params:
        working_dir = "resources/sra/PARIS",
        out_dir = "../../fastq/PARIS"
    conda:
        "../envs/sra_toolkit.yaml"
    log:
        "results/logs/fasterq_dump/PARIS/{accession}.log"
    message:
        "Download PARIS {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "fasterq-dump {wildcards.accession} --outdir {params.out_dir} &> {log}"