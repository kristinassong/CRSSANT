rule prefetch_PARIS:
    output:
        "resources/sra/{accession}/{accession}.sra"
    conda:
        "../envs/sra_toolkit.yaml"
    log:
        "results/logs/prefetch/PARIS_{accession}.log"
    message:
        "Download all necessary files for {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} &> {log}"

"""
rule fasterq_dump_PARIS:
    output:
        "resources/{accession}.fastq"
    log:
        "logs/{accession}.log"
    message:
        "Download PARIS fastq files from GEO."
    threads: 6  # defaults to 6
    shell:
        "fasterq-dump "
"""