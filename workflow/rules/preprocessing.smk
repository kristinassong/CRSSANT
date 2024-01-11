rules preprocess_PARIS:
    input:
        rules.fasterq_dump_PARIS.output
    output:
        "data"
    threads:
        16
    conda:
        "../envs/trimmomatic.yaml"
    log:
    message:
        "Remove 3' adapters using Trimmomatic."
    shell:
        remove 3' adapters using trimmomatic-0.32.jar SE -threads 16 -phred33
removed PCR duplicates using splitFastqLibrary (in the icSHAPE pipeline)
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "&> {log}"