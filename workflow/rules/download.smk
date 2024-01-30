rule prefetch:
    output:
        "resources/sra/{experiment}/{accession}/{accession}.sra"
    params:
        max_size = "30g",
        out_dir = "resources/sra/{experiment}"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Prefetch all necessary files for {wildcards.experiment} {wildcards.accession}."
    shell:
        "prefetch {wildcards.accession} --max-size {params.max_size} -O {params.out_dir}"


rule fastq_dump:
    input:
        rules.prefetch.output
    output:
        "resources/fastq/{experiment}/{accession}.fastq.gz"
    params:
        working_dir = "resources/sra/{experiment}",
        out_dir = "../../fastq/{experiment}"
    conda:
        "../envs/sra_toolkit.yaml"
    message:
        "Download {wildcards.experiment} {wildcards.accession} fastq files from GEO."
    threads:
        6
    shell:
        "cd {params.working_dir} && "
        "parallel-fastq-dump --sra-id {wildcards.accession} --threads {threads} "
        "--skip-technical --outdir {params.out_dir} --split-3 --gzip"


rule download_icSHAPE_git:
    output:
        dir = directory('resources/icSHAPE'),
        script = "resources/icSHAPE/icSHAPE_pipeline.pl"
    conda:
        '../envs/git.yaml'
    message:
        "Download icSHAPE git repository."
    shell:
        "git clone git@github.com:qczhang/icSHAPE.git {output}"