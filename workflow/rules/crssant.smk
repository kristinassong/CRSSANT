rule combine_gap1filter_trans:
    input:
        gap1filter = rules.filter_spliced_short_gaps_gap1.output,
        trans = rules.merge_sam.output.trans
    output:
        "results/crssant/{experiment}/{accession}/{accession}_pri_crssant.sam"
    message:
        "Combine {wildcards.experiment} {wildcards.accession} gap1filter.sam and trans.sam to produce alignfile SAM required for crssant.py."
    script:
        "../scripts/merger.py"


rule create_genesfile:
    input:
        rules.coco_ca.output.gtf_corrected
    output:
        "resources/genes.bed"
    message:
        "Create a bed file that contains a list of all genes in the genome."
    shell:
        "grep -P \'\tgene\t\' {input} | " 
        "cut -f1,4,5,7,9 | sed \'s/[[:space:]]/\t/g\' | "
        "sed \'s/[;|\"]//g\' | awk -F $\'\t\' \'BEGIN {{OFS=FS}} {{print $1,$2-1,$3,$8,\"1000\",$4}}\' | "
        "sed -n \'/^[0-9,X,Y,MT]/Ip\' | awk \'$1 !~ /_/\' | "
        "sort -k1,1 -k2,2n > {output}"


rule sam_to_sorted_bam:
    input:
        rules.combine_gap1filter_trans.output
    output:
        "results/crssant/{experiment}/{accession}/{accession}_pri_crssant_sorted.bam"
    params:
        tmp_bam = "results/crssant/{experiment}/{accession}/{accession}_pri_crssant.bam"
    conda:
        "../envs/mapping.yaml"
    message:
        "Convert {wildcards.experiment} {wildcards.accession} alignfile SAM to sorted BAM."
    shell:
        "samtools view -bS -o {params.tmp_bam} {input} && "
        "samtools sort -o {output} {params.tmp_bam} && "
        "samtools index {output}"


rule create_bedgraphs:
    input:
        rules.sam_to_sorted_bam.output
    output:
        plus_bg = "results/crssant/{experiment}/{accession}/{accession}_plus.bedgraph",
        minus_bg = "results/crssant/{experiment}/{accession}/{accession}_minus.bedgraph"
    params:
        "results/STAR/index/chrNameLength.txt"
    conda:
        "../envs/mapping.yaml"
    message:
        "Create {wildcards.experiment} {wildcards.accession} strand-separated bedgraphs."
    shell:
        "bedtools genomecov -bg -split -strand + -ibam {input} -g {params} | sed -n \'/^[0-9,X,Y,MT]/Ip\' | awk \'$1 !~ /_/\' | sort -k1,1 -k2,2n > {output.plus_bg} && "
        "bedtools genomecov -bg -split -strand - -ibam {input} -g {params} | sed -n \'/^[0-9,X,Y,MT]/Ip\' | awk \'$1 !~ /_/\' | sort -k1,1 -k2,2n > {output.minus_bg}"


rule DG_NG_assembly:
    input:
        alignfile = rules.combine_gap1filter_trans.output,
        genesfile = rules.create_genesfile.output,
        plus_bg = rules.create_bedgraphs.output.plus_bg,
        minus_bg = rules.create_bedgraphs.output.minus_bg
    output:
        bedpe = "results/crssant/{experiment}/{accession}/{accession}_pri_crssant_dg.bedpe"
    params:
        outdir = "results/crssant/{experiment}/{accession}/",
        cluster = "cliques",
        t_o = 0.1
    threads:
        16
    conda:
        "../envs/crssant.yaml"
    message:
        "Assemble {wildcards.experiment} {wildcards.accession} alignments to DGs and NGs."
    shell:
        "python3 workflow/scripts/crssant.py -out {params.outdir} "
        "-cluster {params.cluster} -n {threads} -t_o {params.t_o} "
        "{input.alignfile} {input.genesfile} {input.plus_bg},{input.minus_bg}"