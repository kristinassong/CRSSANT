###############################################################
# Arm length stats
###############################################################

rule merge_gap1_gapm_trans:
    input:
        gap1_filtered = rules.filter_spliced_short_gaps_gap1.output,
        gapm_filtered = rules.filter_spliced_short_gaps_gapm.output,
        trans = rules.merge_sam.output.trans
    output:
        "results/stats/{experiment}/{accession}/{accession}_pri_filtered_arm.sam"
    message:
        "Merge {wildcards.experiment} {wildcards.accession} gap1_filtered.sam, gapm_filtered.sam and trans.sam to generate arm length stats."
    shell:
        "cat {input.gap1_filtered} {input.gapm_filtered} {input.trans} > {output}"


rule segment_stats_I:
    input:
        rules.merge_gap1_gapm_trans.output
    output:
        "results/stats/{experiment}/{accession}/{accession}_seglen.list"
    params:
        file_type = "sam"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} arm length distribution - Part I (Produce a list of numbers)."
    script:
        "../scripts/seglendist.py"


rule segment_stats_II:
    input:
        rules.segment_stats_I.output
    output:
        "results/stats/{experiment}/{accession}/{accession}_seglen.pdf"
    params:
        file_type = "list"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} arm length distribution - Part II (Produce a cumulative distribution histogram)."
    script:
        "../scripts/seglendist.py"


###############################################################
# Gap length stats
###############################################################

rule merge_gap1_gapm:
    input:
        gap1_filtered = rules.filter_spliced_short_gaps_gap1.output,
        gapm_filtered = rules.filter_spliced_short_gaps_gapm.output
    output:
        "results/stats/{experiment}/{accession}/{accession}_pri_filtered_gaps.sam"
    message:
        "Merge {wildcards.experiment} {wildcards.accession} gap1_filtered.sam and gapm_filtered.sam to generate gap length stats."
    shell:
        "cat {input.gap1_filtered} {input.gapm_filtered} > {output}"


rule gap_stats_I:
    input:
        rules.merge_gap1_gapm.output
    output:
        "results/stats/{experiment}/{accession}/{accession}_gaplen.list"
    params:
        file_type = "sam",
        gap = "all"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} gap length distribution - Part I (Produce a list of numbers)."
    script:
        "../scripts/gaplendist.py"


rule gap_stats_II:
    input:
        rules.gap_stats_I.output
    output:
        "results/stats/{experiment}/{accession}/{accession}_gaplen.pdf"
    params:
        file_type = "list",
        gap = "all"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} gap length distribution - Part II (Produce a cumulative distribution histogram)."
    script:
        "../scripts/gaplendist.py"