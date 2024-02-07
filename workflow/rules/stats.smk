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
        "results/stats/{experiment}/{accession}/{accession}_seglen.pdf"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} arm length distribution."
    shell:
        "python3 workflow/scripts/seglendist.py {input} {output} && "
        "rm -f {input}"


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
        "results/stats/{experiment}/{accession}/{accession}_gaplen.pdf"
    conda:
        "../envs/matplotlib.yaml"
    message:
        "Generate {wildcards.experiment} {wildcards.accession} gap length distribution."
    shell:
        "python3 workflow/scripts/gaplendist.py {input} {output} 100 && "
        "rm -f {input}"