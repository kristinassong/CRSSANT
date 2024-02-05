###############################################################
# CoCo -- Corrected annotation
###############################################################

rule coco_ca:
    input:
        gtf = config["genome_gtf"],
        coco_dir = rules.download_CoCo_git.output
    output:
        gtf_corrected = "resources/correct_annotation.gtf"
    conda:
        "../envs/coco.yaml"
    message:
        "Generate corrected annotation using CoCo."
    shell:
        "export PATH=$PWD/{input.coco_dir}/bin:$PATH && "
        "python3 {input.coco_dir}/bin/coco.py ca {input.gtf} -o {output.gtf_corrected}"


###############################################################
# Map reads to the genome -- ROUND1
###############################################################

rule STAR_index:
    input:
        fasta = config["genome_fasta"],
        gtf = rules.coco_ca.output.gtf_corrected
    output:
        directory("results/STAR/index")
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    log:
        "results/logs/STAR/index.log"
    message:
        "Generate genome index files using STAR."
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99 "
        "&> {log}"


rule STAR_align_1:
    input:
        fq = "resources/preprocessed_fastq/{experiment}/{accession}_preprocessed.fastq",
        idx = rules.STAR_index.output
    output:
        "results/STAR/{experiment}/{accession}/{accession}_1_Aligned.sortedByCoord.out.bam"
    params:
        out_prefix = "results/STAR/{experiment}/{accession}/{accession}_1_"
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    log:
        "results/logs/STAR/{experiment}/{accession}_1_align.log"
    message:
        "Align {wildcards.experiment} {wildcards.accession} reads to the reference genome using STAR (ROUND1)."
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.fq} "
        "--outFileNamePrefix {params.out_prefix} "
        "--runThreadN {threads} --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  "
        "--outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 "
        "--outFilterMatchNminOverLread 0 --outSAMattributes All "
        "--outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 "
        "--scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 "
        "--scoreGenomicLengthLog2scale -1 --chimFilter None --chimOutType WithinBAM HardClip "
        "--chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 "
        "-- chimScoreDropMax 80 --chimNonchimScoreDropMin 20 "
        "--limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 "
        "&> {log}"


rule primary_alignments_1:
    input:
        rules.STAR_align_1.output
    output:
        "results/STAR/{experiment}/{accession}/{accession}_1_pri.sam"
    params:
        nonchimeric_temp_sam = "results/STAR/{experiment}/{accession}/{accession}_nonchimeric_temp.sam",
        chimeric_temp_sam = "results/STAR/{experiment}/{accession}/{accession}_chimeric_temp.sam",
        nonchimeric_pri_bam = "results/STAR/{experiment}/{accession}/{accession}_nonchimeric_pri.bam",
        nonchimeric_pri_sam = "results/STAR/{experiment}/{accession}/{accession}_nonchimeric_pri.sam"
    conda:
        "../envs/mapping.yaml"
    message:
        "Extract primary reads for {wildcards.experiment} {wildcards.accession} (ROUND1)."
    shell:
        "samtools view -h {input} | awk \'$1~/^@/ || NF<21\' > {params.nonchimeric_temp_sam} && "
        "samtools view -h {input} | awk \'$1!~/^@/ && NF==21\' > {params.chimeric_temp_sam} && "
        "samtools view -bS -F 0x900 -o {params.nonchimeric_pri_bam} {params.nonchimeric_temp_sam} && "
        "samtools view -h {params.nonchimeric_pri_bam} > {params.nonchimeric_pri_sam} && "
        "cat {params.nonchimeric_pri_sam} {params.chimeric_temp_sam} > {output} && "
        "rm -f {params.nonchimeric_temp_sam} {params.chimeric_temp_sam} {params.nonchimeric_pri_bam} {params.nonchimeric_pri_sam}"


rule classify_alignments_1:
    input:
        rules.primary_alignments_1.output
    output:
        gap1 = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_gap1.sam",
        gapm = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_gapm.sam",
        trans = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_trans.sam",
        homo = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_homo.sam",
        cont = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_cont.sam"
    params:
        output_prefix = "results/gaptypes/{experiment}/{accession}/{accession}_1_pri_",
        glenlog = -1,
        minlen = 15,
        nprocs = 10
    conda:
        "../envs/crssant.yaml"
    message:
        "Filter {wildcards.experiment} {wildcards.accession} alignments to remove low-confidence segments, rearrange and classify into 5 distinct types (ROUND1)."
    script:
        "../scripts/gaptypes.py"


###############################################################
# Rearrange softclipped continuous reads and remap -- ROUND2
###############################################################

rule softreverse:
    input:
        rules.classify_alignments_1.output.cont
    output:
        "resources/preprocessed_fastq/{experiment}/{accession}_preprocessed.softrev.fastq"
    message:
        "Rearrange {wildcards.experiment} {wildcards.accession} softclipped continuous alignments."
    script:
        "../scripts/softreverse.py"


rule STAR_align_2:
    input:
        fq = rules.softreverse.output,
        idx = rules.STAR_index.output
    output:
        "results/STAR/{experiment}/{accession}/{accession}_2_Aligned.sortedByCoord.out.bam"
    params:
        out_prefix = "results/STAR/{experiment}/{accession}/{accession}_2_"
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    log:
        "results/logs/STAR/{experiment}/{accession}_2_align.log"
    message:
        "Align rearranged {wildcards.experiment} {wildcards.accession} reads to the reference genome using STAR (ROUND2)."
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.fq} "
        "--outFileNamePrefix {params.out_prefix} "
        "--runThreadN {threads} --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  "
        "--outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 "
        "--outFilterMatchNminOverLread 0 --outSAMattributes All "
        "--outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 "
        "--scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 "
        "--scoreGenomicLengthLog2scale -1 --chimFilter None --chimOutType WithinBAM HardClip "
        "--chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 "
        "-- chimScoreDropMax 80 --chimNonchimScoreDropMin 20 "
        "--limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 "
        "&> {log}"


rule primary_alignments_2:
    input:
        rules.STAR_align_2.output
    output:
        "results/STAR/{experiment}/{accession}/{accession}_2_pri.sam"
    params:
        nonchimeric_temp_sam = "results/STAR/{experiment}/{accession}/{accession}_2_nonchimeric_temp.sam",
        chimeric_temp_sam = "results/STAR/{experiment}/{accession}/{accession}_2_chimeric_temp.sam",
        nonchimeric_pri_bam = "results/STAR/{experiment}/{accession}/{accession}_2_nonchimeric_pri.bam",
        nonchimeric_pri_sam = "results/STAR/{experiment}/{accession}/{accession}_2_nonchimeric_pri.sam"
    conda:
        "../envs/mapping.yaml"
    message:
        "Extract primary reads for rearranged {wildcards.experiment} {wildcards.accession} (ROUND2)."
    shell:
        "samtools view -h {input} | awk \'$1~/^@/ || NF<21\' > {params.nonchimeric_temp_sam} && "
        "samtools view -h {input} | awk \'$1!~/^@/ && NF==21\' > {params.chimeric_temp_sam} && "
        "samtools view -bS -F 0x900 -o {params.nonchimeric_pri_bam} {params.nonchimeric_temp_sam} && "
        "samtools view -h {params.nonchimeric_pri_bam} > {params.nonchimeric_pri_sam} && "
        "cat {params.nonchimeric_pri_sam} {params.chimeric_temp_sam} > {output} && "
        "rm -f {params.nonchimeric_temp_sam} {params.chimeric_temp_sam} {params.nonchimeric_pri_bam} {params.nonchimeric_pri_sam}"


rule classify_alignments_2:
    input:
        rules.primary_alignments_2.output
    output:
        gap1 = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_gap1.sam",
        gapm = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_gapm.sam",
        trans = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_trans.sam",
        homo = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_homo.sam",
        cont = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_cont.sam"
    params:
        output_prefix = "results/gaptypes/{experiment}/{accession}/{accession}_2_pri_",
        glenlog = -1,
        minlen = 15,
        nprocs = 10
    conda:
        "../envs/crssant.yaml"
    message:
        "Filter rearranged {wildcards.experiment} {wildcards.accession} alignments to remove low-confidence segments, rearrange and classify into 5 distinct types (ROUND2)."
    script:
        "../scripts/gaptypes.py"


###############################################################
# Combine output from both rounds of STAR mapping
###############################################################

rule merge_sam:
    input:
        gap1_1 = rules.classify_alignments_1.output.gap1,
        gap1_2 = rules.classify_alignments_2.output.gap1,
        gapm_1 = rules.classify_alignments_1.output.gapm,
        gapm_2 = rules.classify_alignments_2.output.gapm,
        homo_1 = rules.classify_alignments_1.output.homo,
        homo_2 = rules.classify_alignments_2.output.homo,
        trans_1 = rules.classify_alignments_1.output.trans,
        trans_2 = rules.classify_alignments_2.output.trans
    output:
        gap1 = "results/gaptypes/{experiment}/{accession}/{accession}_pri_gap1.sam",
        gapm = "results/gaptypes/{experiment}/{accession}/{accession}_pri_gapm.sam",
        homo = "results/gaptypes/{experiment}/{accession}/{accession}_pri_homo.sam",
        trans = "results/gaptypes/{experiment}/{accession}/{accession}_pri_trans.sam",
    params:
        gap1_tmp = "results/gaptypes/{experiment}/{accession}/{accession}_pri_gap1.tmp",
        gapm_tmp = "results/gaptypes/{experiment}/{accession}/{accession}_pri_gapm.tmp",
        homo_tmp = "results/gaptypes/{experiment}/{accession}/{accession}_pri_homo.tmp",
        trans_tmp = "results/gaptypes/{experiment}/{accession}/{accession}_pri_trans.tmp",
        header = "results/gaptypes/{experiment}/{accession}/{accession}_header.txt",
        header_bam = rules.STAR_align_1.output
    conda:
        "../envs/mapping.yaml"
    message:
        "Combine outputs from both rounds of STAR mapping for {wildcards.experiment} {wildcards.accession}."
    shell:
        "python3 workflow/scripts/merger_sams.py {input.gap1_1} {input.gap1_2} none gap1 {params.gap1_tmp} && "
        "python3 workflow/scripts/merger_sams.py {input.gapm_1} {input.gapm_2} none gapm {params.gapm_tmp} && "
        "python3 workflow/scripts/merger_sams.py {input.homo_1} {input.homo_2} none homo {params.homo_tmp} && "
        "python3 workflow/scripts/merger_sams.py {input.trans_1} {input.trans_2} none trans {params.trans_tmp} && "
        "samtools view -H {params.header_bam} > {params.header} && "
        "cat {params.header} {params.gap1_tmp} > {output.gap1} && "
        "cat {params.header} {params.gapm_tmp} > {output.gapm} && "
        "cat {params.header} {params.homo_tmp} > {output.homo} && "
        "cat {params.header} {params.trans_tmp} > {output.trans} && "
        "rm -f {params.header} results/gaptypes/{wildcards.experiment}/{wildcards.accession}/*.tmp"


###############################################################
# Filter spliced and short gaps
###############################################################

rule filter_spliced_short_gaps_gap1:
    input:
        rules.merge_sam.output.gap1
    output:
        "results/gaptypes/{experiment}/{accession}/{accession}_pri_gap1_filtered.sam"
    params:
        annotation = rules.coco_ca.output.gtf_corrected,
        idloc = 11,
        short = "yes"
    message:
        "Filter {wildcards.experiment} {wildcards.accession} gap1 alignments that have only splicing junctions and short 1-2 nt gaps due to artifacts."
    script:
        "../scripts/gapfilter.py"


rule filter_spliced_short_gaps_gapm:
    input:
        rules.merge_sam.output.gapm
    output:
        "results/gaptypes/{experiment}/{accession}/{accession}_pri_gapm_filtered.sam"
    params:
        annotation = rules.coco_ca.output.gtf_corrected,
        idloc = 11,
        short = "yes"
    message:
        "Filter {wildcards.experiment} {wildcards.accession} gapm alignments that have only splicing junctions and short 1-2 nt gaps due to artifacts."
    script:
        "../scripts/gapfilter.py"