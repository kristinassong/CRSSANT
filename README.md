# 🥐 NOTE 🥐
This is a modified version of the CRSSANT pipeline implemented using the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/) 🐍. 

Reference: https://github.com/zhipenglu/CRSSANT; https://github.com/whl-usc/rna2d3d

## Installing Snakemake
Please follow the [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install the Snakemake workflow management tool. We recommend using `Conda/Mamba` to install Snakemake.

This Snakemake workflow has been tested with `v7.32.4`.

## Downloading FASTQ files from GEO
This Snakemake workflow includes the following datasets. The datasets/samples to be included in the workflow can be modified in `config/config.yaml`.

**PARIS**
| Accession  | Sample Name   |
| -----------| ------------- |
| SRR2814761 | HeLa Cells    |
| SRR2814762 | HeLa Cells    |
| SRR2814763 | HEK293T Cells |
| SRR2814764 | HEK293T Cells |
| SRR2814765 | HEK293T Cells |

**PARIS2**
| Accession   | Sample Name   |
| ----------- | ------------- |
| SRR11624581 | HEK293T Cells |
| SRR11624582 | HEK293T Cells |
| SRR11624583 | HEK293T Cells |
| SRR11624584 | HEK293T Cells |
| SRR11624585 | HEK293T Cells |
| SRR11624586 | HEK293T Cells |
| SRR11624587 | HEK293T Cells |
| SRR11624588 | HEK293T Cells |

**SPLASH**
| Accession  | Sample Name          |
| ---------- | -------------------- |
| SRR3404924 | Lymphoblastoid Cells |
| SRR3404925 | Lymphoblastoid Cells |
| SRR3404936 | Lymphoblastoid Cells |
| SRR3404937 | Lymphoblastoid Cells |

**LIGR-seq**
| Accession  | Sample Name   |
| ---------- | ------------- |
| SRR3361013 | HEK293T Cells |
| SRR3361017 | HEK293T Cells |

Note that multiple runs of the Snakemake workflow may be required to successfully download all these datasets. We recommend verifying the download of all required datasets before moving onto the next steps of the CRSSANT pipeline.

## Running the Snakemake workflow
For a dry-run of this Snakemake workflow, simply run the following code from `CRSSANT/`.
```
snakemake -n
```
To run this Snakemake workflow, simply run the following code from `CRSSANT/`.
```
snakemake --profile profile_slurm
```

Please see below for more details on the CRSSANT pipeline. If you have any specific questions regarding this Snakemake workflow, please contact [Kristina Sungeun Song](mailto:kristina.song@usherbrooke.ca). Questions on the technicalities of CRSSANT should be addressed to the original authors: [Zhang et al. 2022 Genome Research](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9104705/).

# ⇩⇩⇩ Written by the original authors ⇩⇩⇩

# CRSSANT: Cross-linked RNA Secondary Structure Analysis using Network Techniques
RNA crosslinking, proximity ligation and high throughput sequencing produces non-continuous reads that indicate base pairing and higher order interactions, either in RNA secondary structures or intermolecular complexes. CRSSANT (pronounced 'croissant') is a computational pipeline for analyzing non-continuous/gapped reads from a variety of methods that employ the crosslink-ligation principle, including [PARIS](https://www.ncbi.nlm.nih.gov/pubmed/27180905), [LIGR](https://www.ncbi.nlm.nih.gov/pubmed/27184080), [SPLASH](https://www.ncbi.nlm.nih.gov/pubmed/27184079), [COMRADES](https://www.ncbi.nlm.nih.gov/pubmed/30202058), [hiCLIP](https://www.ncbi.nlm.nih.gov/pubmed/25799984), etc. CRSSANT optimizes short-read mapping, automates alignment processing, and clusters gap1 and trans alignments into duplex groups (DG) and non-overlapping groups (NG). More complex arrangments are assembled into higher level structures. In particular gapm alignments with 2 gaps or 3 segments are assembled into tri-segment groups (TGs). Overlapping alignments are used to discover homotypic interactions (RNA homodimers). 

Briefly, the CRSSANT pipeline operates as follows. First, sequencing reads that have been processed to remove adapters are mapped references with STAR and a new set of optimized options. Second, alignments are filtered, rearranged and classified into different types (gaptypes.py and gapfilter.py). Third, we use network analysis methods to cluster non-continuous alignments into DGs and calculate the confidence for each DG. The DGs are used as the foundation for the assembly of TGs. 

CRSSANT is written in Python and available as source code that you can download and run directly on your own machine (no compiling needed). An earlier version of the DG assembly method is available here: (https://github.com/ihwang/CRSSANT). For more about the CRSSANT pipeline, please refer to this study [Zhang et al. 2022 Genome Research](https://genome.cshlp.org/content/early/2022/03/24/gr.275979.121.abstract).

## Table of contents
* [Step 0: Download and prepare environment](https://github.com/zhipenglu/CRSSANT#step-0-download-and-prepare-environment)
* [Step 1: Preprocessing fastq input files](https://github.com/zhipenglu/CRSSANT#step-1-preprocessing-fastq-input-files)
* [Step 2: Map reads to the genome](https://github.com/zhipenglu/CRSSANT#step-2-map-reads-to-the-genome)
* [Step 3. Rearrange softclipped alignments and remap](https://github.com/zhipenglu/CRSSANT#step-3-rearrange-softclipped-alignments-and-remap)
* [Step 4: Classify alignments](https://github.com/zhipenglu/CRSSANT#step-4-classify-alignments)
* [Step 5: Segment and gap statistics](https://github.com/zhipenglu/CRSSANT#step-5-segment-and-gap-statistics)
* [Step 6: Filter spliced and short gaps](https://github.com/zhipenglu/CRSSANT#step-6-filter-spliced-and-short-gaps)
* [Step 7: Cluster gap1 and trans alignments to DGs](https://github.com/zhipenglu/CRSSANT#step-7-cluster-gap1-and-trans-alignments-to-dgs)
* [Step 8: Cluster gapm alignments to TGs](https://github.com/zhipenglu/CRSSANT#step-8-cluster-gapm-alignments-to-tgs)
* [Step 9: Analysis of RNA homodimers](https://github.com/zhipenglu/CRSSANT#step-9-analysis-of-rna-homodimers)
* [Running CRSSANT as a pipeline](https://github.com/zhipenglu/CRSSANT#running-crssant-as-a-pipeline)

## Step 0: Download and prepare environment
Download the scripts and save it to a known path/location. No special installation is needed, but the python package dependencies need to be properly resolved before use. You will need Python version 3.6+ and the following Python packages. We recommend downloading the latest versions of these packages using the Ananconda/Bioconda package manager. Currently, the NetworkX version only works with python 3.6, but not higher versions.

* [NetworkX v2.1+](https://networkx.github.io/) ([Anaconda link](https://anaconda.org/anaconda/networkx))
* [NumPy](http://www.numpy.org/) ([Anaconda link](https://anaconda.org/anaconda/numpy))
* [SciPy](https://www.scipy.org/) ([Anaconda link](https://anaconda.org/anaconda/scipy))
* [scikit-learn](http://scikit-learn.org/) ([Anaconda link](https://anaconda.org/anaconda/scikit-learn))

Additional tools for used for mapping and general processing of high throughput sequencing data, including STAR, samtools and bedtools. The STAR should be used with optimized parameters shown below. 
* [STAR v2.7.1+](https://github.com/alexdobin/STAR)
* [samtools v1.1+](http://www.htslib.org/)
* [bedtools v2.22+](https://bedtools.readthedocs.io/en/latest/)

For visualization of the results, we recommend IGV, which has features for grouping alignments based on tags, such as DG and NG that we implemented here. IGV can also directly visualize DG summary information and RNA secondary structures, see [Step 4: Cluster alignments to groups](https://github.com/zhipenglu/CRSSANT#step-4-cluster-alignments-to-groups) for details. VARNA is recommended for visualizing RNA secondary structures in a variety of formats, including 
* [IGV v2.4+](https://software.broadinstitute.org/software/igv/)
* [VARNA v1.0+](http://varna.lri.fr/)


### System requirements and tests

The programs are generally run in x86-64 compatible processors, including 64 bit Linux or Mac OS X, if there is enough memory. Read mapping against mammalian genomes using STAR requires at least 30G memory. Alignment classification typically requires 100G memory. As a result, these two steps should be run in a cluster with large memory. 

Test datasets and example output files are provided for all steps except STAR mapping, which is a well maintained package. Test files are located in the `tests` folder. Furthermore, we provided source data for all figures in this paper to help readers reproduce the figures and troubleshoot potential problems (`sourcedata`). The analysis pipeline is preferably run as separate steps to allow maximal control and quality assurance. In addition, shell scripts for a typical pipeline are also provided as examples. 


## Step 1: Preprocessing fastq input files
Sequencing data can be processed using various published tools to demultiplex samples, and remove barcodes. Since each library preparation uses a different approach, we cannot recommend the same method for all of them. Here is an example of preprocessing PARIS data based on our recently published library preparation protocol. 


## Step 2: Map reads to the genome
It is assumed that the reads have been demultiplexed and adapters removed. Before mapping the reads, genome indices should be generated with the same STAR version. Reads in the fastq format are mapped to the genome using STAR and a set of optimized parameters as follows. `runThreadN` and `genomeLoad` should be adjusted based on available resources and running environment.  

```
STAR --runMode alignReads --genomeDir /path/to/index --readFilesIn /path/to/reads/files --outFileNamePrefix /path/to/output/prefix --runThreadN 1 --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  --outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outSAMattributes All --outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreGenomicLengthLog2scale -1 --chimFilter None --chimOutType WithinBAM HardClip --chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 -- chimScoreDropMax 80 --chimNonchimScoreDropMin 20
```
Successful STAR mapping generates the following 7 files: `Aligned.out.bam`, `Aligned.sortedByCoord.out.bam`, `Log.final.out`, `Log.out`, `Log.progress.out`, `SJ.out.tab`, and `Unmapped.out.mate1`. The `bam` file is converted back to `sam` for the next step of processing, keeping the header lines (`samtools view -h`). Sorting is not necessary for the next alignment classification step.  

### --Optimized STAR parameters
Here is a brief explanation of the optimized parameters for non-continuous alignments. See the bioRxiv preprint referenced at the top of this README for more detailed discussion of the optimization. The output contains very short segments, down to 7nt, and the unreliable ones are removed later using alignment-span-based penalty and segment connection dependent filtering (`gaptypes.py`).

* `--outFilterScoreMinOverLread 0` and `--outFilterMatchNminOverLread 0` allows mapping short segments
* `--outSAMattributes All` includes chimeric tags needed for alignment processing
* `--outSAMtype BAM Unsorted SortedByCoordinate` simplifies subsequent SAM processing
* `--alignIntronMin 1` shifts deletions (`D`) to gaps (`N`) to equalize penalty. 
* `--scoreGap* 0` removes all gap open penalty.
* `--scoreGenomicLengthLog2scale -1` increases alignment span-based penalty
* `--chimFilter None` enables detection of chimeric alignments (primarily homotypic) near the 5' and 3' ends of references
* `--chimOutType WithinBAM HardClip` output alignments in one file and removes hardclips. 
* `--chimSegmentMin 5` and `--chimJunctionOverhangMin 5` map chimera more permissively
* `--chimScoreJunctionNonGTAG 0` removes penalty for splicing junctions in chimera
* `-- chimScoreDropMax 80`, a higher value, and `--chimNonchimScoreDropMin 20` ensures that chimera are not produced when normal gapped alignments are possible. 

### --Other relevant parameters
In addition to the parameters listed above, the following ones need to be adjusted according to the data. Larger datasets may produce more gapped alignments (similar to splice junctions) and therefore requires higher limits. Even higher numbers can be tried if the following recommended numbers are still insufficient. 
* `--limitOutSJcollapsed`: recommended 10,000,000 (default 1 million)
* `--limitIObufferSize`: recommended 1,500,000,000 (default 150 million)


## Step 3: Rearrange softclipped alignments and remap
The STAR mapper, even with the optimized parameters discrimates against backward chimeric alignments compared to forward chimera (on the same strand and chromosome). To increase the detection of backward chimera, softclipped continuous alignments from the first round STAR mapping are rearranged so that they can be remapped by STAR in a second round. To rearrange softclipped continuous alignments, use the script softreverse.py as follows. The softrev.fastq file is mapped again using the parameters listed in Step 2. The output from the two rounds of STAR mappings are combined for subsequent analysis (Step 4). 
```
python softreverse.py input.sam softrev.fastq
```


## Step 4: Classify alignments
In this step, alignments in the sam file are filtered to remove low-confidence segments, rearranged and classified into ![5 distinct types](figures/s3_noncon.pdf) using `gaptypes.py`. 
```
python gaptypes.py input.sam output_prefix glenlog nprocs
```
Recommended parameters are as follows. 
* `glenlog`: -1. Scaling factor for gap extension penalty, equivalent to `scoreGenomicLengthLog2scale` in STAR 
* `minlen`: 15. Minimal length for a segment in an alignment to be considered confident for building the connection database
* `nprocs`: 10. Number of CPUs to use for the run, depending availability of resources. 


### --Output from `gaptypes.py`
Successful completion of this step results in 7 files. All of these sam files can be converted to sorted bam for visualization on IGV. 
* `cont.sam`: continuous alignments
* `gap1.sam`: non-continuous alignments, each has 1 gap
* `gapm.sam`: non-continuous alignments, each has more than 1 gaps
* `trans.sam`: non-continuous alignments with the 2 arms on different strands or chromosomes
* `homo.sam`: non-continuous alignments with the 2 arms overlapping each other
* `bad.sam`: non-continuous alignments with complex combinations of indels and gaps 
* `log.out`: log file for the run, including input and output alignment counts (for `gap1`, `gapm`, `trans`, `homo` and `bad`)


### --Testing `gaptypes.py`
Here is an example test of the `gaptypes.py` script for filtering, rearranging and classifying alignments, using data provided in `tests/gaptypes`. The input data are alignments to two RNAs, SNORD118 (U8) chr17:8173453-8173588 and SNORD113 (U13) chr8:33513475-33513577, from HEK293 PARIS reads mapped to hg38pri. The output sam files can be converted to sorted bam for visualization on IGV under hg38 genome reference. 
```
python gaptypes.py PARIS_U8U13_hg38pri_only.sam PARIS_U8U13_hg38pri_only -1 15 8
Output counts: 
Total input alignment number: 10093
Continuous alignments (no gaps, cont.sam): 9179
Two-segment gapped alignments (gap1.sam): 339
Multi-segment gapped alignments (gapm.sam): 1
Other chimeric (different str/chr, trans.sam): 1
Overlapping chimeric (homotypic, homo.sam): 281
Bad alignments (bad.sam): 0
```

## Step 5: Segment and gap statistics
Segment and gap length distribution can be produced using two scripts `seglendist.py` and `gaplendist.py` to help understand the quality and properties of the sequencing data. Both scripts either use the sam file to produce a list of numbers, or use a list of numbers from a file to produce the cumulative distribution histogram. This step is not required for the subsequent steps; it only provides information to assess data quality. 

Example command for creating the gap length list is as follows, where input.sam can be gap1.sam and gapm.sam, or their filtered output where splicing junctions and/or short gaps have been removed. 
```
python gaplendist.py inputfile filetype outputfile gaptype
python seglendist.py inputfile filetype outputfile
```
* `filetype`: 'sam' file or text file containing 'list' of lengths 
* `gaptype`: 'min', only the shortest gap in the alignment, or 'all' for all gaps

When input is sam, output is list file. When input is a file of lengths list, output is pdf figure. The percentage of gaps or segments within a certain range can are output as well when running these scripts. The figure output can be adjusted by changing parameters in the python script. 


## Step 6: Filter spliced and short gaps
Output files `gap1.sam` and `gapm.sam` may contain alignments that have only splicing junctions and short 1-2 nt gaps due to artifacts. These are filtered out using `gapfilter.py` before further processing. Splicing junctions and short gaps in other output files can be safely ignored. The `annotation` file containing the splicing junctions should be in GTF format. `idloc`, location of the transcript_ID field, is usually field 11. `short` is set to either `yes` which means 'remove short 1-2nt gaps', or `no`, which means 'ignore short 1-2nt gaps'.  
```
python gapfilter.py annotation insam outsam idloc short
```
The output from `gap1.sam`, typically named `gap1filter.sam`, only contains alignments that pass the filter. The output from `gapm.sam`, typically named gapmfilter.sam, contain alignments with either 1 or more gaps that pass the filter. The following output is printed to the screen: 

* `Total alignments`
* `All gapped alignments`
* `Alignments with at least 1 good gaps`
* `Alignments with at least 2 good gaps`
* `Number of annotated splicing junctions`

Alignments counts in log.out from step 2 and filtered counts are used to calculate percentage of non-continuous alignments alignments using the following formula: 
```
(gap1filtercount + gapmfiltercount + transcount + homocount)/inputcount
```
### --Testing `gapfilter.py`
Here is an example test of the gapfilter.py script using data provided in `tests/gapfilter` and the output stats. The input and output sam files can be converted to sorted bam for visualization on IGV under hg38 genome reference. 

```
python gapfilter.py ACTB.gtf PARIS_hg38pri_ACTB_gap1.sam PARIS_hg38pri_ACTB_gap1filter.sam 11 yes
Output: 
Total alignments: 352
All gapped alignments: 352
Alignments with at least 1 good gaps: 291
Alignments with at least 2 good gaps: 0
Number of annotated splicing junctions: 5
```

## Step 7: Cluster gap1 and trans alignments to DGs
After filtering alignments, To assemble alignments to DGs and NGs using the crssant.py script, three types of input files are required, `alignfile`, `genesfile` and `bedgraphs`. For more on these parameters, see the explanation below and the bioRxiv preprint referenced at the top of this README. 
```
python crssant.py [-h] [-out OUT] [-cluster CLUSTER] [-n N] [-covlimit COVLIMIT] [-t_o T_O] [-t_eig T_EIG] alignfile genesfile bedgraphs
```

Positional arguments:
* `alignfile`: Path to alignedment file (SAM)
* `genesfile`: Path to gene annotations (BED)
* `bedgraphs`: Path to genome coverage files. Coma separated 2 files for + and - strands

Optional arguments: 
* `-h`: show the help message and exit
* `-out`: path of output folder
* `-cluster`: clustering method, "cliques" or "spectral"(default)
* `-n`: Number of threads. Default is 8
* `-covlimit`: Max coverage to be directly graphed. Default is 1000. 
* `-t_o`: Overlap threshold 0-1 inclusive. Default: 0.5 for "spectral", 0.1 for "cliques"
* `-t_eig`: Eigenratio threshold (positive) for "spectral" only. Default: 5. 


### --Required input files
Alignment files `gap1filter.sam` and `trans.sam` are combined to produce `alignfile` sam, keeping one set of header lines at the beginning. The header lines are passed on to output. `genesfile` is a list of all genes in the genome. The start and end for each gene is used to assign the alignments to gene pairs, and determine which alignments correspond to intramolecular structures or intermolecular interactions. The `bedgraphs` files can be produced using the bedtools package, for example: 
```
bedtools genomecov -bg -split -strand + -ibam alignfile_sorted.bam -g chromosome_sizes > alignfile_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam alignfile_sorted.bam -g chromosome_sizes > alignfile_minus.bedgraph
```

### --Clustering methods
The pipeline uses the spectral clustering method to cluster reads into DGs with overlap threshold parameter of `t_o=0.5` and eigenratio threshold of `t_eig=5`. The default spectral clustering method may be operated with different overlap threshold and eigenratio threshold parameters by specifying one or both with the flags `t_o` and `t_eig`, respectively. `t_o` may be any float between 0 and 1, and `t_eig` may be any positive number. Increasing `t_o` tends to result in more DGs that contain fewer reads, and increasing `t_eig` tends to result in fewer DGs containing more reads.

The user may also specify the cliques-finding method for clustering DGs by specifying the clustering flag `cluster` with the `cliques` option. If the cliques-finding method is specified, `t_o` may also be specified, and again may be any float between 0 and 1. The eigenvalue threshold `t_eig` is not used in the cliques-finding method. By default, for the cliques-finding method the overlap threshold is set to 0.1. 

### --Output from `crssant.py`
After DG clustering, crssant.py verifies that the DGs do not contain any non-overlapping reads, i.e. any reads where the start position of its left arm is greater than or equal to the stop position of the right arm of any other read in the DG. If the DGs do not contain any non-overlapping reads, then the following output files ending in the following are written:

* `.sam`: SAM file containing alignments that were successfully assigned to DGs, plus DG and NG annotations
* `dg.bedpe`: bedpe file listing all duplex groups.

The sam output can be converted to bam for visualization in IGV, where DG and NG tags can be used to sort and group alignments. The bed output file can be visualized in IGV, where the two arms of each DG can be visualized as two 'exons', or as an arc the connects far ends of the DG (http://software.broadinstitute.org/software/igv/node/284). 

DGs with the 2 arms on the same strand and chromosome can be further converted to bed12 for visualization on IGV. The fields of the BED12 file are used according to the [standard definition](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) with the exception of fields 4 and 5. Field 5 `score` (or field 8 in bedpe) is defined as the number of non-continuous alignments in this DG.  Field 4 `name` (or field 7 in bedpe) is defined in the format `gene1,gene2_DGID_covfrac`, where `gene1,gene2` represents the two genes that this DG connects, `DGID` is a numerical ID of the DG, `covfrac` is the confidence of the DG, defined as c / sqrt(a\*b) and
* c = number of reads in a given DG
* a = number of reads overlapping the left arm of the DG
* b = number of reads overlapping the right arm of the DG


### --Testing `crssant.py`
Here is an example test of the `crssant.py` script for DG and NG assembly, using data provided in `tests/crssant`. The output sam files can be converted to sorted bam for visualization on IGV under hg38 genome reference. The DG and NG tags can be used to group and sort alignments. 
```
python crssant.py -out ./ -cluster cliques ACTB.sam ACTB.bed ACTB_plus.bedgraph,ACTB_minus.bedgraph
```

DGs on the same chromosome and strand can be converted to bed12 format for visualization on IGV. First use `bedpetobed12.py` to convert, then use bedtools to sort bed file and finally add the header line `track graphType=arc`. 
```
python bedpetobed12.py ACTB.cliques.t_o0.1_dg.bedpe ACTB.cliques.t_o0.1_dg.bed
sortBed -i ACTB.cliques.t_o0.1_dg.bed > ACTB.cliques.t_o0.1_dg_sorted.bed
```


## Step 8: Cluster gapm alignments to TGs
After assembly of DGs from single-gap or two-segment alignments (including gap1 and trans), the DGs are used as the foundation to assemble tri-segment groups (TGs) from gapm alignments (only 3-segment ones are assembled at the moment, since reads with more than 3 segments are extremely rare). 

### --Required input files
The output bedpe file from the DG assembly for the RNA of interest is used as the foundation to assemble TGs. The gapmfilter.sam file which contains the filtered gapm alignments are assembled into the TGs. 

### --Output from `gapmcluster.py`
TG clustering produces a sam file just like the gapm input sam file, except the addition of a new TG tag at the end of each alignment record. The output sam files can be converted to sorted bam for visualization on IGV under hg38 genome reference. The TG tag can be used to group and sort alignments. 

### --Testing `gapmcluster.py`
Here is an example test of the `gapmcluster.py` script for TG assembly, using data provided in `tests/gapm`. This example is the human 7SK RNA. The correlation between DG and TG coverages can be calculated using the DGTGcorr.py script, for example, using data from the Fig. 5E sourcedata folder. 
```
python gapmcluster.py RN7SK_hg38_manualDGs.bedpe RN7SK_hg38_gapm.sam
```



## Step 9: Analysis of RNA homodimers
The overlapping chimeric alignments indicate homotypic interactions, or RNA homodimers. The homo.sam file is further processed as follows to identify potential homodimers. 
To ensure the identification of RNA homodimers using STAR mapping and gaptypes.py classification, the RNAs of interest must be flanked by additional non-N sequences, or the `--chimFilter None` option needs to be set. Homo alignments (homo.sam) with less than 2nt overlapping between two arms were filtered out to avoid potential artifacts. To cluster overlapping alignments into groups with similar DG/NG tags as the gap1/trans alignments, the crssant.py script is applied as follows. The input and output files are similar to the descriptions above. 

```
python crssant.py [-h] [-out OUT] [-cluster CLUSTER] [-n N] [-covlimit COVLIMIT] [-t_o T_O] [-t_eig T_EIG] alignfile genesfile bedgraphs
```

To plot the homodimers as heatmaps, use the following script:

```
python plot_heatmap_for_homodimer.py
```


## Running CRSSANT as a pipeline
The entire pipeline can be run step by step to give the user maximal control. We prefer this approach for two reasons. First, it is often easier for beginners to locate the problems and fix them. Second, even for experienced users, some of the steps may require changes of parameters, when running the pipeline step by step allows the users to only rerun the necessary steps, and therefore saves time. Once all trouble shooting and parameter testing are completed, a python script is used to integrate the pipeline, steps 4-8, after the STAR mapping step. 

```
python run_CRSSANT.py input_dir output_dir sample_name Gtf idloc genesfile outprefix
```
* `input_dir`: Full path to the directory of input files
* `output_dir`: Full path to the directory of output files
* `sample_name`: Name of the sample
* `Gtf`: The annotation file containing the splicing junctions
* `idloc`: Column number of the transript_ID field in GTF files, typically 11.
* `genesfile`: Gene annotations (BED format)
* `outprefix`: Output prefix

### --Input and output files 
Before running `run_CRSSANT.py`, ensure that the following files are placed in a single location/folder. 
* BAM/SAM file containing aligned reads from the sample of interest using the parameters listed above. 
* Gene bed file should contain six required BED fileds: Chrom, ChromStart, ChromEnd, Name, Score, Strand.
* GTF file, the annotation file containing the splicing junctions should be in GTF format.

The output files from run_CRSSANT.py are gathered in four folders:
* alignments_classify: classified alignments (6 categories)
* alignments_statistics: alignment statistics, containing segment and gap length distributions
* alignments_DGs: created at step of DGs assembly
* alignments_TGs: created at step of TGs assembly

### --Output folder 1: alignments_classify
Successful completion of this step results in 7 files. All of these sam files can be converted to sorted bam for visualization on IGV.
* cont.sam: continuous alignments
* gap1.sam: non-continuous alignments with 1 gap
* gapm.sam: non-continuous alignments with more than 1 gaps
* trans.sam: non-continuous alignments with the 2 arms on different strands or chromosomes
* homo.sam: non-continuous alignments with the 2 arms overlapping each other
* bad.sam: non-continuous alignments with complex combinations of indels and gaps
* log.out: log file for the run, including input and output alignment counts (for cont, gap1, gapm, trans, homo and bad)

### --Output folder 2: alignments_statistics
The following two files are generated to show the distribution of gap and segment lengths. 
* GapLen_distribution.pdf: Gap length distribution
* SegLen_distribution.pdf: Segment length distribution

### --Output folder 3: alignments_DGs
The following 2 files are produced for DG assembly. 
* .sam: SAM file containing alignments that were successfully assigned to DGs, plus DG and NG annotations
* dg.bedpe: bedpe file listing all duplex groups.

### --Output folder 4: alignments_TGs
TG clustering produces a sam file just like the gapm input sam file, except the addition of a new TG tag at the end of each alignment record. The output sam files can be converted to sorted bam for visualization on IGV. The TG tag can be used to group and sort alignments.


## END
