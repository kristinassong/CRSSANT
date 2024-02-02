#!/usr/bin/env python


import pandas as pd
import pyranges as pr


def extract_gtf(gtf_df,info_type):
    """
    Extract gene information from annotation gtf.
    """

    genes_df = gtf_df

    if info_type == 'gene':
        genes_df = gtf_df[gtf_df['Feature']=='gene'].reset_index(drop=True)

    if info_type == 'CdsUtr':
        genes_df = gtf_df[(gtf_df['Feature']=='CDS')|(gtf_df['Feature']=='five_prime_utr')|(gtf_df['Feature']=='three_prime_utr')].reset_index(drop=True)
    
    genes_df['Score'] = 1000
    genes_df = genes_df[['Chromosome','Start','End','gene_id','Score','Strand']]

    return genes_df


def main():
    gtf = snakemake.input[0]
    outfile = snakemake.output[0]
    info_type = snakemake.params.info_type

    gtf_df = pr.read_gtf(gtf).df
    extract_gtf(gtf_df,info_type).to_csv(outfile,sep='\t',index=False,header=False)


if __name__ == '__main__':
    main()