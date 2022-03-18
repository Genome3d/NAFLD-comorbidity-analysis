#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
import argparse
from scipy.stats import hypergeom
import statsmodels.stats.multitest as mt
from tqdm import tqdm


def read_input_file(snps):
    print('Reading input file...')
    colnames = ['snp']
    f = pd.read_csv(snps, sep = "\t", header=None, usecols = [0], names = colnames).drop_duplicates()
    return f

def parse_gwas(gwas):
    print('Parsing GWAS associations...')
    cols = ['SNPS', 'SNP_ID_CURRENT', 'DISEASE/TRAIT', #'MAPPED_GENE',
            'P-VALUE',  'OR or BETA', '95% CI (TEXT)']
    df = pd.read_csv(gwas, sep='\t', usecols=cols, low_memory=False)
    #df['SNP_ID_CURRENT'] = df['SNP_ID_CURRENT'].fillna('').apply(lambda x: clean_snp_ids(x))
    df = df.assign(SNPS=df['SNPS'].str.split(';')).explode('SNPS')
    df = df.assign(SNPS=df['SNPS'].str.split(',')).explode('SNPS')
    df['SNPS'] = df['SNPS'].str.strip()
    return df

def clean_snp_ids(snp_id):
    try:
        return f"rs{int(snp_id)}"
    except:
        if snp_id == '':
            return snp_id
        else:
            snp_id = snp_id.split('-')[0]
            snp_id = snp_id.split('_')[0]
            return f"rs{snp_id}" if not snp_id[:-1].isdigit() else f"rs{snp_id[:-1]}"

def find_disease(gwas, snps, out):
    sig_res = []
    # Total GWAS SNPs.
    M = gwas['SNPS'].nunique()
    overlap = gwas.loc[gwas['SNPS'].isin(snps['snp'])].drop_duplicates()
    pp= overlap['DISEASE/TRAIT'].unique()
    print(pp)
    print(len(pp))
    sys.exit(0)
    N = overlap['SNPS'].nunique()
    probs_df = []
    snp_trait_df = []
    
    if N == 0:
        print('SNP(s) is not found in GWAS catalog...')
        print('Exiting...')
        sys.exit(0)

    for trait in tqdm(overlap['DISEASE/TRAIT'].drop_duplicates()):
        # Total trait-associated SNPs in GWAS Catalog 
        n = gwas[gwas['DISEASE/TRAIT']==trait]['SNPS'].nunique()
        trait_overlap = overlap[overlap['DISEASE/TRAIT'] == trait][['SNPS', 'DISEASE/TRAIT']]
        # Total trait-associated eQTLs
        X = trait_overlap['SNPS'].nunique()
        pval = hypergeom.sf(X-1, M, n, N)
        probs_df.append((trait, M, n, N, X, pval))
        snp_trait_df.append(trait_overlap.drop_duplicates())
    probs_cols = ['trait', 'total_gwas_snps', 'trait_snps',
                      'eqtls_in_catalog', 'trait_eqtls', 'pval']
    probs_df = pd.DataFrame(probs_df, columns=probs_cols)
    probs_df['adj_pval'] = mt.multipletests(probs_df['pval'], method='fdr_bh')[1]
    sig_df = probs_df[probs_df['adj_pval'] < 0.05]
    sig_res.append(sig_df)
    snp_trait_df = pd.concat(snp_trait_df)
    snp_trait_df = snp_trait_df.rename(columns={'DISEASE/TRAIT': 'trait'})
    snp_trait_df = snp_trait_df[snp_trait_df['trait'].isin(sig_df['trait'])]
    write_results(snp_trait_df, 'sig_trait_snps.txt', out)
    write_results(probs_df,  'enrichment.txt', out)
    sig_res = pd.concat(sig_res)
    print(sig_res)
    write_results(sig_res, 'significant_enrichment.txt',  out)

def write_results(res, fp, out):
    print(f'\tWriting {fp}...')
    os.makedirs(out, exist_ok=True)
    res.to_csv(os.path.join(out, fp), sep='\t', index=False)

        
def parse_args():
    parser = argparse.ArgumentParser(
        description='Perform disease enrichment analysis for a list of SNPs.')
    parser.add_argument(
        '-g', '--gwas', required=True,
        help='''A file containing GWAS associations''' )
    parser.add_argument(
        '-o', '--output-dir', required=True, help='Directory to write results.')
    parser.add_argument(
        '-i', '--input-file', required=True,
        help='Text file containing rsIDs in the first column.')

    return parser.parse_args()


if __name__=='__main__':
    args = parse_args()
    inputs = read_input_file(args.input_file)
    gwas = parse_gwas(args.gwas)
    find_disease(gwas, inputs, args.output_dir)
