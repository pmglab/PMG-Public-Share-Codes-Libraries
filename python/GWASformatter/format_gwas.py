#!/usr/bin/env python

import os
import sys
import time
import argparse
try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit('\nNumPy and Pandas are required. If output an LDSC format file, SciPy is also required.\n')

description = '''
This Python script formats GWAS summary statistics to KGGSEE compatible format, GCTA format and
LDSC format. A SNP reference file contains columns of chromosome, basepair coordinate, SNP ID,
allele 1, allele 2, and allele 1 frequency will be loaded. The input GWAS summary statistics will
be mapped to the SNP reference by either SNP coordinates or SNP IDs. SNPs with any confliction in
coordinate, SNP ID or alleles between the GWAS summary statistics and the SNP reference will be
removed. Allele frequencies and effect sizes will be flipped to match the allele specified by the
--ref-a1-col flag. All columns of the output files will be compatible to tht SNP reference file.

The --sum-file-out flag specified file will contain the following columns:
CHR    chromosome
BP     basepair coordinate
SNP    SNP ID
A1     the effect allele
A2     the other allele
FRQ    frequency of A1
BETA   effect size of A1
SE     stderr of BETA
P      p-value
Neff   effective sample size

The --gcta-out flag specified output file will contain columns of SNP, A1, A2, FRQ, BETA, SE, P and
Neff. The --ldsc-out flag specified output file will contain columns of SNP, A1, A2, N and Z. These
two outputs can be masked by booleanizable-value columns (e.g., 0 for False and 1 for True) in the
SNP reference file.
'''

if len(sys.argv) == 1:
    print(description)
    sys.exit(f'Show help message: {sys.argv[0]} -h\n')



parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='This Python script formats GWAS summary statistics to KGGSEE compatible format, GCTA format and LDSC format.')

parser.add_argument('--sum-file-in',  default='NA',  type=str,   help='Input file. This flag is required.', required=True)
parser.add_argument('--sum-file-out', default='NA',  type=str,   help='Output file of a KGGSEE compatible format. This flag is required.', required=True)
parser.add_argument('--gcta-out',     default='NA',  type=str,   help='Output file of the GCTA format. This file should have an extension of ".ma". This flag is optional.')
parser.add_argument('--ldsc-out',     default='NA',  type=str,   help='Output file of the LDSC format. This flag is optional.')
parser.add_argument('--chrom-col',    default='NA',  type=str,   help='Name of the chromosome column. Either the (--chrom-col and --pos-col) or --snp-col flag is required.')
parser.add_argument('--pos-col',      default='NA',  type=str,   help='Name of the SNP coordinate column.')
parser.add_argument('--snp-col',      default='NA',  type=str,   help='Name of the SNP ID column. If --chrom-col and --pos-col are specified, this flag is ignored.')
parser.add_argument('--a1-col',       default='A1',  type=str,   help='Name of the effect allele column.')
parser.add_argument('--a2-col',       default='A2',  type=str,   help='Name of the non-effect allele column.')
parser.add_argument('--freq-a1-col',  default='NA',  type=str,   help='Name of the effect allele frequency column. If not specified, the reference file FRQ column will be used.')
parser.add_argument('--beta-col',     default='NA',  type=str,   help='Name of the effect size column. Either the --beta-col or --or-col flag is required.')
parser.add_argument('--beta-se-col',  default='NA',  type=str,   help='Name of the beta (log odds ratio) standard error column. Either the --beta-se-col or (--or-95l-col and--or-95u-col) flag is required.')
parser.add_argument('--or-col',       default='NA',  type=str,   help='Name of the odds ratio column. If --beta-col is specified, this flag is ignored.')
parser.add_argument('--or-95l-col',   default='NA',  type=str,   help='Name of the odds ratio 0.95 lower confidence limit column. If --beta-se-col is specified, this flag is ignored.')
parser.add_argument('--or-95u-col',   default='NA',  type=str,   help='Name of the odds ratio 0.95 upper confidence limit column. If --beta-se-col is specified, this flag is ignored.')
parser.add_argument('--p-col',        default='P',   type=str,   help='Name of the p-value column.')
parser.add_argument('--nmiss-col',    default='NA',  type=str,   help='Name of the sample size column. Either the --nmiss-col or --n flag is required.')
parser.add_argument('--neff',         default=0.0,   type=float, help='GWAS sample size. If --nmiss-col is specified, this flag is ignored.')
parser.add_argument('--info-col',     default='NA',  type=str,   help='Name of the imputation information score column. This argument is optional.')
parser.add_argument('--info-min',     default=0.9,   type=float, help='Minimum INFO score. If --info-col is not specifed, this flag is ignored.')

default_ref_file_in = f'{os.path.dirname(os.path.realpath(__file__))}/1kgEURmac5_UKBinfo9_BothAllelesMatched_hg19_dbSNPb151.lite.tsv.gz'
parser.add_argument('--ref-file-in',     default=default_ref_file_in,   type=str, help='The SNP reference file.')
parser.add_argument('--ref-chrom-col',   default='hg19_chr',            type=str, help='Name of the chromosome column of the reference file.')
parser.add_argument('--ref-pos-col',     default='hg19_bp',             type=str, help='Name of the SNP coordinate column of the reference file.')
parser.add_argument('--ref-snp-col',     default='rsid',                type=str, help='Name of the SNP ID column of the reference file.')
parser.add_argument('--ref-a1-col',      default='hg19_alt',            type=str, help='Name of the allele column of the reference file, by which the effect allele will be mapped.')
parser.add_argument('--ref-a2-col',      default='hg19_ref',            type=str, help='Name of the allele column of the reference file, by which the non-effect allele will be mapped')
parser.add_argument('--ref-freq-a1-col', default='hg19_alt_1kgEUR_frq', type=str, help='Name of the effect allele frequency column of the reference file.')
parser.add_argument('--ref-gcta-col',    default='NA',                  type=str, help='Name of a booleanizable-value column, by which the program masks the GCTA-format output. Specify "NA" to disable mask.')
parser.add_argument('--ref-ldsc-col',    default='NA',                  type=str, help='Name of a booleanizable-value column, by which the program masks the LDSC-format output. Specify "NA" to disable mask.')



def check_file_writable(fnm):
    if os.path.exists(fnm):
        # The path (file or dir) exists.
        if os.path.isfile(fnm):
            # The path is a file.
            return os.access(fnm, os.W_OK) # Return if the file is writable.
        else:
            return False # The path is a dir, so cannot write as a file.
    else:
        # The path does not exist.
        pdir = os.path.dirname(fnm)
        if not pdir:
            pdir = '.'
        return os.access(pdir, os.W_OK) # Return if the parent dir is writable.



def check_argument(args):
    args = pd.Series(args.__dict__)
    try:
        df = pd.read_csv(args.sum_file_in, sep='\s+', nrows=0)
    except FileNotFoundError:
        sys.exit(f'\nError: "{args.sum_file_in}" is not readable.\n')
    try:
        ref_df = pd.read_csv(args.ref_file_in, sep='\s+', nrows=0)
    except FileNotFoundError:
        sys.exit(f'\nError: "{args.ref_file_in}" is not readable.\n')
    if not check_file_writable(args.sum_file_out):
        sys.exit(f'\nError: "{args.sum_file_out}" is not writable.\n')

    # Initialize a dict of valid arguments. The following arguments are either required or prespecified.
    essential_col = ['ref_file_in', 'ref_chrom_col', 'ref_pos_col', 'ref_snp_col', 'ref_a1_col', 'ref_a2_col', 'ref_freq_a1_col',
                     'sum_file_in', 'sum_file_out', 'a1_col', 'a2_col', 'p_col']
    valid_args = args[essential_col]

    # Determin wether the coordinate or the ID of SNPs will be used to map alleles to the SNP reference.
    if args.chrom_col != 'NA' and args.pos_col != 'NA':
        valid_args['chrom_col'] = args.chrom_col
        valid_args['pos_col'] = args.pos_col
    else:
        if args.snp_col != 'NA':
            valid_args['snp_col'] = args.snp_col
        else:
            print('\nEither the (--chrom-col and --pos-col) or --snp-col flag is required.')
            sys.exit(f'Show help message: {sys.argv[0]} -h\n')

    log1 = [f'\n"{args.sum_file_out}" of KGGSEE compatible format will be output.']
    # Determin wether the FRQ of GWAS file or the FRQ of SNP reference will be used.
    if args.freq_a1_col != 'NA':
        valid_args['freq_a1_col'] = args.freq_a1_col
        log1.append(f'The output FRQ column will be from the "{args.freq_a1_col}" column of "{args.sum_file_in}".')
    else:
        log1.append(f'The output FRQ column will be from the "{args.ref_freq_a1_col}" column of "{args.ref_file_in}".')

    # Determin wether BETA or OR of the GWAS file will be used.
    if args.beta_col != 'NA':
        valid_args['beta_col'] = args.beta_col
        log1.append(f'The output BETA column will be from the "{args.beta_col}" column of "{args.sum_file_in}".')
    elif args.or_col != 'NA':
        valid_args['or_col'] = args.or_col
        log1.append(f'The output BETA column will be from the natural logarithm of the "{args.or_col}" column of "{args.sum_file_in}".')
    else:
        print('\nEither the --beta-col or --or-col flag is required.')
        sys.exit(f'Show help message: {sys.argv[0]} -h\n')

    # Determin wether the coordinate or the ID of SNPs will be used to map alleles to the SNP reference.
    if args.beta_se_col != 'NA':
        valid_args['beta_se_col'] = args.beta_se_col
    else:
        if args.or_95l_col != 'NA' and args.or_95u_col != 'NA':
            valid_args['or_95l_col'] = args.or_95l_col
            valid_args['or_95u_col'] = args.or_95u_col
        else:
            print('\nEither the --beta-se-col or (--or-95l-col and--or-95u-col) flag is required.')
            sys.exit(f'Show help message: {sys.argv[0]} -h\n')

    # Determin wether the N of GWAS file or the N specified by user will be used.
    if args.nmiss_col != 'NA':
        valid_args['nmiss_col'] = args.nmiss_col
        log1.append(f'The output Neff column will be from the "{args.nmiss_col}" column of "{args.sum_file_in}".')
    elif args.neff != 0:
        valid_args['neff'] = args.neff
        log1.append(f'The output Neff column will be set to {round(args.neff)}.')
    else:
        print('\nEither the --nmiss-col or --n flag is required.')
        sys.exit(f'Show help message: {sys.argv[0]} -h\n')

    # Determin wether filter SNPs by INFO.
    if args.info_col != 'NA':
        valid_args['info_col'] = args.info_col
        valid_args['info_min'] = args.info_min

    # Determin wether a file of GCTA format will be output.
    if args.gcta_out != 'NA':
        if not check_file_writable(args.gcta_out):
            sys.exit(f'\nError: "{args.gcta_out}" is not writable.\n')
        valid_args['gcta_out'] = args.gcta_out
        if args.ref_gcta_col != 'NA':
            valid_args['ref_gcta_col'] = args.ref_gcta_col
        log1.append(f'"{args.gcta_out}" of GCTA format will be output.')

    # Determin wether a file of LDSC format will be output.
    if args.ldsc_out != 'NA':
        if not check_file_writable(args.ldsc_out):
            sys.exit(f'\nError: "{args.ldsc_out}" is not writable.\n')
        try:
            from scipy.stats import norm
        except ImportError:
            sys.exit('\nSciPy is required to output an LDSC format file.\n')
        valid_args['ldsc_out'] = args.ldsc_out
        if args.ref_ldsc_col != 'NA':
            valid_args['ref_ldsc_col'] = args.ref_ldsc_col
        log1.append(f'"{args.ldsc_out}" of LDSC format will be output.')

    # Check if all specifed columns of --sum-file-in exist.
    possible_sum_col = ['chrom_col', 'pos_col', 'snp_col', 'a1_col', 'a2_col', 'freq_a1_col', 'beta_col', 'or_col', 'beta_se_col', 'or_95l_col', 'or_95u_col', 'p_col', 'nmiss_col', 'info_col']
    used_sum_col = valid_args[valid_args.index.isin(possible_sum_col)]
    err_sum_col = used_sum_col[~used_sum_col.isin(df.columns)]
    if err_sum_col.shape[0] != 0:
        sys.exit(f'\nError: {err_sum_col.to_list()} is not in the header of "{args.sum_file_in}".\n')

    # Check if all specifed columns of --ref-file-in exist.
    possible_ref_col = ['ref_chrom_col', 'ref_pos_col', 'ref_snp_col', 'ref_a1_col', 'ref_a2_col', 'ref_freq_a1_col', 'ref_gcta_col', 'ref_ldsc_col']
    used_ref_col = valid_args[valid_args.index.isin(possible_ref_col)]
    err_ref_col = used_ref_col[~used_ref_col.isin(ref_df.columns)]
    if err_ref_col.shape[0] != 0:
        sys.exit(f'\nError: {err_ref_col.to_list()} is not in the header of "{args.ref_file_in}".\n')

    # Print valid arguments and determined options to STDOUT.
    print('\nEffective settings:')
    for k in valid_args.index:
        print('--'+k.replace('_','-'), valid_args[k])
    print('\n'.join(log1), flush=True)

    return valid_args



def munge_sumstats(x):
    chr_dtype = pd.api.types.CategoricalDtype(categories=[str(a) for a in range(1,23)] + ['X'], ordered=True)

    # Read SNP reference.
    print(f'Reading {x.ref_file_in}...', flush=True)
    ref_df = pd.read_csv(x.ref_file_in, sep='\s+', dtype={x.ref_chrom_col:chr_dtype, x.ref_pos_col:pd.Int32Dtype()})
    rename_ref_col = {x.ref_chrom_col:'CHR', x.ref_pos_col:'BP', x.ref_snp_col:'SNP', x.ref_a1_col:'A1', x.ref_a2_col:'A2'}
    if 'freq_a1_col' not in x.index:
        rename_ref_col[x.ref_freq_a1_col] = 'FRQ'
    if 'ref_gcta_col' in x.index:
        rename_ref_col[x.ref_gcta_col] = 'GCTA'
    if 'ref_ldsc_col' in x.index:
        rename_ref_col[x.ref_ldsc_col] = 'LDSC'
    ref_df = ref_df[rename_ref_col.keys()].rename(rename_ref_col, axis=1)

    # Avoid the potential overriding of A1, A2 and FRQ columns of GWAS file by the SNP reference
    rename_gwas_col = {x.a1_col:'a1', x.a2_col:'a2', x.p_col:'P'}
    if 'beta_col' in x.index:
        rename_gwas_col[x.beta_col] = 'BETA'
    if 'nmiss_col' in x.index:
        rename_gwas_col[x.nmiss_col] = 'Neff'
    if 'freq_a1_col' in x.index:
        rename_gwas_col[x.freq_a1_col] = 'FRQ'
    if 'beta_se_col' in x.index:
        rename_gwas_col[x.beta_se_col] = 'SE'

    print(f'Reading {x.sum_file_in}...', flush=True)
    # Read the GWAS file and map to the SNP reference by either coordinate or ID.
    if 'chrom_col' in x.index:
        rename_gwas_col.update({x.chrom_col:'CHR', x.pos_col:'BP'})
        ref_df.set_index(['CHR','BP'], inplace=True)
        df = pd.read_csv(x.sum_file_in, sep='\s+', dtype={x.chrom_col:chr_dtype, x.pos_col:pd.Int32Dtype()})
        print(f'Read {df.shape[0]} SNPs.', flush=True)
        df = df.rename(rename_gwas_col, axis=1).set_index(['CHR','BP'])
        df = df[~df.index.duplicated(keep=False)]
        print(f'{df.shape[0]} SNPs with unique coordinate.', flush=True)
        if 'SNP' in df.columns:
            df.drop('SNP', axis=1, inplace=True)
    else:
        rename_gwas_col.update({x.snp_col:'SNP'})
        ref_df.set_index('SNP', inplace=True)
        df = pd.read_csv(x.sum_file_in, sep='\s+')
        print(f'Read {df.shape[0]} SNPs.', flush=True)
        df = df.rename(rename_gwas_col, axis=1).set_index('SNP')
        df = df[~df.index.duplicated(keep=False)]
        print(f'{df.shape[0]} SNPs with unique ID.', flush=True)
        if 'CHR' in df.columns:
            df.drop('CHR', axis=1, inplace=True)
        if 'BP' in df.columns:
            df.drop('BP', axis=1, inplace=True)


    # df = pd.concat([df, ref_df], axis=1, join='inner') # This function result in a bug. Only chr1-9 remained in some files.
    df = df[df.index.isin(ref_df.index)]
    ref_df = ref_df[ref_df.index.isin(df.index)]
    df = pd.concat([df, ref_df], axis=1)
    df = df.reset_index().sort_values(['CHR','BP'])
    print(f'{df.shape[0]} SNPs in reference.', flush=True)

    # Filter SNPs by INFO.
    if 'info_col' in x.index:
        df = df[df[x.info_col] > x.info_min]
        print(f'{df.shape[0]} SNPs with INFO > {x.info_min}')
    # Fix zero p-values.
    df.loc[df['P']<1e-300, 'P'] = 1e-300
    # Convert OR values to BETA values.
    if 'or_col' in x.index:
        df['BETA'] = np.log(df[x.or_col])
    if 'or_95l_col' in x.index:
        df['SE'] = np.abs(np.log(df[x.or_95u_col]/df[x.or_95l_col]) / 3.92)
    # Add an N column.
    if 'neff' in x.index:
        df['Neff'] = round(x.neff)

    # Filter SNPs with conflicting alleles between the GWAS file and the SNP reference.
    df['a1'] = df['a1'].str.upper()
    df['a2'] = df['a2'].str.upper()
    df['m1'] = np.all(df[['a1', 'a2']] == df[['A1', 'A2']].values, axis=1)
    df['m2'] = np.all(df[['a1', 'a2']] == df[['A2', 'A1']].values, axis=1)
    df = df[np.any(df[['m1','m2']], axis=1)]
    print(f'{df.m1.sum()} SNPs with matched alleles to reference.', flush=True)
    print(f'{df.m2.sum()} SNPs with flipped alleles to reference.', flush=True)

    # Map the FRQ and BETA of the GWAS file to match the A1 allele of the SNP reference.
    df['BETA'] = df['BETA'] * df['m1'].map({True:1, False:-1})
    if 'freq_a1_col' in x.index:
        df['FRQ'] = np.abs(df['m2'].astype(int) - df['FRQ'])


    # Output a file of KGGSEE compatibale format.
    df_kggsee = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'Neff']]
    print(f'Writing {df_kggsee.shape[0]} SNPs to {x.sum_file_out}...', flush=True)
    # The np.savetxt function works much faster than the df.to_csv method.
    np.savetxt(x.sum_file_out, df_kggsee.values, fmt='%s\t%d\t%s\t%s\t%s\t%.6f\t%.4e\t%.4e\t%.4e\t%d', header='\t'.join(df_kggsee.columns), comments='')

    # Output a file of GCTA format.
    if 'gcta_out' in x.index:
        if 'ref_gcta_col' in x.index:
            df_gcta = df.loc[df['GCTA'].astype(bool), ['SNP', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'Neff']]
        else:
            df_gcta = df[['SNP', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'Neff']]
        print(f'Writing {df_gcta.shape[0]} SNPs to {x.gcta_out}...', flush=True)
        np.savetxt(x.gcta_out, df_gcta.values, fmt='%s\t%s\t%s\t%.6f\t%.4e\t%.4e\t%.4e\t%d', header='\t'.join(df_gcta.columns), comments='')

    # Output a file of LDSC format.
    if 'ldsc_out' in x.index:
        from scipy.stats import norm
        if 'ref_ldsc_col' in x.index:
            df_ldsc = df.loc[df['LDSC'].astype(bool), ['SNP', 'A1', 'A2', 'Neff', 'P', 'BETA']]
        else:
            df_ldsc = df[['SNP', 'A1', 'A2', 'Neff', 'P', 'BETA']]
        df_ldsc.columns = ['SNP', 'A1', 'A2', 'N', 'P', 'BETA']
        df_ldsc['Z'] = norm.isf(df_ldsc['P']/2) * np.sign(df_ldsc['BETA'])
        df_ldsc = df_ldsc[['SNP', 'A1', 'A2', 'Z', 'N']]
        print(f'Writing {df_ldsc.shape[0]} SNPs to {x.ldsc_out}...', flush=True)
        np.savetxt(x.ldsc_out, df_ldsc.values, fmt='%s\t%s\t%s\t%.4f\t%d', header='\t'.join(df_ldsc.columns), comments='')



if __name__ == '__main__':
    start_time = time.time()
    print(f'\nBegin at {time.ctime()}\nChecking arguments...')

    args = parser.parse_args()
    valid_args = check_argument(args) # return a pd.Series object
    munge_sumstats(valid_args)

    time_elapsed = round(time.time() - start_time)
    print(f'\nFinish at {time.ctime()}\nTime elapsed: {time_elapsed} seconds.\n')
