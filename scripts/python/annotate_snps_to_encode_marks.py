import pandas as pd
import argparse
import sys
import os


def split_df(encode, input_file, chunk_size):
    chunks = [input_file[i:i+chunk_size] for i in range(0,len(input_file),chunk_size)]
    snps_merged = []
    for i in chunks:
        chunks_df = pd.DataFrame(i)
        snps_processed = process_snps(encode, chunks_df)
        snps_merged.append(snps_processed)
        #sys.exit(0)
    out = pd.concat(snps_merged)
    return out

def read_encode_file(encode_file):
    colnames = ["chr", "start", "end", "histone"]
    encode = pd.read_csv(encode_file, sep = "\t", header=None, names = colnames)
    return encode

def read_input_files(input_file):
    inputs = pd.read_csv(input_file, sep = "\t", header=0)
    return inputs

def process_snps(encode, snps):
    snps_merged = pd.merge(encode, snps, left_on = "chr", right_on = "snp_chr")
    snps_res = snps_merged.loc[(snps_merged.snp_locus >= snps_merged.start) & (snps_merged.snp_locus <= snps_merged.end)]
    return snps_res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-e", "--encode_file", required=True, help="A tab-delimited text file containing four columns", dest="encode")
    parser.add_argument("-s", "--input_snps", required=True, help="A tab_delimited text file containing three columns: rsid's, chromosome(e.g., 'chrX') and position(e.g.,'100345')(REQUIRED)", dest= "snps")
    parser.add_argument("-o", "--out_dir", required=True, help="directory to write results (REQUIRED)", dest = "out")

    global args

    args=parser.parse_args()

    if args.snps == 'none':
        message = '''
        ERROR: you must specify input file with rsid's, chr, and location.
        For more details: use -h flag 
        '''
    if args.encode == 'none':
        message = '''
        ERROR: you must specify encode file '''
    if args.out == 'none':
        message = '''
        ERROR: you must specify out directory'''
    if not os.path.isdir(args.out):
        os.makedirs(args.out)
    ## read files

    input_file = read_input_files(args.snps)
    encode = read_encode_file(args.encode)
    
    if input_file.shape[0] > 10000:
        snps_res = split_df(encode, input_file, 10000)
    else:
        snps_res = process_snps(encode, input_file)
    
    
    ## process data
    #snps_res = process_snps(encode, chunks)

    ## write outputs 
    snps_ofile = os.path.join(args.out, "snp_annotation.txt")

    snps_res.to_csv(snps_ofile, sep = "\t", header=True, index=False)

