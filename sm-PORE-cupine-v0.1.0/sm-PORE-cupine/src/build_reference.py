import sys
import numpy as np
import pandas as pd
import argparse
from pyfaidx import Fasta
from collections import OrderedDict
import h5py
import scipy.stats as spstats

def load_kmer_model(f_kmer):
    dict_kmer5 = {}
    with open(f_kmer, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            gene, mean = row[0], float(row[1])  # to meet new format
            dict_kmer5[gene] = {"mean": mean}  # remove std to meet RNA004

    return dict_kmer5

def get_expected_signal_from_kmer_model(f_fasta, f_kmer, f_output, 
                                        p_kmer=5, p_to_zscore=False):
    
    dict_kmer5 = load_kmer_model(f_kmer)
    mean_kmer = np.array([ v["mean"] for v in dict_kmer5.values()]).mean()
    
    o_fasta = Fasta(f_fasta)
    p_genes = o_fasta.keys()

    o_h5 = h5py.File(f_output, "w")

    for gene in p_genes:

        seq = str(o_fasta[gene]).upper()
        
        expected_signal = []
        for i in range(len(seq)-p_kmer):
            kmer = seq[i:i+p_kmer]
            if "N" in kmer:
                mean = mean_kmer  #np.nan
            else:
                mean = dict_kmer5[kmer]["mean"]
            expected_signal.append(mean)
    
        if p_to_zscore:
            o_h5.create_dataset(f'{gene}', data=spstats.zscore(expected_signal).tolist())
        else:
            o_h5.create_dataset(f'{gene}', data=expected_signal)

def main():
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", default=None, help="fasta file")
    parser.add_argument("-k", default=None, help="kmer model txt")
    parser.add_argument("-o", default=None, help="output file")
    parser.add_argument("-n", default=5, type=int, help="kmer size")
    parser.add_argument("-z", default=False, action="store_true", 
                        help="to z-score normalized k-mer model")
    args = parser.parse_args()
    
    f_fasta = args.i
    f_kmer = args.k
    f_output = args.o
    p_kmer = args.n
    p_to_zscore = args.z

    get_expected_signal_from_kmer_model(f_fasta, f_kmer, f_output, 
                                        p_kmer=p_kmer, p_to_zscore=p_to_zscore)


if __name__ == "__main__":

    main()