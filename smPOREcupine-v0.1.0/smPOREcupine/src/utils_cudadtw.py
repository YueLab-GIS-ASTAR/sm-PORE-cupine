from pyfaidx import Fasta
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spstats
import h5py
from tslearn import metrics
import gzip
import pickle as pkl
from scipy.stats import zscore

import ruptures as rpt
from multiprocessing import Pool
import array as arr
import argparse
from collections import OrderedDict
import sys
import subprocess
import glob

def get_properties(o_h5, read_id):

    dset = o_h5["/read_%s/Raw/Signal" % read_id]
    seq = o_h5["/read_%s/Analyses/Basecall_1D_000/BaseCalled_template/Fastq" % read_id]
    seq = seq[()].decode("utf-8").split("\n")[1].replace("U", "T")

    key_id = "/read_%s/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment" % read_id
    
    raw_signal = dset[:].astype(int)[::-1]
    
    return raw_signal, seq

def get_info(dict_mapping, read_id):
    f_h5 = dict_mapping[read_id]
    o_h5 = h5py.File(f_h5, "r")
    raw_signal, seq = get_properties(o_h5, read_id)
    o_h5.close()
    return (None, None, None, None, raw_signal, seq)

def create_dataset(f_pkl, p_samples=100, p_kmer=64,
                   f_output_dir=None, p_norm=False, 
                   p_outlier_min=None, p_outlier_max=None):

    with gzip.open(f_pkl, "rb") as f:
        dict_results = pkl.load(f)

    read_ids = list(dict_results.keys())
    if (p_samples is not None) or (p_samples != -1):
        read_ids = read_ids[:p_samples]

    # Write out individual segemented signal 
    for read_id in read_ids:
        dict_preprocess = dict_results[read_id]

        if not dict_preprocess["success"]:
            continue

        changepoint_signal = dict_preprocess["changepoint_signal"]
        index_of_adapter = dict_preprocess["index_of_adapter"]
        if p_norm:
            changepoint_signal = zscore(changepoint_signal)
        final_signal = changepoint_signal[:index_of_adapter]

        
        if p_outlier_min is not None:
            norm_signal = zscore(final_signal)
            final_signal = final_signal[norm_signal>p_outlier_min]
        if p_outlier_max is not None:
            norm_signal = zscore(final_signal)
            final_signal = final_signal[norm_signal<p_outlier_max]
        
        if len(final_signal) < p_kmer:
            continue

        f_iid_last_bin = "%s/%s.bin" % (f_output_dir.rstrip("/"), read_id)
        with open(f_iid_last_bin, "wb") as f:
            f.write(arr.array("f", list(final_signal[-p_kmer:])))
            #f.write(arr.array("f", list(final_signal[:p_kmer])))
            #f.write(arr.array("f", list(final_signal)))

def create_dataset_multi(f_pkl, p_samples=100, p_kmer=64,
                         f_output=None, f_iids=None, p_norm=False, 
                         p_outlier_min=None, p_outlier_max=None):

    list_of_signals = []
    list_of_iids = []

    with gzip.open(f_pkl, "rb") as f:
        dict_results = pkl.load(f)

    read_ids = list(dict_results.keys())
    if (p_samples is not None) or (p_samples != -1):
        read_ids = read_ids[:p_samples]

    # Write out individual segemented signal 
    for read_id in read_ids:
        dict_preprocess = dict_results[read_id]

        if not dict_preprocess["success"]:
            continue

        changepoint_signal = dict_preprocess["changepoint_signal"]
        index_of_adapter = dict_preprocess["index_of_adapter"]
        if p_norm:
            changepoint_signal = zscore(changepoint_signal)
        final_signal = changepoint_signal[:index_of_adapter]

        
        if p_outlier_min is not None:
            if p_norm == False:
                norm_signal = zscore(final_signal)
                final_signal = final_signal[norm_signal>p_outlier_min]
            else:
                final_signal = final_signal[final_signal>p_outlier_min]
        if p_outlier_max is not None:
            if p_norm == False:
                norm_signal = zscore(final_signal)
                final_signal = final_signal[norm_signal<p_outlier_max]
            else:
                final_signal = final_signal[final_signal<p_outlier_max]
        
        if len(final_signal) < p_kmer:
            continue

        output_signal = list(final_signal[-p_kmer:])
        list_of_signals.extend([float(len(output_signal))] + output_signal)
        list_of_iids.append(read_id)

    with open(f_output, "wb") as f:
        f.write(arr.array("f", list_of_signals))
    with open(f_iids, "w") as f:
        for no, iid in enumerate(list_of_iids):
            f.write("%s\t%s\n" % (no, iid))

def create_reference(f_h5, f_fasta, f_bin, f_index, p_length,
                     p_partition=1000000, p_first=None, p_last=None, p_norm=False):
    
    o_h5 = h5py.File(f_h5, "r")
    o_fa = Fasta(f_fasta)
    p_genes = o_h5.keys()
    
    if p_length == "full":
        p_first, p_last = None, None
    elif p_length.startswith("first"):
        p_first = int(p_length.split("-")[1])
    elif p_length.startswith("last"):
        p_last = int(p_length.split("-")[1])

    if p_partition == -999 and p_first is None and p_last is None:  # full length seq fit into 1 partition
        p_partition = sum([ len(v) for v in o_fa.values() ])
    elif p_partition == -1 and p_last is not None:  # size of last partial seq fit into 1 partition
        p_partition = sum([ p_last for v in o_fa.values() ])
    elif p_partition == -9 and p_first is not None:  # size of first partial seq fit into 1 partition
        p_partition = sum([ p_first for v in o_fa.values() ])
    else:
        assert False

    list_of_indices = []
    dict_concat_signal = OrderedDict()
    curr_global_index = 0
    curr_local_index = 0
    curr_partition = 0
    for no, p_gene in enumerate(p_genes):

        if p_norm:
            signal = zscore(o_h5[p_gene][:]).tolist()
        else:
            signal = list(o_h5[p_gene[:]])
        
        if p_first is not None and p_last is not None:
            print("first and last parameter can't be both initialized")
            assert False
        elif p_first is None and p_last is None:
            signal = signal
        elif p_first is not None:
            signal = signal[:p_first]
        elif p_last is not None:
            signal = signal[-p_last:]
        
        if (curr_local_index + len(signal)) > p_partition:
            curr_local_index = 0
            curr_partition += 1
            
        # store signal into partitions
        try:
            dict_concat_signal[curr_partition].extend(signal)
        except KeyError:
            dict_concat_signal[curr_partition] = signal

        # store information for index
        list_of_indices.append((no, curr_partition, p_gene,
                                curr_global_index, curr_global_index+len(signal),
                                curr_local_index, curr_local_index+len(signal)))
        
        curr_global_index = curr_global_index + len(signal)
        curr_local_index = curr_local_index + len(signal)

    for curr_part, signals in dict_concat_signal.items():
        
        with open(f_bin + ("_%s.bin" % curr_part), "wb") as f:
            f.write(arr.array("f", signals))

    with open(f_index, "w") as f:
        for l in list_of_indices:
            print("\t".join(map(str, l)), file=f)

def main():
    
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command",
                                       help="sub-command help")

    # Arguments for Reference
    parser_create_ref = subparsers.add_parser("reference", help="Create Reference")
    parser_create_ref.add_argument("-r", default=None, help="hdf5 ref")
    parser_create_ref.add_argument("-f", default=None, help="fast ref")
    parser_create_ref.add_argument("-b", default=None, help="output prefix")
    parser_create_ref.add_argument("-i", default=None, help="output index")
    parser_create_ref.add_argument("-l", default=None, help="length")
    parser_create_ref.add_argument("-p", default=1000000, type=int, help="partition size")
    parser_create_ref.add_argument("-u", default=False, action="store_true", help="normalize?")

    # Arguments for Dataset
    parser_create_dataset = subparsers.add_parser("dataset", help="Create Dataset")
    parser_create_dataset.add_argument("-d", default=None, help="input pkl file")
    parser_create_dataset.add_argument("-o", default=None, help="output directory")
    parser_create_dataset.add_argument("-s", default=None, type=int, help="No. of reads to randomly get")
    parser_create_dataset.add_argument("-k", default=256, type=int, help="K-mer size / Length of Seed")
    parser_create_dataset.add_argument("-u", default=False, action="store_true", help="normalize?")
    parser_create_dataset.add_argument("-x", default=None, help="remove outlier: max normalized signal")
    parser_create_dataset.add_argument("-n", default=None, help="remove outlier: min normalized signal")

    # Arguments for Dataset Multi
    parser_create_multi = subparsers.add_parser("dataset-multi", help="Create Multiple Signals in Dataset")
    parser_create_multi.add_argument("-d", default=None, help="input pkl file")
    parser_create_multi.add_argument("-o", default=None, help="output bin")
    parser_create_multi.add_argument("-q", default=None, help="output iid")
    parser_create_multi.add_argument("-s", default=None, type=int, help="No. of reads to randomly get")
    parser_create_multi.add_argument("-k", default=256, type=int, help="K-mer size / Length of Seed")
    parser_create_multi.add_argument("-u", default=False, action="store_true", help="normalize?")
    parser_create_multi.add_argument("-x", default=None, help="remove outlier: max normalized signal")
    parser_create_multi.add_argument("-n", default=None, help="remove outlier: min normalized signal")
    args = parser.parse_args()

    if args.command == "reference":

        f_h5 = args.r
        f_fasta = args.f
        f_bin = args.b
        f_index = args.i
        p_length = args.l
        p_partition = args.p
        p_norm = args.u

        create_reference(f_h5, f_fasta, f_bin, f_index, p_length,
                         p_partition=p_partition, p_norm=p_norm)


    elif args.command == "dataset":

        f_pkl = args.d
        f_output_dir = args.o
        p_samples = args.s
        p_kmer = args.k
        p_norm = args.u
        p_outlier_min = args.n
        p_outlier_max = args.x

        if p_outlier_min is not None:
            p_outlier_min = float(p_outlier_min)
        if p_outlier_max is not None:
            p_outlier_max = float(p_outlier_max)

        create_dataset(f_pkl, p_samples=p_samples, p_kmer=p_kmer,
                       p_norm=p_norm,
                       f_output_dir=f_output_dir, 
                       p_outlier_min=p_outlier_min, 
                       p_outlier_max=p_outlier_max)
        

    elif args.command == "dataset-multi":

        f_pkl = args.d
        f_output = args.o
        f_iid = args.q
        p_samples = args.s
        p_kmer = args.k
        p_norm = args.u
        p_outlier_min = args.n
        p_outlier_max = args.x

        if p_outlier_min is not None:
            p_outlier_min = float(p_outlier_min)
        if p_outlier_max is not None:
            p_outlier_max = float(p_outlier_max)

        create_dataset_multi(f_pkl, p_samples=p_samples, p_kmer=p_kmer,
                             p_norm=p_norm,
                             f_output=f_output, 
                             f_iids=f_iid,
                             p_outlier_min=p_outlier_min, 
                             p_outlier_max=p_outlier_max)
        

if __name__ == "__main__":
    
    main()
