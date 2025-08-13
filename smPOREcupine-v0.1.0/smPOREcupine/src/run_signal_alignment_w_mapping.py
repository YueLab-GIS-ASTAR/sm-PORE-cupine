import sys
import pandas as pd
import ruptures as rpt
import argparse
import h5py
import numpy as np
from tslearn import metrics
import scipy.stats as spstats
import pickle as pkl
from pyfaidx import Fasta
from itertools import product
import gzip
from collections import OrderedDict

from multiprocessing import Pool

def get_mapping(f_mapping):
    dict_mapping = {}
    with open(f_mapping, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            dict_mapping[row[0]] = [row[1]] #.split(",")
    return dict_mapping

def load_signal_reference(f_reference):
    signal_references = {}
    
    o_h5 = h5py.File(f_reference)
    p_genes = list(o_h5.keys())

    for p_gene in p_genes:
        try:
            dset = o_h5[p_gene]
            signal_references[p_gene] = dset[:]
        except KeyError:
            pass
        
    return signal_references

def get_list_of_reads(f_h5):
    with h5py.File(f_h5, "r") as o_fast5:
        list_of_reads = [ iid.split("_")[1] for iid in o_fast5["/"] ]
    return list_of_reads

def get_raw_signal(o_h5, read_id):
    try:
        dset = o_h5["/read_%s/Raw/Signal" % read_id]
        raw_signal = dset[:].astype(int)[::-1]
    except KeyError:
        raw_signal = np.array([], dtype=int)
        
    return raw_signal

def align_with_subseq_DTW(signal_ref, signal_query, 
                          p_ref_use_norm, p_query_use_norm):
    if p_ref_use_norm:
        signal_ref = spstats.zscore(signal_ref)
    if p_query_use_norm:
        signal_query = spstats.zscore(signal_query)
    path, simi = metrics.dtw_subsequence_path(signal_query, signal_ref)
    final_length = path[-1][1]-path[0][1]
    return path, simi, final_length

def align_to_reference(datum):
    
    read_id, ref_genes, dict_results, p_ref_use_norm, p_query_use_norm = datum
    """
    dict_alignment_results = {"success": False,
                              "ref_gene": [],
                              "path": [],
                              "score_edist": [],
                              "query_align_len": []}
    """    
    dict_results["success_align"] = False
    dict_results["ref_gene"] = []
    dict_results["path"] = []
    dict_results["score_edist"] = []
    dict_results["query_align_len"] = []
    
    
    changepoint_signal = dict_results["changepoint_signal"]
    index_of_adapter = dict_results["index_of_adapter"]
    
    if dict_results["success"]:
        
        for ref_gene in ref_genes:
        
            signal_reference = signal_references[ref_gene]

            try:
                path, score, align_len = align_with_subseq_DTW(signal_reference,
                                                               changepoint_signal[:index_of_adapter],
                                                               p_ref_use_norm=p_ref_use_norm,
                                                               p_query_use_norm=p_query_use_norm)
                norm_score = score/len(changepoint_signal[:index_of_adapter])
                dict_results["success_align"] = True
            except IndexError:
                path, score, align_len, norm_score = None, None, None, None
            
            dict_results["ref_gene"].append(ref_gene)
            dict_results["path"].append(path)
            dict_results["score_edist"].append(score)
            dict_results["query_align_len"].append(align_len)
        
    
    return (read_id, dict_results)


def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l", default=1000000, type=int, help="Maximum Length of Trace")
    parser.add_argument("-d", default=-1, type=int, help="Depth")
    parser.add_argument("-p", default=8, type=int, help="No. of Threads/Cores")
    parser.add_argument("-m", default=None, help="list of read id to mapping file")
    parser.add_argument("-r", default=None, help="signal reference hdf5 file")
    parser.add_argument("-i", default=None, help="input segmented pkl file")
    parser.add_argument("-o", default=None, help="output file")
    parser.add_argument("-n", default=False, action="store_true", 
                        help="z-score normalize ref signal before subseq dtw alignment")
    parser.add_argument("-q", default=False, action="store_true", 
                        help="z-score normalize query signal before subseq dtw alignment")
    args = parser.parse_args()

    
    p_length = args.l
    p_depth = args.d
    p_cores = args.p
    f_reference = args.r
    f_mapping = args.m
    f_pkl = args.i
    f_output = args.o
    p_ref_use_norm = args.n
    p_query_use_norm = args.q

    # Load Data
    # ----------

    with gzip.open(f_pkl, "rb") as f:
        dict_results = pkl.load(f)

    # Get Reads
    # ----------
    read_ids = list(dict_results.keys())
    if (p_depth is not None) or (p_depth != -1):
        read_ids = read_ids[:p_depth]
    dict_rid_to_mapping = get_mapping(f_mapping)

    # Get Reference Identifiers
    # --------------------------
    sr = load_signal_reference(f_reference)
    all_ref_genes = list(sr.keys())

    # Add Data
    # ---------
    data = []
    for read_id in read_ids:
        dict_preprocess = dict_results[read_id]
        try:
            ref_genes = dict_rid_to_mapping[read_id]
        except KeyError:
            continue

        # If unmapped, use all genes
        if ref_genes[0] == "*":
            ref_genes = all_ref_genes

        data.append((read_id, ref_genes, dict_preprocess, 
                     p_ref_use_norm, p_query_use_norm))

    
    # Load Signal Reference
    # ----------------------
    def init_vars(f_reference):
        global signal_references
        signal_references = load_signal_reference(f_reference)
        

    # Run Analysis in Chunks
    # ------------------------
    with Pool(processes=p_cores, initializer=init_vars, initargs=(f_reference,)) as p:
        results = p.map(align_to_reference, data)

    dict_output = OrderedDict()
    for read_id, dict_results in results:
        dict_output[read_id] = dict_results

    # Write out as a Pickled File first
    # ----------------------------------
    with open(f_output, "wb") as f:
        pkl.dump(dict_output, f)
    
        
if __name__ == "__main__":
    
    main()
