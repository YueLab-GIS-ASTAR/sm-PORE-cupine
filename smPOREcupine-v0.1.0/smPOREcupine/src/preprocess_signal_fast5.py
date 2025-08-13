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

def get_sizes(f_sizes):
    dict_sizes = {}
    with open(f_sizes, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            dict_sizes[row[0]] = int(row[1])
    return dict_sizes


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

def segment_raw_signals(datum, p_method="pelt-rbf-window", p_window=5000):
    read_id, raw_signal = datum
    if p_method == "pelt-rbf-window":
        segments = []
        num_intervals = len(raw_signal)//p_window
        for i in range(num_intervals):
            if i == num_intervals-1:
                model = rpt.KernelCPD(kernel="rbf").fit(raw_signal[i*p_window:])
            else:
                model = rpt.KernelCPD(kernel="rbf").fit(raw_signal[i*p_window:(i+1)*p_window])
            segments.extend([ (i*p_window) + j for j in model.predict(pen=2) ])
    elif p_method == "pelt-rbf":
        model = rpt.KernelCPD(kernel="rbf").fit(raw_signal)
        segments = model.predict(pen=2)
    elif p_method == "pelt-linear":
        model = rpt.KernelCPD(kernel="linear").fit(raw_signal)
        segments = model.predict(pen=40000)
        
    return (read_id, raw_signal, segments)

def collapse_signals(raw_signal, segments):
    changepoint_signal = []
    changepoint_length = []
    changepoint_order = []
    changepoint_index = []
    for i in range(1, len(segments)):
        if i == 1:
            changepoint_signal.append(raw_signal[:segments[i]].mean())
            changepoint_length.append(segments[i]-0)
            changepoint_order.append(segments[i-1])
            changepoint_index.append(i-1)
        else:
            changepoint_signal.append(raw_signal[segments[i-1]:segments[i]].mean())
            changepoint_length.append(segments[i]-segments[i-1])
            changepoint_order.append(segments[i-1])
            changepoint_index.append(i-1)
    
    changepoint_signal = np.array(changepoint_signal)
    changepoint_length = np.array(changepoint_length)
    changepoint_order = np.array(changepoint_order)
    changepoint_index = np.array(changepoint_index)

    return changepoint_signal, changepoint_length, changepoint_order, changepoint_index

'''
def remove_polyA_adapter(segments, changepoint_length, changepoint_signal):
    # Find poly-A adapter
    length_of_adapter = 0
    index = None
    changepoint_signal = spstats.zscore(changepoint_signal)
    for index in np.argsort(changepoint_length)[::-1]:
        if changepoint_signal[index] >= 0 and changepoint_signal[index] <= 2: #600-850:
            #DEBUG print(segments[index], segments[index+1], abs(segments[index]-segments[index+1]), index)
            length_of_adapter = abs(segments[index]-segments[index+1])
            break
        
    return length_of_adapter, index


def is_polyA_detection_correct(segments, length_of_adapter, index,
                               p_min_polyA=350, p_min_trace=1000):
    if length_of_adapter <= p_min_polyA:  # remove short polyA traces (real polyAs)
        return False
    if (index is None) or (segments[index] <= p_min_trace): # remove overly short remainder traces
        return False
    
    return True
'''

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)

def remove_polyA_adapter(segments, changepoint_length, changepoint_signal, 
                         changepoint_order, changepoint_index,
                         polyA_lower=0.5077654502073372,
                         polyA_upper=1.172873383054103,
                         p_min_polyA=10, 
                         p_translocation_rate=70, 
                         p_freq=3000,
                         p_approx_half_length_of_adapter=30,
                         p_stepsize=3):
    
    changepoint_signal = spstats.zscore(changepoint_signal)

    candidates_length = changepoint_length[(changepoint_signal >= polyA_lower) & (changepoint_signal <= polyA_upper)]
    candidates_order = changepoint_order[(changepoint_signal >= polyA_lower) & (changepoint_signal <= polyA_upper)]
    candidates_index = changepoint_index[(changepoint_signal >= polyA_lower) & (changepoint_signal <= polyA_upper)]

    candidates_length_group = []
    candidates_index_group = []
    for c in consecutive(candidates_index, stepsize=p_stepsize):
        if len(c) > 0:
            if c[-1] < (len(segments) - p_approx_half_length_of_adapter):
                candidates_index_group.append(c)
                candidates_length_group.append(changepoint_length[c].sum())

    candidates_length_group = np.array(candidates_length_group)
    candidates_index_group = np.array(candidates_index_group, dtype=object)

    p_threshold = ((p_min_polyA-5+1)/p_translocation_rate) * p_freq

    if len(candidates_index_group[candidates_length_group >= p_threshold]) > 0:
        consec_order = candidates_index_group[candidates_length_group >= p_threshold][-1]
        length_of_polyA = candidates_length_group[candidates_length_group >= p_threshold][-1]
        index_of_polyA = consec_order[0]
    else:
        consec_order = None
        length_of_polyA = None
        index_of_polyA = None

    return length_of_polyA, index_of_polyA

def is_polyA_detection_correct(segments, index_of_polyA,
                               p_min_trace=1000):
    if (index_of_polyA is None) or (segments[index_of_polyA] <= p_min_trace): # remove overly short remainder traces
        return False
    return True

def preprocess_nanopore_reads(datum):
    
    dict_preprocess_results = {"success": False,
                               "segments": None,
                               "changepoint_signal": None,
                               "changepoint_length": None,
                               "length_of_adapter": None,
                               "index_of_adapter": None}
    
    read_id, raw_signal, segments = segment_raw_signals(datum)
    
    (changepoint_signal, changepoint_length, 
     changepoint_order, changepoint_index) = collapse_signals(raw_signal,
                                                              segments)
    
    length_of_adapter, index_of_adapter = remove_polyA_adapter(segments,
                                                    changepoint_length,
                                                    changepoint_signal, 
                                                    changepoint_order, 
                                                    changepoint_index)

    is_polyA = is_polyA_detection_correct(segments, index_of_adapter)

    if is_polyA:
        dict_preprocess_results["success"] = is_polyA
        dict_preprocess_results["segments"] = segments
        dict_preprocess_results["length_of_adapter"] = length_of_adapter
        dict_preprocess_results["index_of_adapter"] = index_of_adapter
        dict_preprocess_results["changepoint_signal"] = changepoint_signal
        dict_preprocess_results["changepoint_length"] = changepoint_length
    
    return (read_id, dict_preprocess_results)

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l", default=1000000, type=int, help="Maximum Length of Trace")
    parser.add_argument("-d", default=-1, type=int, help="Depth")
    parser.add_argument("-p", default=8, type=int, help="No. of Threads/Cores")
    parser.add_argument("-i", default=None, help="input file")
    parser.add_argument("-o", default=None, help="output file")
    args = parser.parse_args()

    
    p_length = args.l
    p_depth = args.d
    p_cores = args.p
    f_h5 = args.i
    f_output = args.o
    
    # Load Data
    # ----------
    read_ids = get_list_of_reads(f_h5)
    if (p_depth is not None) or (p_depth != -1):
        read_ids = read_ids[:p_depth]
    
    with h5py.File(f_h5, "r") as o_h5:
        data = [ (read_id, get_raw_signal(o_h5, read_id)) for read_id in read_ids ]
        data = [ (read_id, raw_signal) for read_id, raw_signal in data
                    if (len(raw_signal) > 0) and (len(raw_signal) <= p_length) ]

        new_data = []
        for read_id, raw_signal in data:
            try:
                new_data.append((read_id, raw_signal))
            except KeyError:
                pass
        data = new_data

    # Run Analysis in Chunks
    # ------------------------
    with Pool(processes=p_cores) as p:
        results = p.map(preprocess_nanopore_reads, data)

    dict_results = OrderedDict()
    for r in results:
        dict_results[r[0]] = r[1]

    # Write out as a Pickled File first
    # ----------------------------------
    with open(f_output, "wb") as f:
        pkl.dump(dict_results, f)

if __name__ == "__main__":
    main()
