import sys
import h5py
import pickle as pkl
import gzip
import argparse
from pyfaidx import Fasta
import numpy as np
import scipy.stats as spstats
from itertools import product


def get_raw_signal(o_h5, read_id):
    try:
        dset = o_h5["/read_%s/Raw/Signal" % read_id]
        raw_signal = dset[:].astype(int)[::-1]
    except KeyError:
        raw_signal = np.array([], dtype=int)
        
    return raw_signal

def get_event_level_stats(raw_signal, changepoint_signal, changepoint_length, index_of_adapter):
    cpt_mean, cpt_std = changepoint_signal[:index_of_adapter].mean(), changepoint_signal[:index_of_adapter].std()

    curr_index = 0
    raw_means, raw_stds = [], []
    event_level_mean, event_level_std = [], []
    for length in changepoint_length[:index_of_adapter]:
        segment = raw_signal[curr_index:curr_index+length]
        raw_mean, raw_std = segment.mean(), segment.std()
        norm_segment = (segment - cpt_mean)/cpt_std
        norm_mean, norm_std = norm_segment.mean(), norm_segment.std()
        raw_means.append(raw_mean)
        raw_stds.append(raw_std)
        event_level_mean.append(norm_mean)
        event_level_std.append(norm_std)
        curr_index += length
    event_level_mean = np.array(event_level_mean)
    event_level_std = np.array(event_level_std)
    return raw_means, raw_stds, event_level_mean, event_level_std

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", default=None, help="output")
    parser.add_argument("-i", default=None, help="pickled files")
    parser.add_argument("-f", default=None, help="fasta file")
    parser.add_argument("-x", default=None, help="fast5 file")
    parser.add_argument("-t", default=0.08, type=float, help="threshold margin")
    args = parser.parse_args()
    
    f_output = args.o
    f_pkl = args.i
    f_fasta = args.f
    f_fast5 = args.x
    p_threshold_margin = args.t #0.08
    
    #signal_references = get_signal_references(f_fasta)
    o_fa = Fasta(f_fasta)
    o_h5 = h5py.File(f_fast5, "r")

    with gzip.open(f_pkl, "rb") as f:

        dict_results = pkl.load(f)

    for rid, dict_result in dict_results.items():

        success = dict_result["success"]

        if not success:
            continue

        raw_signal = get_raw_signal(o_h5, rid)
        segments = np.array(dict_result["segments"])
        changepoint_signal = np.array(dict_result["changepoint_signal"])
        changepoint_length = np.array(dict_result["changepoint_length"])
        length_of_adapter = np.array(dict_result["length_of_adapter"])
        index_of_adapter = dict_result["index_of_adapter"]
        ref_gene = dict_result["ref_gene"]
        score_edist = np.array(dict_result["score_edist"])
        query_align_len = np.array(dict_result["query_align_len"])
        paths = dict_result["path"]
        
        if len(score_edist) > 1:

            score_order = np.argsort(score_edist)
            index_min = score_order[0]
            index_rest = score_order[1:]
            
            score_min = score_edist[index_min]
            #score_avg = score_edist[index_rest].mean()
            score_avg = np.percentile(score_edist[index_rest], 50)
            
            if np.abs(score_min - score_avg)/score_avg < p_threshold_margin:
                continue

        # need to deal with this for equal secondary alignments
        best = np.argmin(score_edist)
        path = paths[best]
        chrom = ref_gene[best]

        if path is None:
            continue

        #dict_positions = dict([ (k, v) for k, v in 
        #                        enumerate(signal_references[chrom]["positions"]) ])

        (raw_mean, raw_stdv,
            event_level_mean, event_level_stdv) = get_event_level_stats(raw_signal,
                                                                        changepoint_signal,
                                                                        changepoint_length,
                                                                        index_of_adapter)
        
        dict_merge = {}
        for p_query, p_ref in path:

            try:
                dict_merge[p_ref]
            except KeyError:
                dict_merge[p_ref] = {
                    "event_level_mean": [], 
                    "event_level_stdv": [],
                    "raw_mean": [],
                    "raw_stdv": [],
                    "dwell_time": []
                }
                
            try:
                dict_merge[p_ref]["raw_mean"].append(raw_mean[p_query])
                dict_merge[p_ref]["raw_stdv"].append(raw_stdv[p_query])
                dict_merge[p_ref]["event_level_mean"].append(event_level_mean[p_query])
                dict_merge[p_ref]["event_level_stdv"].append(event_level_stdv[p_query])
                dict_merge[p_ref]["dwell_time"].append(changepoint_length[:index_of_adapter][p_query])
            except IndexError:
                continue

        for k, v in dict_merge.items():
            
            #position = int(dict_positions[k])
            position = int(k)
            raw_mean = np.array(v["raw_mean"]).mean()
            raw_stdv = np.array(v["raw_stdv"]).mean()
            norm_mean = np.array(v["event_level_mean"]).mean()
            norm_stdv = np.array(v["event_level_stdv"]).mean()
            dwell_time = np.array(v["dwell_time"]).sum()
            strand = "+"
            base = str(o_fa[chrom][position+2]).upper()

            print(rid, position, chrom, strand, norm_stdv, norm_mean, raw_stdv, raw_mean, dwell_time, base)

    o_fa.close()
    o_h5.close()

if __name__ == "__main__":

    main()
