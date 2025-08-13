#!/usr/bin/env python3

import sys
import argparse
import gzip

import numpy as np
import scipy.stats as spstats
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import OneClassSVM
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.isotonic import IsotonicRegression
from collections import Counter
import h5py
from scipy import sparse
import scipy.io as sio
import tabix
from multiprocessing import Pool
from itertools import islice

def parse_list(f_list):
    list_of_loci = []
    with open(f_list, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            p_chrom, p_start, p_end, p_strand, p_name = row[0], int(row[1]), int(row[2]), row[3], row[4]
            list_of_loci.append((p_chrom, p_start, p_end, p_strand, p_name))
    return list_of_loci

"""
def get_threshold(f_name, p_chrom):
    with h5py.File(f_name, "r") as f:
        counter_per_read = list(f[f'{p_chrom}/read_sizes'])
        filter_by_read_counts = max(counter_per_read) * 0.5
    return filter_by_read_counts
"""

"""
def get_read_info(f_name, p_contigs):
    with h5py.File(f_name, "r") as f:
        for p_contig in p_contigs:
            read_names = [ rid.decode("utf8") for rid in f[f'{p_chrom}/read_names'][:] ]
            read_index = f[f'{p_contig}/read_index'][:]
            read_sizes = f[f'{p_contig}/read_sizes'][:]
            read_start = f[f'{p_contig}/read_start'][:]
            try:
                dict_info[p_contig]
            except KeyError:
                dict_info[p_contig] = {}
            for n, i, s, p in zip(read_names, read_index, read_sizes, read_start):
                dict_info[p_contig][i] = (n, s, p)
    return dict_info
"""
"""
def get_positions(f_input, p_chrom):
    with h5py.File(f_input, "r") as f:
        pos = sorted(map(int, list(f[f'{p_chrom}/position'].keys())))
    return pos
"""

"""
def get_dataframe_by_position(f_name, p_contig, pos):
    with h5py.File(f_name, "r") as f:
        block_norm_mean = f[f'{p_contig}/position/{pos}/event_level_mean'][:]
        block_norm_stdv = f[f'{p_contig}/position/{pos}/event_level_stdv'][:]
        block_raw_mean = f[f'{p_contig}/position/{pos}/raw_mean'][:]
        block_raw_stdv = f[f'{p_contig}/position/{pos}/raw_stdv'][:]
        block_index = f[f'{p_contig}/position/{pos}/read_index'][:]
        block_dwell = f[f'{p_contig}/position/{pos}/dwell_time'][:]
        
    df = pd.DataFrame({"event_level_stdv": block_norm_stdv, "event_level_mean": block_norm_mean,
                       "raw_stdv": block_raw_stdv, "raw_mean": block_raw_mean,
                       "count": block_dwell, "index": block_index})
    df = df.set_index("index")
    return df
"""

def get_dataframe_by_position(o_tbxs, p_contig, partition):
    
    colnames = ["iid", "pos", "contig", "strand", "event_level_mean", 
                "event_level_stdv", "raw_mean", "raw_stdv", "dwell_time", "base"]
    queries = []
    for o_tbx in o_tbxs:
        queries.extend(o_tbx.query(p_contig, partition[0], partition[-1]+1))
    df = pd.DataFrame(queries, columns=colnames)
    df = df.set_index("iid")

    df["pos"] = pd.to_numeric(df["pos"])
    df["raw_stdv"] = pd.to_numeric(df["raw_stdv"], errors='coerce').fillna(0)
    df["event_level_stdv"] = pd.to_numeric(df["event_level_stdv"], errors='coerce').fillna(0)
    df["event_level_mean"] = pd.to_numeric(df["event_level_mean"], errors='coerce').fillna(0)
    return df

def get_sizes(f_sizes):
    dict_sizes = {}
    with open(f_sizes, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            dict_sizes[row[0]] = int(row[1])
    return dict_sizes

def parse_alignment_summary(f_alignments):
    dict_align = {}
    for f_alignment in f_alignments:
        with gzip.open(f_alignment, "rt") as f:
            for no, line in enumerate(f):
                if no == 0:
                    continue
                row = line.strip("\r\n").split("\t")
                start, end = int(row[4]), int(row[5])
                length = end - start
                dict_align[row[0]] = (length, start)
    return dict_align

def call_porecupine(datum):

    pos_start, pos_end, pos_mod_block, pos_unmod_block = datum

    poses = list(range(pos_start, pos_end+1))
    indices, positions = [], []
    distances, binaries = [], []
    calibrateds, uncalibrateds = [], []
    modrate = []
    all_indices = []

    for pos in poses:

        pos_mod = pos_mod_block.loc[pos_mod_block.pos==pos,:]
        pos_unmod = pos_unmod_block.loc[pos_unmod_block.pos==pos,:]
        
        # Skip machine learning if insufficient depth
        if pos_mod.shape[0] <= 5 or pos_unmod.shape[0] <= 5:
            continue
        
        if p_raw:
            pos_mod.loc[:,p_stdv] = np.log(pos_mod.loc[:,p_stdv] + 1e-3)
            pos_unmod.loc[:,p_stdv] = np.log(pos_unmod.loc[:,p_stdv] + 1e-3)

        # Train unmodified 
        model = OneClassSVM(kernel="rbf", 
                            nu=p_nu,
                            gamma=p_gamma)
        model.fit(X=pos_unmod.loc[:,[p_mean, p_stdv]])
        
        # Test modified sample for novelty detection
        labels = model.predict(pos_mod.loc[:,[p_mean, p_stdv]])

        # Relabel to convention | 0: unmod, 1: modified
        labels[labels==1] = 0
        labels[labels==-1] = 1

        #DEBUG print(Counter(labels))
        
        # Output modrate
        count_mod = np.count_nonzero(labels==1)
        count_all = labels.shape[0]
        modrate.append([pos, count_mod, count_all, count_mod/count_all])


        # distance from decision boundary
        distance = model.decision_function(pos_mod.loc[:,[p_mean, p_stdv]])
        # P(mod=1|distance) rather than P(mod=0|distance)
        uncalibrated = 1.0 - (distance - distance.min()) / (distance.max() - distance.min())

        
        # calibrate probabilities using logistic/isotonic regression
        if np.count_nonzero(labels==1) == 0 or np.count_nonzero(labels==0) == 0:
            calibrated = np.zeros(len(labels))
        else:
            if p_method == "Logistic":
                model = LogisticRegression(random_state=386)
                model.fit(X=distance.reshape(-1, 1),
                        y=labels)
                pred_proba = model.predict_proba(X=distance.reshape(-1, 1))
                calibrated = pred_proba[:,1]
            elif p_method == "LogisticCV":
                model = LogisticRegressionCV(cv=5, random_state=386)
                model.fit(X=distance.reshape(-1, 1),
                        y=labels)
                pred_proba = model.predict_proba(X=distance.reshape(-1, 1))
                calibrated = pred_proba[:,1]
            elif p_method == "Isotonic":
                model = IsotonicRegression()
                model.fit(X=distance.reshape(-1, 1),
                        y=labels)
                calibrated = model.predict(T=distance.reshape(-1, 1))
        
        # Keep only 1s
        # --------------
        indices_label, positions_label = [], []
        binaries_label, distance_label = [], []
        calibrated_label, uncalibrated_label = [], []
        all_indices_label = []

        for o_index, label, cal, uncal, dist in zip(pos_mod.index.values,
                                                    labels, calibrated,
                                                    uncalibrated,
                                                    distance):
            
            if label == 1:
                
                indices_label.append(o_index)
                positions_label.append(pos)
                binaries_label.append(label)
                calibrated_label.append(cal)
                uncalibrated_label.append(uncal)
                distance_label.append(dist)
                
            all_indices_label.append(o_index)

        indices.extend(indices_label)
        positions.extend(positions_label)
        binaries.extend(binaries_label)
        calibrateds.extend(calibrated_label)
        uncalibrateds.extend(uncalibrated_label)
        distances.extend(distance_label)
        all_indices.extend(all_indices_label)

    return (poses, indices, positions, binaries, calibrateds, uncalibrateds, distances, modrate, all_indices)

 
def chunk(arr_range, arr_size):
    arr_range = iter(arr_range)
    return iter(lambda: tuple(islice(arr_range, arr_size)), ())

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", default=None, help="chrom:start:end:strand:name (0-based-inclusive, 0-based-exclusive)")
    parser.add_argument("-l", default=None, help="BED file: list of coords")
    parser.add_argument("-m", default=None, help="mod tabix file")
    parser.add_argument("-u", default=None, help="unmod tabix file")
    parser.add_argument("-r", default=False, action="store_true", help="use raw and log?")
    parser.add_argument("-o", default=None, help="output prefix for mtx")
    parser.add_argument("-nu", default=0.001, type=float, help="nu for SVM")
    parser.add_argument("-g", default=0.1, type=float, help="gamma for SVM")
    parser.add_argument("-a", default="calibrated", 
                        help="app: calibrated | uncalibrated | distance | binary | all")
    parser.add_argument("-t", default="Logistic", help="method: Logistic | LogisticCV | Isotonic")
    parser.add_argument("-s", default=None, help="size threshold")
    parser.add_argument("-z", default=None, help="sizes file")
    parser.add_argument("-e", default=None, help="final alignment file")
    parser.add_argument("-p", default=1, type=int, help="No. of cores")
    parser.add_argument("-b", default=200, type=int, help="Block size (No. of Positions to keep in Memory)")
    args = parser.parse_args()

    p_locus = args.c
    f_list = args.l
    f_sizes = args.z
    f_alignments = args.e.split(",")

    dict_sizes = get_sizes(f_sizes)
    
    if f_alignments is not None:
        dict_align = parse_alignment_summary(f_alignments)

    list_of_loci = []
    if p_locus is not None:
        #p_attr = p_locus.split(":")
        #p_chrom, p_start, p_end, p_strand, p_name = p_attr[0], int(p_attr[1]), int(p_attr[2]), p_attr[3], p_attr[4]

        p_chrom = p_locus
        p_start = 0
        p_end = dict_sizes[p_chrom]
        p_strand = "+"
        p_name = p_locus
        
        list_of_loci.append((p_chrom, p_start, p_end, p_strand, p_name))
    else:
        if f_list is not None:
            list_of_loci = parse_list(f_list)

    p_nu = args.nu
    p_gamma = args.g
    p_app = args.a
    p_method = args.t
    p_raw = args.r
    p_threshold = args.s
    p_cores = args.p
    p_block = args.b

    f_mods = args.m.split(",")
    f_unmods = args.u.split(",")

    o_mods = [ tabix.open(f_mod) for f_mod in f_mods ]
    o_unmods = [ tabix.open(f_unmod) for f_unmod in f_unmods ]

    p_contigs = [ "%s" % locus[0] for locus in list_of_loci ]
    #dict_info = get_read_info(f_mod, p_contigs)

    print(p_contigs)
    print(list_of_loci)

    print(f_mods)
    print(f_unmods)
    
    for (p_chrom, p_start, p_end, p_strand, p_name) in list_of_loci:

        #p_contig = "%s_%s" % (p_chrom, p_strand)
        p_contig = p_chrom

        f_prefix = args.o.rstrip("/") + "/" + f'{p_name}'

        f_calibrated = f'{f_prefix}_calibrated.mtx'
        f_uncalibrated = f'{f_prefix}_uncalibrated.mtx'
        f_distance = f'{f_prefix}_distance.mtx'
        f_binary = f'{f_prefix}_binary.mtx'
        f_nid = f'{f_prefix}_binary.iids.gz'
        f_modrate = f'{f_prefix}_modrate.txt'

        if p_raw:
            p_mean, p_stdv = "raw_mean", "raw_stdv"
        else:
            p_mean, p_stdv = "event_level_mean", "event_level_stdv"

        
        ## TODO
        
        def init_vars(args):
            p_raw = args["p_raw"]
            p_nu = args["p_nu"]
            p_gamma = args["p_gamma"]
            p_mean = args["p_mean"]
            p_stdv = args["p_stdv"]

        # container to pass to each subprocess
        args = {}
        args["p_raw"] = p_raw
        args["p_nu"] = p_nu
        args["p_gamma"] = p_gamma
        args["p_mean"] = p_mean
        args["p_stdv"] = p_stdv

        

        # Calculate all necessary parameters
        # -----------------------------------

        #dict_uniq_iids, counter_uniq_iids = {}, 0


        # Start running SVM for each position
        # ------------------------------------    
        indices, positions = [], []
        calibrateds, uncalibrateds = [], []
        distances, binaries = [], []
        modrate = []

        partitions = chunk(list(range(p_start, p_end)), p_block)
        #print(p_start, p_end)
        #print(partitions)
        
        all_indices = []

        for partition in partitions:
            
            data = []
            
            pos_mod_block = get_dataframe_by_position(o_mods, p_contig, partition)
            pos_unmod_block = get_dataframe_by_position(o_unmods, p_contig, partition)

            pos_start = partition[0]
            pos_end = partition[-1]
            print(pos_start, pos_end)
            print(pos_mod_block.shape, pos_unmod_block.shape)
            
            data.append((pos_start, pos_end, pos_mod_block, pos_unmod_block))

            #print(pos_mod_block)
            #print(pos_unmod_block)
            #assert False

            """
            with Pool(processes=p_cores, initializer=init_vars, initargs=(args,)) as p:
                results = p.map(call_porecupine, data)
            """
            results = call_porecupine(data[0])
            
            #for result in sorted(results, key=lambda q: q[0][0]):
            indices.extend(results[1])
            positions.extend(results[2])
            binaries.extend(results[3])
            calibrateds.extend(results[4])
            uncalibrateds.extend(results[5])
            distances.extend(results[6])
            modrate.extend(results[7])
            all_indices.extend(list(set(results[8])))
        
        # Re-number the index
        list_old_indices = list(set(all_indices))
        dict_old_to_new_index = dict([ (i, no) for no, i in enumerate(list_old_indices) ])
        new_indices = [ dict_old_to_new_index[o_index] for o_index in indices ]

        # Write modrate file
        # -------------------
        with open(f_modrate, "w") as o_modrate:
            for mr in modrate:
                print("\t".join(map(str, mr)), file=o_modrate)


        # Write MTX files
        # ----------------
        sio.mmwrite(f_calibrated, sparse.coo_matrix((calibrateds, (new_indices, positions))))
        sio.mmwrite(f_uncalibrated, sparse.coo_matrix((uncalibrateds, (new_indices, positions))))
        sio.mmwrite(f_distance, sparse.coo_matrix((distances, (new_indices, positions))))
        sio.mmwrite(f_binary, sparse.coo_matrix((binaries, (new_indices, positions))))

        
        # Write iids file
        # ----------------
        with gzip.open(f_nid, "wt") as f:
            for n_index, o_index in enumerate(list_old_indices):
                iid = o_index
                #size, start = dict_mapping[iid][4], dict_mapping[iid][2]

                if f_alignments is None:
                    size, start = 0, 0
                else:
                    try:
                        size, start = dict_align[iid]
                    except KeyError:
                        size, start = 0, 0

                f.write("%s,%s,%s,%s\n" % (n_index, iid, size, start))
        
    
