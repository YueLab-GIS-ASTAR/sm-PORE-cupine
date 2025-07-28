import sys
import argparse
import numpy as np
from intervaltree import IntervalTree, Interval
import gzip

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", default=None, help="mapping index")
    parser.add_argument("-d", default=None, help="signal iid")
    parser.add_argument("-s", default=None, help="cuDTW index and scores")
    parser.add_argument("-k", default=20, type=int, help="next K transcripts")
    args = parser.parse_args()

    f_index = args.i #"bruno_et_al_2010_tx_loci.index"
    f_iid = args.d
    f_scores = args.s #"batch_0_norm_drop.txt" #"check.txt"
    
    p_num_scores = args.k

    dict_results = {}
    dict_qid_to_rid = {}

    tree = IntervalTree()

    with open(f_index, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")

            iid, start, end = row[2], int(row[3]), int(row[4])
            tree[start:end] = iid

    with open(f_iid, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            qid, rid = int(row[0]), row[1]
            dict_qid_to_rid[qid] = rid
    
    with gzip.open(f_scores, "rt") as f:
        for line in f:
            row = line.strip("\r\n").split()
            qid, gid, score = int(row[0]), int(row[1]), float(row[2])

            try:
                results = tree[gid]
                iid = list(results)[0].data
            except KeyError:
                continue

            try:
                dict_results[qid]
            except KeyError:
                dict_results[qid] = {}

            try:
                prev_score = dict_results[qid][iid]
                if score < prev_score:
                    dict_results[qid][iid] = score
            except KeyError:
                dict_results[qid][iid] = score

    for qid, d in dict_results.items():

        list_of_hits = list(d.keys())
        list_of_scores = np.array(list(d.values()))

        sorted_index = np.argsort(list_of_scores)
        best_index = sorted_index[0]

        best_score = list_of_scores[best_index]
        avg_score = list_of_scores[sorted_index[1:(p_num_scores+1)]].mean()

        print("\t".join(map(str, [dict_qid_to_rid[qid],
                                  list_of_hits[best_index], 
                                  "%.2f" % best_score,
                                  "%.2f" % avg_score,
                                  "%.3f" % ((avg_score-best_score)/avg_score),
                                  len(sorted_index)-1])))


if __name__ == "__main__":
    main()