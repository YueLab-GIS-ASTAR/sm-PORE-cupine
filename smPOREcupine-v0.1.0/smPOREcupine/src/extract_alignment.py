import numpy as np
import pandas as pd
import pickle as pkl
import gzip
import h5py
import argparse
import math
from collections import OrderedDict, Counter

"""
def parse_list(f_list_of_pass_reads, p_col=9):
    dict_list_of_pass_reads = {}
    with open(f_list_of_pass_reads, "r") as f:
        for no, line in enumerate(f):
            row = line.strip("\r\n").split("\t")
            iid = row[3]
            is_pass = True if row[p_col] == "TRUE" else False
            dict_list_of_pass_reads[iid] = is_pass
    return dict_list_of_pass_reads
"""

def parse_list(f_list_of_pass_reads):
    df = pd.read_csv(f_list_of_pass_reads, sep="\t")
    dict_list_of_pass_reads = dict(df[["read_id", "passes_filtering"]].values.tolist())
    return dict_list_of_pass_reads

def collapse_cigar(cigar_raw):
    cigar = []
    prev_state = cigar_raw[0]
    prev_count = 1
    for state in cigar_raw[1:]:
        if state == prev_state:
            prev_count += 1
        else:
            cigar.append("%s%s" % (prev_count, prev_state))
            prev_state = state
            prev_count = 1
    
    # Last
    cigar.append("%s%s" % (prev_count, prev_state))
    return "".join(cigar)
            

def path_to_cigar(path):
    dict_ref = {}
    dict_query = {}
    cigar_raw = []
    for q, r in path:

        try:
            dict_ref[r] += 1
        except KeyError:
            dict_ref[r] = 1

        try:
            dict_query[q] += 1
        except KeyError:
            dict_query[q] = 1

        if dict_ref[r] == 1 and dict_query[q] == 1:
            cigar_raw.append("M")
        elif dict_ref[r] == 1 and dict_query[q] > 1:
            cigar_raw.append("D")
        elif dict_ref[r] > 1:
            cigar_raw.append("I")
    cigar = collapse_cigar(cigar_raw)
    return cigar

def main():
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", default=None, help="output file")
    parser.add_argument("-s", default=None, help="sequencing summary file")
    parser.add_argument("-c", default=9, type=int, help="column of passes_filtering in sequencing_summary.txt")
    parser.add_argument("input", nargs="+", help="all input alignment pickled files")

    args = parser.parse_args()
    f_alignments = args.input 
    f_list_of_pass_reads = args.s  #"sequencing_summary_FAR18166_6c2beaa1.txt"
    p_col = args.c

    #f_alignment = "/home/ubuntu/projects/proj_het/analysis/sanitycheck-dtw-ribosxitch/data/alignment/NAIN3/WT/50mM_NAIN3_WT/first_pass_minimap2/batch_0.pkl.gz"
    #f_alignments = [f_alignment]

    p_sampling_frequency = 3000
    p_translocation_rate = 70
    p_kmer = 5

    dict_list_of_pass_reads = parse_list(f_list_of_pass_reads)

    print("\t".join(map(str, ["iid", "length_of_read", "length_of_polyA",
                              "chrom", "ref_start", "ref_end", "cigar", 
                              "score_min", "score_avg", "diff", 
                              "is_pass", "is_minimap2", "batch_prefix"])))
    
    for f_alignment in f_alignments:

        batch_prefix = f_alignment.split("/")[-1].replace(".pkl.gz", "")

        with gzip.open(f_alignment, "rb") as f:
            results = pkl.load(f)

        iids = list(results.keys())
        for iid in iids:

            success = results[iid]["success"]
            success_align = results[iid]["success_align"]

            # Check whether this FASTA sequence was alignable by minimap2
            is_minimap2 = False
            if "first_pass_minimap2"in f_alignment:
                is_minimap2 = True
            elif "first_pass_dtw" in f_alignment:
                is_minimap2 = False

            # Check whether this read belonged to the 
            try:
                is_pass = dict_list_of_pass_reads[iid]
            except KeyError:
                assert False

            # Check whether this read was alignable
            if success and success_align:
                
                ref_gene = results[iid]["ref_gene"]
                score_edist = np.array(results[iid]["score_edist"])
                paths = results[iid]["path"]
                index_of_polyA = results[iid]["index_of_adapter"]
                length_of_polyA = math.ceil(results[iid]["length_of_adapter"] / (p_sampling_frequency/p_translocation_rate)) + (p_kmer - 1)
                length_of_read = index_of_polyA

                if len(score_edist) > 1:
                    try:
                        score_order = np.argsort(score_edist)
                    except TypeError:
                        print(score_edist)
                        assert False
                    index_min = score_order[0]
                    index_rest = score_order[1:]

                    score_min = score_edist[index_min]
                    score_avg = score_edist[index_rest].mean()
                    diff = np.abs(score_min - score_avg)/score_avg
                else:
                    score_min = score_edist.min()
                    score_avg = -1
                    diff = -1

                # need to deal with this for equal secondary alignments
                best = np.argmin(score_edist)
                path = paths[best]
                chrom = ref_gene[best]

                if path is None:
                    ref_start = -1
                    ref_end = -1
                    cigar = None
                else:
                    path_dtw = np.array(path).reshape(-1, 2)
                    ref_start = path_dtw[:,1][0]
                    ref_end = path_dtw[:,1][-1]
                    cigar = path_to_cigar(path_dtw)

            else:
                length_of_read = -1
                length_of_polyA = -1
                chrom = "*"
                ref_start = -1
                ref_end = -1
                cigar = None
                score_min = -1
                score_avg = -1
                diff = -1

            print("\t".join(map(str, [iid, length_of_read, length_of_polyA,
                                      chrom, ref_start, ref_end, cigar, 
                                      score_min, score_avg, diff, 
                                      is_pass, is_minimap2, batch_prefix])))

if __name__ == "__main__":
    main()

