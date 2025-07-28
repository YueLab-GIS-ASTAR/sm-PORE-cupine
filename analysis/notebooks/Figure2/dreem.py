import sys
import argparse
import time
import os
import glob
import numpy as np

# Add Path for DREEM
sys.path.append("/mnt/projects/chengyza/projects/proj_het/src/dreem/code")

import BitVector_Functions
import EM_Plots
import EM_CombineRuns
import Run_EMJobs
import EM_Files

class DREEM(object):
    """
    Class wrapper for DREEM

    Attributes
    ----------
    MIN_ITS : int
        Min number of iterations per EM run (300)
    INFO_THRESH : float
        Threshold for informative bits (0.05)
    CONV_CUTOFF : float
        Diff in log like for convergence (0.5)
    NUM_RUNS : int
        Number of independent EM runs per K (10)
    MAX_K : int
        Max K to work on (2)
    SIG_THRESH : float
        Threshold to distinguish signal from noise (0.001)
    BV_THRESH  : int
        Number of consecutive 1s before filtering
    NORM_PERC_BASES : int
        Perc of bases to use for normalization (10)
    inc_TG : bool
        Include Ts and Gs?
    """
    
    def __init__(self, p_sid, input_data, input_dir, output_dir, outplot_dir,
                 NUM_CPUS=1, MIN_ITS=300, INFO_THRESH=0.05, CONV_CUTOFF=0.5, NUM_RUNS=10,
                 MAX_K=2, SIG_THRESH=0.001, BV_THRESH=4, NORM_PERC_BASES=10, inc_TG=True, p_seed=386):
        
        self.p_sid = p_sid
        self.input_data = input_data  # numpy array
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.outplot_dir = outplot_dir
        self.NUM_CPUS = NUM_CPUS
        self.MIN_ITS = MIN_ITS
        self.INFO_THRESH = INFO_THRESH
        self.BV_THRESH = BV_THRESH
        self.CONV_CUTOFF = CONV_CUTOFF
        self.NUM_RUNS = NUM_RUNS
        self.MAX_K = MAX_K
        self.SIG_THRESH = SIG_THRESH
        self.NORM_PERC_BASES = NORM_PERC_BASES
        self.inc_TG = inc_TG
        self.p_seed = p_seed
        

    def convert_to_dreem(self, data, sid):
        
        length = data.shape[1]
        seq = "".join([ "A" for _ in range(length)])
        
        dreem_data = []
        dreem_data.append("\t".join(["@ref", "%s;%s" % (sid, sid), seq]))
        dreem_data.append("\t".join(["@coordinates:length", "1,%s:%s" % (length, length)]))
        dreem_data.append("\t".join(["Query_name", "Bit_vector", "N_Mutations"]))
        
        for no, d in enumerate(data):
            no_of_mutations = int(d.sum())
            dreem_data.append("\t".join(["%s.1" % no,
                                         "".join(map(str, map(int, d))),
                                         str(no_of_mutations)]))
        return dreem_data
    
    def fit(self):

        np.random.seed(self.p_seed)

        dreem_data = self.convert_to_dreem(self.input_data, self.p_sid)
        X = EM_Files.Load_BitVectors(None, self.INFO_THRESH, self.SIG_THRESH, self.BV_THRESH,
                                     self.inc_TG, self.output_dir, bvfile_list=dreem_data)
    
        CPUS = self.NUM_CPUS #1 #4 # 24
        wind_size = 100
        norm_bases = int((wind_size * self.NORM_PERC_BASES) / 100)
        struct = False

        bvfile_basename = self.p_sid
        
        cur_BIC = float('inf')  # Initialize BIC
        BIC_failed = False  # While test is not passed
        #while not BIC_failed and K <= MAX_K:
        for K in [self.MAX_K]:
            
            print('Working on K =', K)
            RUNS = self.NUM_RUNS
            ITS = self.MIN_ITS if K != 1 else 10  # Only 10 iters for K=1
            for run in range(1, RUNS + 1):
                print('Run number:', run)
                
                Run_EMJobs.Run_EMJob(X, bvfile_basename, ITS, self.INFO_THRESH,
                                     self.CONV_CUTOFF, self.SIG_THRESH,
                                     self.outplot_dir, K, CPUS, run)
            # Processing of results from the EM runs
            EM_CombineRuns.Post_Process(bvfile_basename, K, RUNS,
                                        cur_BIC, norm_bases, struct,
                                        self.input_dir, self.outplot_dir)
            # Check BIC
            latest_BIC = EM_CombineRuns.Collect_BestBIC(bvfile_basename, K,
                                                        self.outplot_dir)
            if latest_BIC > cur_BIC:  # BIC test has failed
                BIC_failed = True
            cur_BIC = latest_BIC  # Update BIC
        
        time_taken = 0
        
        # Write params to log file
        EM_Plots.Log_File(bvfile_basename, self.NUM_RUNS, self.MIN_ITS,
                          self.CONV_CUTOFF, self.INFO_THRESH, self.SIG_THRESH, self.inc_TG,
                          norm_bases, K - 2, time_taken, self.outplot_dir)

        return self

    
    def parse_dreem_results(self, library, truths, p_threshold=0.5):
        
        p_cluster = self.MAX_K
        glob_res = glob.glob(self.outplot_dir + "K_%s" % p_cluster + "/run_*-best/Responsibilities.txt")
        f_result = glob_res[0]
        
        datum_pred = []
        dict_pred = {}
        with open(f_result, "r") as f:
            for no, line in enumerate(f):
                if no % 2 == 1:
                    row = line.strip("\r\n").split("\t")
                    rid, posteriors, bitvector = row[0], list(map(float, row[1:1+p_cluster])), row[-1]
                    cluster_no = np.argmax(posteriors)  # max posterior
                    datum_pred.append(cluster_no)
                    dict_pred[bitvector] = cluster_no

        new_datum_truth = []
        new_datum_pred = []
        new_data = []
        new_index = []
        for no, (d, t) in enumerate(zip(library, truths)):
            try:
                p = dict_pred["".join(map(str, d.tolist()))]
                new_datum_truth.append(t)
                new_datum_pred.append(p)
                new_data.append(d)
                new_index.append(no)
            except KeyError:
                pass
        new_datum_truth = np.array(new_datum_truth)
        new_datum_pred = np.array(new_datum_pred)
        new_data = np.array(new_data)
        new_index = np.array(new_index)
        
        self.labels_ = new_datum_pred
        self.truths = new_datum_truth
        self.data = new_data
        self.index = new_index

        return self    


def main(p_sid, input_file, input_dir, output_dir, outplot_dir,
         NUM_CPUS=1, MIN_ITS=300, INFO_THRESH=0.05, CONV_CUTOFF=0.5, NUM_RUNS=10,
         MAX_K=2, SIG_THRESH=0.001, NORM_PERC_BASES=10, inc_TG=True, p_seed=386):

    model = DREEM(p_sid, input_file, input_dir, output_dir, outplot_dir,
                  NUM_CPUS=1, MIN_ITS=300, INFO_THRESH=0.05, CONV_CUTOFF=0.5, NUM_RUNS=10,
                  MAX_K=2, SIG_THRESH=0.001, NORM_PERC_BASES=10, inc_TG=True, p_seed=p_seed)
    model.fit()
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", default=None, help="name of experiment")
    parser.add_argument("-p", default=None, help="prefix")
    parser.add_argument("-i", default=None, help="input file")
    parser.add_argument("-di", default=None, help="input dir")
    parser.add_argument("-do", default=None, help="output dir")
    parser.add_argument("-r", default=10, type=int, help="Number of Runs")
    parser.add_argument("-t", default=0.5, type=float, help="Convergence Threshold")
    parser.add_argument("-min_iters", default=300, type=int, help="Min number of Iterations")
    parser.add_argument("-c", default=3, type=int, help="Number of Clusters")
    parser.add_argument("-threads", default=1, type=int, help="Number of CPU threads")
    parser.add_argument("-seed", default=386, type=int, help="seed")
    parser.add_argument("-bv", default=4, type=int, help="Bit Vector threshold")
    args = parser.parse_args()
    
    dir_prefix = args.p
    dir_output = "./"
    dir_outplot = args.do
    dir_input = args.di
    f_input = args.i
    
    p_name = args.n
    p_clusters = args.c
    p_no_of_runs = args.r
    p_threshold = args.t
    p_min_iters = args.min_iters
    p_threads = args.threads
    p_seed = args.seed
    p_bv = args.bv
    
    """
    dir_prefix = "/home/chengyza/projects/proj_het/scratch/"
    dir_output = "./"
    dir_outplot = dir_prefix + "dreem_ec/output/MRPS21_WT_1AI/"
    dir_input = dir_prefix + "dreem_ec/input/"
    f_input = dir_input + "MRPS21_WT_1AI.dreem"
    """

    """ DREEM parameters
    MIN_ITS=300, INFO_THRESH=0.05, CONV_CUTOFF=0.5, NUM_RUNS=10,
    MAX_K=2, SIG_THRESH=0.001, NORM_PERC_BASES=10, inc_TG=True
    """
    
    main(p_name, f_input, dir_input, dir_output, dir_outplot,
         MAX_K=p_clusters, MIN_ITS=p_min_iters, CONV_CUTOFF=p_threshold,
         BV_THRESH=p_bv, NUM_RUNS=p_no_of_runs, NUM_CPUS=p_threads, p_seed=p_seed)
    
