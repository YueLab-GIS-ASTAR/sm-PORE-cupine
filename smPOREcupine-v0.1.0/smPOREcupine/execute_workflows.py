import argparse
import sys
import yaml
import snakemake
import os

from . import _program
thisdir = os.path.abspath(os.path.dirname(__file__))

from .src.build_reference import get_expected_signal_from_kmer_model

def main():

    parser = argparse.ArgumentParser(prog=_program,
                                     description=__doc__)
    subparsers = parser.add_subparsers(dest="subcommand",
                                       help="sub-command help")
    
    # Add subparsers
    # ---------------
    parser_run_pipeline = subparsers.add_parser("run_pipeline", help="Run Pipeline")
    parser_run_porecupine = subparsers.add_parser("run_porecupine", help="Run Porecupine")
    parser_run_signal_dtw = subparsers.add_parser("run_signal_dtw", help="Run Signal DTW")
    parser_build_signal_ref = subparsers.add_parser("build_signal_ref", help="Build Signal Ref")
    
    # Add argument description for the following
    # -------------------------------------------
    for sparser in (parser_run_pipeline, 
                    parser_run_porecupine,
                    parser_run_signal_dtw):  # Same helper description
        
        sparser.add_argument("-c", "--config",
                             required=True, 
                             default="config/config_example.yaml",
                             help="Path to config file")
        sparser.add_argument("-s", "--snakefile",
                             default=None,
                             help="Path to Snakefile (Optional)")
        sparser.add_argument("-q", "--cluster",
                             default=None,
                             help="cluster template")
        sparser.add_argument("-r", "--cluster_config", 
                             default=None, 
                             help="cluster configuration file")
        sparser.add_argument("-j", "--no_of_jobs", default=1, type=int,
                             help="No. of Jobs submitted if run in a cluster")
        sparser.add_argument("-p", "--no_of_cores", default=1, type=int,
                             help="No. of Max Cores if run locally")
        sparser.add_argument("-n", "--dry-run", action="store_true", 
                             help="Perform a snakemake dry run")
        sparser.add_argument("-f", "--force", action="store_true",
                             help="Force snakemake to ignore errors and run")
        sparser.add_argument("-w", "--latency-wait", default=120, type=int,
                             help="Latency wait in seconds for files to turn up before declaring an error")

    # Add argument description for build_signal_ref
    # ----------------------------------------------
    parser_build_signal_ref.add_argument("-k", "--kmer_model", 
                                         default=None, help="k-mer model txt (Optional)")
    parser_build_signal_ref.add_argument("-i", "--input", 
                                         default=None, required=True, help="FASTA reference")
    parser_build_signal_ref.add_argument("-o", "--output",
                                         default=None, required=True, help="Signal reference (HDF5)")
    parser_build_signal_ref.add_argument("-n", "--kmer_size",
                                         default=5, help="k-mer size")
    parser_build_signal_ref.add_argument("-z", "--to_zscore", 
                                         default=False, action="store_true", 
                                         help="to z-score normalized k-mer model")

    # Let's Go
    # ---------
    args = parser.parse_args()

    
    if args.subcommand in ("run_pipeline", "run_porecupine", "run_signal_dtw"):

        subcommand_to_targets = {"run_pipeline": "targets_run_porecupine",
                                 "run_porecupine": "targets_run_porecupine",
                                 "run_signal_dtw": "targets_run_signal_dtw"}

        p_subcommand = args.subcommand
        p_targets = subcommand_to_targets[p_subcommand]
        p_config_core_file = "%s/%s" % (thisdir.rstrip("/"), "config/config_core.yaml")
        p_config_user_file = args.config
        p_snakefile = ("%s/%s" % (thisdir.rstrip("/"), 
                                  "run_pipeline.smk") if args.snakefile is None else args.snakefile )
        p_cluster = args.cluster
        p_cluster_config = args.cluster_config
        p_no_of_jobs = args.no_of_jobs
        p_no_of_cores = args.no_of_cores
        sm_latency_wait = args.latency_wait
        sm_dry_run = args.dry_run
        sm_force = args.force

        if os.path.isfile(p_config_core_file):
            with open(p_config_core_file, "r") as f:
                p_config_core = yaml.safe_load(f)

            # source the following
            for key, value in p_config_core.items():
                if (key.startswith("s_") or key.startswith("f_")):
                    p_config_core[key] = "%s/%s" % (thisdir.rstrip("/"), value)
            
        if os.path.isfile(p_config_user_file):
            with open(p_config_user_file, "r") as f:
                p_config_user = yaml.safe_load(f)

            # toggle overrides
            for key, value in p_config_user.items():
                if value is None:
                    p_config_user[key] = p_config_core[key]
            
        else:
            assert("[Required] Configuration file doesn't exist")
        
        # Merge dictionary
        p_config = {**p_config_core, **p_config_user}
        
        if p_cluster == "slurm":  # override
            p_cluster = "sbatch --mail-user={cluster.email} --nodes=1 --ntasks=1 --cpus-per-task={cluster.cpu} --time={cluster.time} --partition={cluster.queue} --job-name={cluster.name} --output={cluster.name}.o%j --error={cluster.name}.e%j"

        # Run snakemake
        status = snakemake.snakemake(
            p_snakefile,
            config=p_config,
            cluster=p_cluster,
            cluster_config=p_cluster_config,
            cores=p_no_of_cores,
            nodes=p_no_of_jobs,
            targets=[p_targets],
            printshellcmds=True,
            dryrun=sm_dry_run,
            forceall=sm_force,
            latency_wait=sm_latency_wait
        )

        return 0 if status else 1

    elif args.subcommand == "build_signal_ref":

        f_fasta = args.input
        f_hdf5 = args.output
        f_kmer = args.kmer_model
        p_kmer = args.kmer_size
        p_to_zscore = args.to_zscore

        get_expected_signal_from_kmer_model(f_fasta, f_kmer, f_output, p_kmer, p_to_zscore)

if __name__ == "__main__":

    main()