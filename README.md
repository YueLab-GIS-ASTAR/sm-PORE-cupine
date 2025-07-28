sm-PORE-cupine
============

sm-PORE-cupine ingests Oxford Nanopore's Direct RNA sequencing data to automatically detect synthetic chemical probes used to identify RNA nucleotides that are not paired in a secondary structure fashion. These unstructured bases are biochemically reactive to NAIN3, a SHAPE-derivative that is routinely used in the lab for structure probing. sm-PORE-cupine first uses mappings from both minimap2 and direct signal alignment mapping (subseqDTW) to directly obtain the signal alignment. A previously described outlier-detection algorithm utilizing one-class SVMs, POREcupine, is then used to call these synthetic RNA modifications to obtain a bit-vector read-out. The final output data is presented using a COO sparse matrix format with an accompanying metadata describing the original read identifier, mapping start position and mapped length. sm-PORE-cupine is a pip distributed snakemake pipeline to help process raw upstream data. It is highly recommended for the user to run the pipeline on a high-performance cluster.

Downstream analysis can be found in `analysis/notebooks`. If you would like to run the analysis, you can install the libraries within `analysis/lib/libsigpore` using instructions found [here](#installation)

Software Pre-requisites
------------------------

  1. Minimap2
  2. SAMtools
  3. Python 3
  4. GNU/Linux sort
  5. Oxford Nanopore's FAST5 API
  6. HTSlib bgzip and tabix
  7. CUDA Toolkit >= 10.0

Installation Requirements
--------------------------

   1. snakemake == 7.15.2
   2. pyfaidx >= 0.6.1
   3. scipy >= 1.5.4
   4. tslearn >= 0.5.2
   5. h5py == 3.8.0
   6. matplotlib >= 3.5.1
   7. numpy <= 1.21.5
   8. pandas >= 1.1.5
   9. ruptures >= 1.1.6
   10. seaborn >= 0.11.2
   11. intervaltree >= 3.1.0

This will be automatically done during `pip install`

TL;DR
------
~~~
conda create -n sm-PORE-cupine python=3.9
conda activate sm-PORE-cupine
git clone --recurse-submodules https://github.com/YueLab-GIS-ASTAR/sm-PORE-cupine
cd sm-PORE-cupine/sm-PORE-cupine-0.1.0
pip install -e .
cd ..
export CPATH=/apps/cuda/targets/x86_64-linux/include:${CPATH}
export LD_LIBRARY_PATH=/apps/cuda/targets/x86_64-linux/lib:${LD_LIBRARY_PATH}
export PATH=/apps/cuda/bin:${PATH}
cd src/cuDTW
make
cd ../..
cd analysis/lib/libsigpore
conda install -c bioconda viennarna=2.4.18
pip install -e .
cd ../../../
sm-PORE-cupine run_pipeline -c config/config_example.yaml -q slurm -r config/cluster_example.json -j 4
~~~


Installation
------------
### Conda/Python ###
Create a new conda environment. Clone the repository. Enter the directory and perform editable installation using pip. The requirments will automatically download the dependencies using pip.
~~~
conda create -n sm-PORE-cupine python=3.9
git clone --recurse-submodules https://github.com/YueLab-GIS-ASTAR/sm-PORE-cupine
cd sm-PORE-cupine/sm-PORE-cupine-0.1.0
pip install -e .
~~~
### Configuration of CUDA Toolkit to user paths and Compilation of cuDTW ###
Example configuration of CUDA LIB paths (specific to A*STAR pluto) before compiling cuDTW.
~~~
export CPATH=/apps/cuda/targets/x86_64-linux/include:${CPATH}
export LD_LIBRARY_PATH=/apps/cuda/targets/x86_64-linux/lib:${LD_LIBRARY_PATH}
export PATH=/apps/cuda/bin:${PATH}
cd src/cuDTW
make
~~~
If an nvidia GPU is unavailable or cuDTW doesn't compile successfully, you can default to using the CPU version of subseqDTW by triggering the config parameter `p_cudtw_use: false` (Default is False). However, this is only recommended for small targeted references (<10Kb).
### Configuration of Downstream Analysis Library (libsigpore) ###
Enter the directory and perform editable installation using pip. pip will download and install dependencies from the `requirements.txt` file 
~~~
cd analysis/lib/libsigpore
conda install -c bioconda viennarna=2.4.18
pip install -e .
~~~

Usage
-----
Local execution
~~~
sm-PORE-cupine run_pipeline -c <config file> -p 16
~~~
Before executing, we can utilize snakemake's dry-run feature to generate the DAG to check whether the sequence of tasks that will be executed is correctly formed
~~~
sm-PORE-cupine run_pipeline -c <config file> -p 16 --dry-run
~~~

To run on a slurm cluster
~~~
sm-PORE-cupine run_pipeline -c <config yaml file> -q slurm -r <config cluster json> -j 4
~~~

Output folder `tree` description
~~~
├── config
│   ├── cluster.json
│   ├── config_core.yaml
│   └── config_pipeline.yaml
├── data
│   ├── alignment
│   ├── bam
│   ├── fastq
│   ├── mapdtw
│   ├── bam
│   ├── preprocess
│   ├── raw
│   └── subset
└── results
~~~

| Folder | Description |
| --- | --- |
| subset | after a 1st pass mapping by minimap2, original fast5 files will be reordered into first_pass_minimap2 or first_pass_dtw. Reads that fail to be mapped by minimap2 will enter the first_pass_dtw bin |
| preprocess | Read signal traces are segmented into events |
| mapdtw | if p_cudtw_use is true, a round of signal alignment using subseqDTW will be done to rescue and obtain the mappings |
| alignment | Read signal events are signal aligned (re-squiggled) using the mappings | 

Config File
------------
The configuration file describes all the parameters and input folders or files that are required to execute the workflow. If a required file is missing, sigpore will output a snakemake error. Hence, it is recommended to run the workflow using `--dry-run` before executing the whole pipeline.

### User-defined Required Paramters ###
| Parameters | Description | Default Value
| --- | --- | --- |
| c_conda_environment | name of the conda environment being used | |
| s_minimap2 | path to minimap2 binary executable | minimap2 |
| s_samtools | path to samtools binary executable | samtools |
| s_python | path to python binary executable | python |
| s_sort | path to sort binary executable | sort |
| s_fast5_subset | path to ONT's fast5_subset | fast5_subset |
| s_bgzip | path to HTSlib's bgzip | bgzip |
| s_tabix | path to HTSlib's tabix | tabix |

### User-defined Optional Parameters ###

| Parameters | Description | Default Value
| --- | --- | --- |
| f_kmer | Path to k-mer model (Optional) | |
| p_kmer_size | k-mer size for Direct RNA | 5 |
| p_threshold_mm2 | Mapping threshold of minimap2 | 0.0 |
| p_threshold_dtw | Mapping threshold of dtw | 0.09 |
| p_mapq_modified | MAPQ threshold for modified samples | 1 |
| p_mapq_unmodified | MAPQ threshold for unmodified samples | 60 |
| p_norm_zscore_ref | Whether to normalize the reference sequence | true |
| p_norm_zscore_query | Whether to normalize the query sequence | true |
| p_cudtw_use | Use cuDTW for mapping (generally for large transcriptome reference) | false |
| p_cudtw_last_k_signal | Last k signal events to use for cuDTW mapping (from the 3' end) | 256 |
| p_cudtw_norm_signal_cutoff | Removal of outlier signal | -2.0 |
| p_cudtw_top_k_contigs | Top 20 possible contigs to use for background | 20 |
| p_cudtw_top_k_mappings | Top 200 possible mappings (contigs and positions) to use for background | 200 |
| p_sort_merge_batch_size | Maximum number of files to merge sort at once | 1000 |
