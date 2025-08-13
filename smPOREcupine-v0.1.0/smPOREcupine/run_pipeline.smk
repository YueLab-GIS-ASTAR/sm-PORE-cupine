import glob
from itertools import product


# source conda environment (if any)
# ----------------------------------
if config["c_conda_environment"] is not None:
    C_ENV  = 'eval "$(conda shell.bash hook)"\n'
    C_ENV += "conda activate " + config["c_conda_environment"]

# internally sourced files
# -------------------------
S_HDF5         = config["s_hdf5"]
S_BUILD_REF    = config["s_build_ref"]
S_CONVERT      = config["s_convert"]
S_EXTRACT      = config["s_extract"]
S_PREPROCESS   = config["s_preprocess"]
S_MAPPING      = config["s_mapping"]
S_CUDTW        = config["s_cudtw"]
S_ANALYZE      = config["s_analyze"]
S_UTILS        = config["s_utils"]
S_PORECUPINE   = config["s_porecupine"]

# user-defined binary path
# -------------------------
S_SAMTOOLS     = config["s_samtools"]
S_MINIMAP2     = config["s_minimap2"]
S_PYTHON       = config["s_python"]
S_SORT         = config["s_sort"]
S_FAST5_SUBSET = config["s_fast5_subset"]
S_BGZIP        = config["s_bgzip"]
S_TABIX        = config["s_tabix"]

# user-defined parameters
# ------------------------
DICT_DIR_NANOPORES = dict([(sample["name"], sample["dir_nanopore"]) 
                           for sample in config["samples"]])
DICT_DIR_REFERENCE = dict([(sample["name"], config["references"][sample["reference"]]["f_reference"]) 
                           for sample in config["samples"]])
DICT_DIR_HDF5      = dict([(sample["name"], config["references"][sample["reference"]]["f_hdf5"]) 
                           for sample in config["samples"]])

# user-defined optional parameters
# ---------------------------------
F_KMER                  = config["f_kmer"]
P_KMER_SIZE             = config["p_kmer_size"]
P_THRESH_MM2            = config["p_threshold_mm2"]
P_THRESH_DTW            = config["p_threshold_dtw"]
P_MAPQ_MOD              = config["p_mapq_modified"]
P_MAPQ_UNMOD            = config["p_mapq_unmodified"]

P_NORM_ZSCORE_REF       = config["p_norm_zscore_ref"]
P_NORM_ZSCORE_QUERY     = config["p_norm_zscore_query"]

P_CUDTW_USE             = config["p_cudtw_use"]
P_CUDTW_LAST_K_SIG      = config["p_cudtw_last_k_signal"]
P_CUDTW_NORM_SIG_CUTOFF = config["p_cudtw_norm_signal_cutoff"]
P_CUDTW_TOP_K_CONTIGS   = config["p_cudtw_top_k_contigs"]
P_CUDTW_TOP_K_MAPPINGS  = config["p_cudtw_top_k_mappings"]

P_SORT_MERGE_BATCH_SIZE = config["p_sort_merge_batch_size"]


# Targets
# ---------
targets_pcp_template = "results/%s/%s/%s/%s/nu_%s_gamma_%s/%s_dtw+porecupine.log"
targets_mod_template = "data/alignment/%s/%s/%s/%s/final_events.tbx.gz"
targets_umd_template = "data/alignment/%s/%s/%s/%s/final_events.tbx.gz"
targets_aln_template = "data/alignment/%s/%s/%s/%s/final_alignment.txt.gz"

targets_for_run_pcp = [ targets_pcp_template % (sample["mod"],
                                                sample["reference"],
                                                sample["name"],
                                                "final",
                                                config["porecupine"][sample["name"]]["nu"],
                                                config["porecupine"][sample["name"]]["gamma"],
                                                contig) 
                        for sample in config["samples"] 
                        for contig in config["porecupine"][sample["name"]]["contig"] ]

rule targets_run_porecupine:
    input: targets_for_run_pcp


targets_for_run_sdtw = [ "data/alignment/%s/%s/%s/final/final_alignment.txt.gz" % (sample["mod"], 
                                                                                   sample["reference"], 
                                                                                   sample["name"])
                         for sample in config["samples"] ]

rule targets_run_signal_dtw:
    input: targets_for_run_sdtw


# Supporting Function Definitions 
# ---------------------------------
def get_mod(wildcards):
    return targets_mod_template % (wildcards.mod,
                                   wildcards.reference,
                                   wildcards.sample,
                                   "final")

def get_unmod(wildcards):
    return targets_umd_template % ("Unmod",
                                   wildcards.reference,
                                   config["porecupine.sample_control_pair"][wildcards.sample],
                                   "final")

def get_sizes(wildcards):
    return config["references"][wildcards.reference]["f_sizes"]

def get_alignment(wildcards):
    return targets_aln_template % (wildcards.mod,
                                   wildcards.reference,
                                   wildcards.sample,
                                   "final")

def get_diff_threshold(wildcards):
    if wildcards.firstpass == "first_pass_minimap2":
        return P_THRESH_MM2
    elif wildcards.firstpass == "first_pass_dtw":
        return P_THRESH_DTW

def get_quality(wildcards):
    if wildcards.mod == "Unmod":
        return P_MAPQ_UNMOD
    else:
        return P_MAPQ_MOD

def get_seq_sum(wildcards):
    DIR_NANOPORE = DICT_DIR_NANOPORES[wildcards.sample]
    return " ".join(glob.glob("%s/**/sequencing_summary_*.txt" % DIR_NANOPORE, recursive=True))

def get_fast5_dir(wildcards):
    DIR_NANOPORE = DICT_DIR_NANOPORES[wildcards.sample]
    if wildcards.mod == "Unmod":
        return " ".join(glob.glob("%s/**/%s/" % (DIR_NANOPORE, "fast5_pass"), recursive=True))
    else:
        return "%s" % (DIR_NANOPORE)

def get_fastq_dir(wildcards):
    DIR_NANOPORE = DICT_DIR_NANOPORES[wildcards.sample]
    if wildcards.mod == "Unmod":
        return " ".join(glob.glob("%s/**/fastq_pass/*.fastq.gz" % DIR_NANOPORE, recursive=True))
    else:
        return " ".join(glob.glob("%s/**/*.fastq.gz" % DIR_NANOPORE, recursive=True))
        #return glob.glob("%s/**/*.fastq.gz" % DIR_NANOPORE, recursive=True)

def get_hdf5(wildcards):
    return DICT_DIR_HDF5[wildcards.sample]

def get_reference(wildcards):
    return DICT_DIR_REFERENCE[wildcards.sample]


# Rules
# ------
rule porecupine:
    input:
        alignment = get_alignment,
        mod = get_mod,
        unmod = get_unmod
    output: "results/{mod}/{reference}/{sample}/{fp}/nu_{nu}_gamma_{gamma}/{chrom}_dtw+porecupine.log"
    params:
        chrom = "{chrom}",
        gamma = "{gamma}",
        nu = "{nu}",
        sizes = get_sizes,
        prefix = "results/{mod}/{reference}/{sample}/{fp}/nu_{nu}_gamma_{gamma}/"
    shell:
        """
        {C_ENV}
        {S_PYTHON} {S_PORECUPINE} \
            -c {params.chrom} \
            -m {input.mod} \
            -u {input.unmod} \
            -o {params.prefix} \
            -nu {params.nu} \
            -g {params.gamma} \
            -s 0 \
            -z {params.sizes} \
            -e {input.alignment} \
            -b 10
        touch {output}
        """

rule summarize_final_alignment:
    input: "data/alignment/{mod}/{reference}/{sample}/final/final_events.tbx.gz"
    output: "data/alignment/{mod}/{reference}/{sample}/final/final_alignment.txt.gz"
    params:
        seqsum = get_seq_sum,
        dir_mm2 = "data/alignment/{mod}/{reference}/{sample}/first_pass_minimap2",
        dir_dtw = "data/alignment/{mod}/{reference}/{sample}/first_pass_dtw",
        r = "data/alignment/summarize_final_alignment.log"
    shell:
        """
        {C_ENV}
        python {S_EXTRACT} \
            -s {params.seqsum} \
            {params.dir_mm2}/batch_*.pkl.gz \
            {params.dir_dtw}/batch_*.pkl.gz \
            | gzip -c > {output}
        touch {params.r}
        """

rule merge_tbx:
    input: 
        mm2 = "data/alignment/{mod}/{reference}/{sample}/first_pass_minimap2/combined.log",
        dtw = "data/alignment/{mod}/{reference}/{sample}/first_pass_dtw/combined.log",
    output:
        tbx = "data/alignment/{mod}/{reference}/{sample}/final/final_events.tbx.gz"
    params:
        mm2 = "data/alignment/{mod}/{reference}/{sample}/first_pass_minimap2",
        dtw = "data/alignment/{mod}/{reference}/{sample}/first_pass_dtw",
        batch = P_SORT_MERGE_BATCH_SIZE
    shell:
        """
        {C_ENV}
        {S_SORT} --batch-size={params.batch} -m -t " " \
            -k3,3 -k2,2n -k1,1 {params.mm2}/batch_*.txt {params.dtw}/batch_*.txt \
            | tr ' ' '\t' \
            | {S_BGZIP} -c > {output.tbx}
        {S_TABIX} -s3 -b2 -e2 -0 {output.tbx}
        """

def get_txts(wildcards):
    checkpoint_output = checkpoints.subset_reads.get(**wildcards).output[0]
    return expand("data/alignment/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.txt",
           mod=wildcards.mod,
           reference=wildcards.reference,
           sample=wildcards.sample,
           firstpass=wildcards.firstpass,
           part=glob_wildcards(os.path.join(checkpoint_output, "batch_{part}.fast5")).part)

rule transitional_combine:
    input: get_txts
    output: "data/alignment/{mod}/{reference}/{sample}/{firstpass}/combined.log"
    shell:
        """
        {C_ENV}
        touch {output}
        """

rule convert_pkl_to_txt:
    input: "data/alignment/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl.gz"
    output: "data/alignment/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.txt"
    params:
        fasta = get_reference,
        fast5 = "data/subset/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.fast5",
        threshold = get_diff_threshold,
        tmp = "data/alignment/{mod}/{reference}/{sample}/{firstpass}/tmp/batch_{part}"
    shell:
        """
        {C_ENV}
        export HDF5_PLUGIN_PATH={S_HDF5}
        mkdir -p {params.tmp}
        {S_PYTHON} {S_CONVERT} \
            -i {input} \
            -f {params.fasta} \
            -x {params.fast5} \
            -t {params.threshold} \
            | {S_SORT} -T {params.tmp} -t" " -k3,3 -k2,2n -k1,1 > {output}
        """

def get_map_targets(wildcards):
    if wildcards.firstpass == "first_pass_minimap2":
        return "data/subset/%s/%s/%s/%s.map" % (wildcards.mod,
                                                wildcards.reference,
                                                wildcards.sample,
                                                wildcards.firstpass)
    elif wildcards.firstpass == "first_pass_dtw":
        if P_CUDTW_USE:
            return "data/subset/%s/%s/%s/%s.remap" % (wildcards.mod, 
                                                      wildcards.reference, 
                                                      wildcards.sample, 
                                                      wildcards.firstpass)
        else:
            return "data/subset/%s/%s/%s/%s.map" % (wildcards.mod,
                                                    wildcards.reference,
                                                    wildcards.sample,
                                                    wildcards.firstpass)

# "data/subset/{mod}/{reference}/{sample}/{firstpass}.map"
rule align_reads:
    input:
        pkl = "data/preprocess/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl.gz",
        aln = get_map_targets,
        ref = get_hdf5
    output: "data/alignment/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl.gz"
    params:
        pkl = "data/alignment/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl",
        query_use_norm = "-q" if P_NORM_ZSCORE_QUERY else ""
    threads: 4
    shell:
        """
        {C_ENV}
        export HDF5_PLUGIN_PATH={S_HDF5}
        export NUMBA_NUM_THREADS={threads}
        {S_PYTHON} {S_MAPPING} \
            -p {threads} \
            -i {input.pkl} \
            -o {params.pkl} \
            -m {input.aln} \
            -r {input.ref} {params.query_use_norm}
        gzip {params.pkl}
        """


def get_remaps(wildcards):
    glob_pattern = "data/subset/%s/%s/%s/first_pass_dtw/batch_{part}.fast5" %(wildcards.mod, wildcards.reference, wildcards.sample)
    parts = glob_wildcards(glob_pattern).part
    return [ "data/mapdtw/%s/%s/%s/first_pass_dtw/batch_%s.map" % (wildcards.mod, wildcards.reference, wildcards.sample, p) for p in parts]
    

rule merge_remap:
    input: get_remaps
    output: "data/subset/{mod}/{reference}/{sample}/first_pass_dtw.remap"
    shell:
        """
        cat {input} > {output}
        """

rule prepare_cudtw_index:
    input: 
        fasta = "data/reference/{rid}.fa",
        hdf5 = "data/reference/{rid}.hdf5"
    output: 
        index = "data/reference/{rid}.index"
    params:
        prefix = "data/reference/{rid}"
    shell:
        """
        {C_ENV}
        {S_PYTHON} {S_UTILS} reference \
            -r {input.hdf5} \
            -f {input.fasta} \
            -b {params.prefix} \
            -i {output.index} \
            -l full \
            -p -999 \
        """

def get_index(wildcards):
    fasta = DICT_DIR_REFERENCE[wildcards.sample]
    return fasta.rstrip(".fa").rstrip(".fasta") + ".index"

def get_bin(wildcards):
    fasta = DICT_DIR_REFERENCE[wildcards.sample]
    return fasta.rstrip(".fa").rstrip(".fasta") + "_0.bin"

#"data/reference/{rid}.index"
rule map_reads_by_dtw:
    input: 
        pkl = "data/preprocess/{mod}/{reference}/{sample}/{fp}/batch_{part}.pkl.gz",
        index = get_index
    output:
        remap = "data/mapdtw/{mod}/{reference}/{sample}/{fp}/batch_{part}.map"
    params:
        ref = get_bin,
        array_bin = "data/mapdtw/{mod}/{reference}/{sample}/{fp}/batch_{part}.bin",
        iid = "data/mapdtw/{mod}/{reference}/{sample}/{fp}/batch_{part}.iid",
        txt = "data/mapdtw/{mod}/{reference}/{sample}/{fp}/batch_{part}.txt.gz",
        top_k = P_CUDTW_TOP_K_CONTIGS,
        next_k = P_CUDTW_TOP_K_MAPPINGS,
        last_k = P_CUDTW_LAST_K_SIG,  # 256
        norm_sig_cutoff = P_CUDTW_NORM_SIG_CUTOFF,  # -2.0,
        query_use_norm = "-u" if P_NORM_ZSCORE_QUERY else ""
    shell:
        """
        {C_ENV}
        {S_PYTHON} {S_UTILS} dataset-multi \
            -d {input.pkl} \
            -o {params.array_bin} \
            -q {params.iid} \
            -k {params.last_k} \
            -n {params.norm_sig_cutoff} {params.query_use_norm}
        {S_CUDTW} {params.next_k} {params.ref} {params.array_bin} | gzip -c > {params.txt}
        {S_PYTHON} {S_ANALYZE} -i {input.index} -d {params.iid} -s {params.txt} -k {params.top_k} > {output.remap}
        """

rule preprocess_signal:
    input: "data/subset/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.fast5"
    output: "data/preprocess/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl.gz"
    params:
        pkl = "data/preprocess/{mod}/{reference}/{sample}/{firstpass}/batch_{part}.pkl"
    threads: 16
    shell:
        """
        {C_ENV}
        export HDF5_PLUGIN_PATH={S_HDF5}
        export NUMBA_NUM_THREADS={threads}
        {S_PYTHON} {S_PREPROCESS} \
            -p {threads} \
            -i {input} \
            -o {params.pkl}
        gzip {params.pkl}
        """

checkpoint subset_reads:
    input: "data/subset/{mod}/{reference}/{sample}/{firstpass}.map"
    output: 
        clusters = directory("data/subset/{mod}/{reference}/{sample}/{firstpass}/")
    params:
        dir_in = get_fast5_dir,
        mapping = "data/subset/{mod}/{reference}/{sample}/{firstpass}.txt"
    threads: 8
    shell:
        """
        {C_ENV}
        cat {input} | cut -f1 > {params.mapping}
        {S_FAST5_SUBSET} -f batch_ -i {params.dir_in} -s {output.clusters} -l {params.mapping} -t {threads} -r
        """

rule get_map:
    input: "data/bam/{mod}/{reference}/{sample}.sorted.bam"
    output:
        minimap2 = "data/subset/{mod}/{reference}/{sample}/first_pass_minimap2.map",
        dtw = "data/subset/{mod}/{reference}/{sample}/first_pass_dtw.map"
    params:
        quality = get_quality
    shell:
        """
        {C_ENV}
        {S_SAMTOOLS} view -F260 {input} | awk '$5>={params.quality}' | cut -f1,3 > {output.minimap2}
        {S_SAMTOOLS} view -F260 {input} | awk '$5<{params.quality}' | cut -f1,3 >> {output.dtw}
        {S_SAMTOOLS} view -f4 {input} | cut -f1,3 >> {output.dtw}
        """

# -F 2064
# remove negative strand
# and supplementary sequences
rule minimap2_align_reads:
    input:
        fastq = "data/fastq/{mod}/{sample}/final_pass+fail.fastq.gz"
    output: "data/bam/{mod}/{reference}/{sample}.sorted.bam"
    params:
        ref = get_reference
    threads: 4
    shell:
        """
        {C_ENV}
        {S_MINIMAP2} -ax map-ont -t {threads} {params.ref} {input.fastq} \
            | {S_SAMTOOLS} sort \
            | {S_SAMTOOLS} view -F2064 -bS - \
            > {output}
        {S_SAMTOOLS} index {output}
        """

rule combine_fastq:
    output:
        fastq = "data/fastq/{mod}/{sample}/final_pass+fail.fastq.gz"
    params:
        fastq_dir = get_fastq_dir
    shell:
        """
        {C_ENV}
        cat {params.fastq_dir} > {output.fastq}
        """

rule build_signal_ref:
    input:
        fasta = "data/reference/{rid}.fa",
        kmer = F_KMER
    output:
        hdf5 = "data/reference/{rid}.hdf5",
    params:
        kmer_size = P_KMER_SIZE,
        normalize = "-z" if P_NORM_ZSCORE_REF else ""
    shell:
        """
        {C_ENV}
        {S_PYTHON} {S_BUILD_REF} \
            -i {input.fasta} \
            -k {input.kmer} \
            -o {output.hdf5} \
            -n {params.kmer_size} {params.normalize}
        """
