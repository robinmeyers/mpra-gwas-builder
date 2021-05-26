# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


# import glob
# import pandas as pd
from snakemake.utils import validate, min_version
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq



##### load config and sample sheets #####


configfile: "config.yaml"
# report: "report/workflow.rst"

validate(config, schema="schemas/config.schema.yaml")

# samples = pd.read_csv(config['samplesheet']).set_index("Sample", drop=False)
# # validate(samples, schema="schemas/samples.schema.yaml")

# SAMPLES = list(samples["Sample"])
# NAMES = list(samples["ID"])
# INDICES1 = list(samples["Index1"]) if "Index1" in samples else []
# INDICES2 = list(samples["Index2"]) if "Index2" in samples else []

# LIBRARY_BASENAME = "ref/barcodes"
# LIBRARY_INDEX_DONEFILE = LIBRARY_BASENAME + ".done"

# THREADS = config['threads']

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


# include: "rules/other.smk"


# wildcard_constraints:
#     directory=".+\/",
#     sample="[^\/]+",
#     fqsuffix="_[RI]\d_\d{3}"


rule all:
    input:
       # "outs/{config[study]}_variant_ref.csv",
       # "outs/{config[study]}_oligo_list.csv",
       # "outs/" + config['study'] + "_index_snps.tsv",
       # "outs/" + config['study'] + "_gwas.tsv",
       # "outs/" + config['study'] + "_ld_snps.tsv",
       "outs/" + config['study'] + "_ld_snps_filtered.tsv",
       "outs/" + config['study'] + "_oligo_list.csv",
       "outs/" + config['study'] + "_lib_stats.csv",
       
       "outs/token"
    run:
        print("hello world! we did it!")

rule test_ld_token:
    output: "outs/token"
    script:
        "scripts/test_token.R"

rule get_gwas_snps:
    output:
        gwas = "outs/{study}_gwas.tsv"
    log: "logs/{study}_gwas_snps.log"
    script:
        "scripts/GetSNPsFromGWAS.R"

rule get_snps_in_LD:
    input:
        gwas = "outs/{study}_gwas.tsv"
    output:
        index_snps = "outs/{study}_index_snps.csv",
        ld_snps = "outs/{study}_ld_snps.tsv"
    log: "logs/{study}_ld_snps.log"
    script: "scripts/GetSNPsInLD.R"

rule intersect_epigenome:
    input:
        ld_snps = "outs/{study}_ld_snps.tsv"
    output:
        epigenome = "outs/{study}_ld_snps_epigenome.tsv",
        peak_stats = "outs/{study}_epigenome_peak_stats.csv"
    log: "logs/{study}_epigenome.log"
    script: "scripts/SNPsEpigenomeIntersect.R"

rule filter_snps:
    input:
        epigenome = "outs/{study}_ld_snps_epigenome.tsv"
    output:
        filtered_snps = "outs/{study}_ld_snps_filtered.tsv"
    log: "logs/{study}_filter_snps.log"
    script: "scripts/FilterSNPs.R"

rule build_library:
    input:
        filtered_snps = "outs/{study}_ld_snps_filtered.tsv"
    output:
        oligos = "outs/{study}_oligo_list.csv",
        barcode_ref = "outs/{study}_barcode_ref.csv",
        variant_ref = "outs/{study}_variant_ref.csv",
        random_ctrl_ref = "outs/{study}_randomcontrol_ref.csv"
    log: "logs/{study}_build_mpra_library.log"
    script: "scripts/BuildMPRALib.R"


rule figures_and_stats:
    input:
        gwas_snps = "outs/{study}_gwas.tsv",
        index_snps = "outs/{study}_index_snps.csv", 
        ld_snps = "outs/{study}_ld_snps.tsv",
        epigenome_snps = "outs/{study}_ld_snps_epigenome.tsv",
        peak_stats = "outs/{study}_epigenome_peak_stats.csv",
        filtered_snps = "outs/{study}_ld_snps_filtered.tsv",
        variant_ref = "outs/{study}_variant_ref.csv"
    output:
        stats = "outs/{study}_lib_stats.csv"
    log: "logs/{study}_figs_stats.log"
    script: "scripts/LibDesignFiguresTables.R"





