# Before use, the following config variable needs to be configured properly.
# 1).The `samplepath` of the `config` configuration generally refers to the nearest level directory in our data handover directory that contains *.fg.gz. Usually,
#   there are different sample name folders below. The *.fq.gz file is located under these sample directory names.
#   Just configure the `samplepath` directory to the one where the sample name folders are located.
# 2).`workpath` refers to the directory where the current script is located, which is the top-level directory for the entire process output.
# 3).`genome` refers to the directory where the reference genome index is located.
# 4).`threads` threads configure the number of threads required for a single task, and the combination of `snakemake -c N` determines how many tasks are parallel.
#   N/threads equals the number of parallel tasks.

config = {
                "samplepath": "rawdata/Fq",
                "workpath": "rawdata/Fq/0_analysis/0_align",
                "genome": "/database/REFDB_OF_DNBC4TOOLS/GRCm39",
                "threads": 12
                }

SAMPLES = ["CTLR_A","CTLR_B","CTRL_C","CTRL_D","KO_A","KO_B","KO_C","KO_D"]

rule all:
        input:
                 expand("{workpath}/{sample}/{sample}_dnbc4tools_rna_run_done.flag",
                                 workpath = config["workpath"],
                                 sample = SAMPLES)

rule dnbc4tools_rna_run:
        input:
                sampledir=directory(config["samplepath"] + "/{sample}_cDNA")

        threads:config["threads"]

        params:
                genome=config["genome"]

        output:
                flag_done = "{workpath}/{sample}/{sample}_dnbc4tools_rna_run_done.flag"

        log:
                dnbc4toolslog="{workpath}/{sample}/{sample}_dnbc4tools_rna_run.log"

        shell:
                "./dnbc4tools_align_for_snamake.sh {input.sampledir} {params.genome} {threads} > {log.dnbc4toolslog} 2>&1 \
                                && touch {output.flag_done} && echo {wildcards.sample}: dnbc4tools rna run all done! > {output.flag_done}"
