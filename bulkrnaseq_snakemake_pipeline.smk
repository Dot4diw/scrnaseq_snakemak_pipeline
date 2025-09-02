import os

config = {
                "sample_dir": "/bulkRNAseq/X101SC25079220-Z01-J001/01.RawData",
                "fastqc_threads": 2,
                "GENOMEDIR": "/database/GCF_000001405.40_GRCh38.p14",
                "annogtf": "/database/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf",
                "star_thread": 10,
                "featurecounts_thread": 10
}

#SAMPLES = [d for d in os.listdir("./") if os.path.isdir(f"./{d}") and d != ".snakemake"]
SAMPLES = ["NC-1D","NC-2D","NC-3D","KI-1D","KI-2D","KI-3D","S-1D","S-2D","S-3D","T-1D","T-2D","T-4D"]
#SAMPLES = ["NC-1D","NC-2D","NC-3D"]

rule all:
        input:
                featurecounts_rawout = expand("{sample_dir}/5_featurecounts/featurecounts_rawout_counts.txt",sample_dir = config["sample_dir"]),
                filter_rawcounts = expand("{sample_dir}/5_featurecounts/featurecounts_rawcounts.csv",sample_dir = config["sample_dir"])
                ## expand("{sample_dir}/star/{sample}_Aligned.sortedByCoord.out.bam",
                ##              sample_dir = config["sample_dir"],
                ##              sample = SAMPLES)

rule fastqc:
        input:
                expand("{sample_dir}/{sample}/{sample}_{read}.fq.gz",
                                sample_dir = config["sample_dir"],
                                sample=SAMPLES,
                                read=["1","2"])

        output:
                directory("{sample_dir}/1_fastqc")

        log:
                "{sample_dir}/1_fastqc/fastqc_run.log"

        threads: config["fastqc_threads"]

        shell:
                """
                mkdir -p {output}
                touch {log}

        fastqc -o {output} -t {threads} {input} > {log} 2>&1
        """

rule multiqc:
        input:
                "{sample_dir}/1_fastqc"
        output:
                directory("{sample_dir}/2_multiqc"),
                flag_multiqc_done = "{sample_dir}/2_multiqc/multiqc_done.flag"

        log:
                "{sample_dir}/2_multiqc/multiqc_run.log"
        shell:
                """
                mkdir -p {output}
                touch {log}

                multiqc {input} --outdir {output} > {log} 2>&1
                echo 'multiqc done!' > {flag_multiqc_done}
                """

rule fastp:
        input:
                in_r1fq = "{sample_dir}/{sample}/{sample}_1.fq.gz",
                in_r2fq = "{sample_dir}/{sample}/{sample}_2.fq.gz",
                flag_multiqc_done = "{sample_dir}/2_multiqc/multiqc_done.flag"
        output:
                out_r1fq = "{sample_dir}/3_fastp/{sample}_R1.trimmed.fq.gz",
                out_r2fq = "{sample_dir}/3_fastp/{sample}_R2.trimmed.fq.gz",
                out_json = "{sample_dir}/3_fastp/{sample}_fastp.json",
                out_html = "{sample_dir}/3_fastp/{sample}_fastp.html",
        log:
                "{sample_dir}/3_fastp/{sample}_fastp.log"
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/3_fastp"

                fastp -i {input.in_r1fq}  -I {input.in_r2fq} \
                                -o {output.out_r1fq} -O {output.out_r2fq} \
                                -j {output.out_json} -h {output.out_html} \
                                > {log} 2>&1
                """

rule star_align:
        input:
                genome_dir = config["GENOMEDIR"],
                in_r1fq = "{sample_dir}/3_fastp/{sample}_R1.trimmed.fq.gz",
                in_r2fq = "{sample_dir}/3_fastp/{sample}_R2.trimmed.fq.gz"
        output:
                outbam = "{sample_dir}/4_star/{sample}_Aligned.sortedByCoord.out.bam"
        log:
                "{sample_dir}/4_star/{sample}_runstar.log"
        params:
                out_prefix = "{sample_dir}/4_star/{sample}_"

        threads: config["star_thread"]

        shell:
                """
                mkdir -p "{wildcards.sample_dir}/4_star"

                STAR --runThreadN {threads} \
                                --genomeDir {input.genome_dir} \
                                --readFilesIn {input.in_r1fq} {input.in_r2fq} \
                                --outFileNamePrefix {params.out_prefix} \
                                --readFilesCommand zcat \
                                --outSAMtype BAM SortedByCoordinate \
                                --outReadsUnmapped Fastx \
                                --quantMode TranscriptomeSAM GeneCounts \
                                --outSAMattributes All \
                                > {log} 2>&1
                """

rule featurecounts:
        input:
                inbams = expand("{sample_dir}/4_star/{sample}_Aligned.sortedByCoord.out.bam", sample_dir = config["sample_dir"], sample = SAMPLES),
                inannogtf = config["annogtf"]

        output:
                featurecounts_rawout = "{sample_dir}/5_featurecounts/featurecounts_rawout_counts.txt",
                filter_rawcounts = "{sample_dir}/5_featurecounts/featurecounts_rawcounts.csv"
        log:
                "{sample_dir}/5_featurecounts/featurecounts_run.log"
        threads: config["featurecounts_thread"]
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/5_featurecounts"

                featureCounts -p \
                                --countReadPairs \
                                -t exon -g gene \
                                -a {input.inannogtf} \
                                --primary \
                                -T {threads} \
                                -o {output.featurecounts_rawout} {input.inbams} > {log} 2>&1

                awk -F'\t' '{{printf "%s,", $1; for(i=6; i<=NF; i++) printf "%s%s", $i, (i==NF?"\\n":",")}}' {output.featurecounts_rawout}  > {output.filter_rawcounts}
                """
