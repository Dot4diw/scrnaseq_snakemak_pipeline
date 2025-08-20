import os

config = {
                "sample_dir": "/X101SC25079220-Z01-J001/01.RawData",
                "result_dir": "qc_results",
                "log_dir": "logs",
                "fastqc_threads": 4,
                "multiqc_title": "qc_summary",
                "GENOMEDIR": "/database/GCF_000001405.40_GRCh38.p14",
                "annogtf": "/database/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf",
                "star_thread": 5,
                "featurecount_thread": 10
}

#SAMPLES = [d for d in os.listdir("./") if os.path.isdir(f"./{d}") and d != ".snakemake"]
SAMPLES = ["NC1D","NC2D","NC3D","KI1D","KI2D","KI3D","S1D","S2D","S3D","T1D","T2D","T3D"]

rule all:
        input:
                expand("{sample_dir}/star/{sample}_Aligned.sortedByCoord.out.bam",
                                sample_dir = config["sample_dir"],
                                sample = SAMPLES)

rule fastqc:
        input:
                expand("{sample_dir}/{sample}/{sample}_{read}.fq.gz",
                                sample_dir = config["sample_dir"],
                                sample=SAMPLES,
                                read=["1","2"])

        output:
                directory("{sample_dir}/{result_dir}/fastqc")

        log:
                "{sample_dir}/{result_dir}/fastqc/fastqc_run.log"

        threads: config["fastqc_threads"]

        shell:
                """
                mkdir -p {output}
                touch {log}

        fastqc -o {output} -t {threads} {input} > {log} 2>&1
        """

rule multiqc:
        input:
                "{sample_dir}/{result_dir}/fastqc"
        output:
                directory("{sample_dir}/{result_dir}/multiqc")

        log:
                "{sample_dir}/{result_dir}/multiqc/multiqc_run.log"
        shell:
                """
                mkdir -p {output}
                touch {log}

                multiqc {input} --outdir {output} > {log} 2>&1
                """

rule fastp:
        input:
                in_r1fq = "{sample_dir}/{sample}/{sample}_1.fq.gz",
                in_r2fq = "{sample_dir}/{sample}/{sample}_2.fq.gz"
        output:
                out_r1fq = "{sample_dir}/fastp/{sample}_R1.trimmed.fq.gz",
                out_r2fq = "{sample_dir}/fastp/{sample}_R2.trimmed.fq.gz",
                out_json = "{sample_dir}/fastp/{sample}_fastp.json",
                out_html = "{sample_dir}/fastp/{sample}_fastp.html",
        log:
                "{sample_dir}/fastp/{sample}_fastp.log"
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/fastp"

                fastp -i {input.in_r1fq}  -I {input.in_r2fq} \
                                -o {output.out_r1fq} -O {output.out_r2fq} \
                                -j {output.out_json} -h {output.out_html} \
                                > {log} 2>&1
                """

rule star_align:
        input:
                genome_dir = config["GENOMEDIR"],
                in_r1fq = "{sample_dir}/fastp/{sample}_R1.trimmed.fq.gz",
                in_r2fq = "{sample_dir}/fastp/{sample}_R2.trimmed.fq.gz"
        output:
                outbam = "{sample_dir}/star/{sample}_Aligned.sortedByCoord.out.bam"
        log:
                "{sample_dir}/star/{sample}_runstar.log"
        params:
                out_prefix = "{sample_dir}/star/{sample}_"

        threads: config["star_thread"]

        shell:
                """
                mkdir -p "{wildcards.sample_dir}/star"

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
                inbams = "{sample_dir}/star/{sample}_Aligned.sortedByCoord.out.bam",
                inannogtf = config["aootgtf"]

        output:
                featurecounts_rawout = "{sample_dir}/featurecounts/featurecounts_rawout_counts.txt",
                filter_rawcounts = "{sample_dir}/featurecounts/featurecounts_rawcounts.csv"
        log:
                "{sample_dir}featurecounts/featurecounts_run.log"
        threads: config["featurecounts_thread"]
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/featurecounts"

                featureCounts -p \
                                --countReadPairs \
                                -t exon -g gene \
                                -a {input.inannogtf} \
                                --primary \
                                -T {threads} \
                                -o {output.featurecounts_rawout} \
                                {inbams}

                awk -F'\t' '{printf "%s,", $1; for(i=6; i<=NF; i++) printf "%s%s", $i, (i==NF?"\n":",")}' {output.featurecounts_rawout}  > {output.filter_rawcounts}
                """
