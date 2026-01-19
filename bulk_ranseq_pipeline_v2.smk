import os

config = {
                "sample_dir": "/projects/bulkRNAseq/A10003198_rawdata",
                "software_dir":"/software/miniconda3/envs/rna-seq/bin",
                "log_dir": "logs",
                "fastqc_threads": 5,
                "GENOMEDIR": "/database/GCF_000001635.27_GRCm39",
                "annogtf": "/database/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf",
                "star_thread": 10,
                "featurecounts_thread": 5
}


SAMPLES = ["KO19SLOMP","KO25SLOMP","KO31SLOMP","KO51SLOMP","KO52SLOMP","KO55SLOMP","WT710SLOMP","WT726SLOMP","WT727SLOMP","WT737SLOMP","WT744SLOMP","WT745SLOMP"]

rule all:
        input:
                featurecounts_rawout = expand("{sample_dir}/5_featurecounts/featurecounts_rawout_counts.txt",sample_dir = config["sample_dir"]),
                filter_rawcounts = expand("{sample_dir}/5_featurecounts/featurecounts_rawcounts.csv",sample_dir = config["sample_dir"]),
                final_rawcounts = expand("{sample_dir}/5_featurecounts/featurecounts_rawcounts_final.csv",sample_dir = config["sample_dir"])
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
        params:
                fastqc_fullpath=os.path.join(config["software_dir"],"fastqc")

        threads: config["fastqc_threads"]

        shell:
                """
                mkdir -p {output}
                touch {log}

                {params.fastqc_fullpath} -o {output} -t {threads} {input} > {log} 2>&1
        """

rule multiqc:
        input:
                "{sample_dir}/1_fastqc"
        output:
                outdir = directory("{sample_dir}/2_multiqc"),
                flag_multiqc_done = "{sample_dir}/2_multiqc/multiqc_done.flag"

        log:
                "{sample_dir}/2_multiqc/multiqc_run.log"
        params:
                multiqc_fullpath=os.path.join(config["software_dir"],"multiqc")
        shell:
                """
                mkdir -p {output.outdir}
                touch {log}

                {params.multiqc_fullpath} {input} --outdir {output.outdir} > {log} 2>&1
                echo 'multiqc done!' > {output.flag_multiqc_done}
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
        params:
                fastp_fullpath=os.path.join(config["software_dir"],"fastp")
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/3_fastp"

                {params.fastp_fullpath} -i {input.in_r1fq}  -I {input.in_r2fq} \
                                -o {output.out_r1fq} -O {output.out_r2fq} \
                                --detect_adapter_for_pe \
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
                out_prefix = "{sample_dir}/4_star/{sample}_",
                star_fullpath=os.path.join(config["software_dir"],"STAR")
        threads: config["star_thread"]

        shell:
                """
                mkdir -p "{wildcards.sample_dir}/4_star"

                {params.star_fullpath} --runThreadN {threads} \
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
                filter_rawcounts = "{sample_dir}/5_featurecounts/featurecounts_rawcounts.csv",
                final_rawcounts = "{sample_dir}/5_featurecounts/featurecounts_rawcounts_final.csv"
        log:
                "{sample_dir}/5_featurecounts/featurecounts_run.log"
        params:
                featurecounts_fullpath=os.path.join(config["software_dir"],"featureCounts")
        threads: config["featurecounts_thread"]
        shell:
                """
                mkdir -p "{wildcards.sample_dir}/5_featurecounts"

                {params.featurecounts_fullpath} -p \
                                --countReadPairs \
                                -t exon -g gene \
                                -a {input.inannogtf} \
                                --primary \
                                -T {threads} \
                                -o {output.featurecounts_rawout} {input.inbams} > {log} 2>&1

                awk -F'\t' '{{printf "%s,", $1; for(i=6; i<=NF; i++) printf "%s%s", $i, (i==NF?"\\n":",")}}' {output.featurecounts_rawout}  > {output.filter_rawcounts}

                awk -F ',' '{{
                    if (NR == 1) {{
                        printf "%s,%s", $2, $3;
                        for (i = 4; i <= NF; i++) {{
                            split($i, parts, "/");
                            sample_name = parts[length(parts)];
                            gsub("_Aligned.*", "", sample_name);
                            printf "%s%s", (i == 3 ? "" : ","), sample_name;
                        }}
                        print "";
                    }} else {{
                        print $0;
                    }}
                }}' {output.filter_rawcounts}  > {output.final_rawcounts}

                """
