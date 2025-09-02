#!/usr/bin/bash

dnbc4tools='/software/dnbc4tools2.1.3/dnbc4tools'

samplepath=${1}
genomedir=${2}
threads=${3}

fq1=$(ls ${samplepath}/*1.fq.gz | paste -sd ',')
fq2=$(ls ${samplepath}/*2.fq.gz | paste -sd ',')

oligopath=$(echo "${curentpath%_cDNA*}_oligo")

samplename=$(basename ${samplepath} | sed 's/_cDNA//g')

oligo1="${oligopath}/*oligo_1.fq.gz"
oligo2="${oligopath}/*oligo_2.fq.gz"

$dnbc4tools rna run \
        --cDNAfastq1 ${fq1} \
        --cDNAfastq2 ${fq2} \
        --oligofastq1 ${oligo1} \
        --oligofastq2 ${oligo2} \
        --name ${samplename} \
        --genomeDir ${genomedir} \
        --threads ${threads}
