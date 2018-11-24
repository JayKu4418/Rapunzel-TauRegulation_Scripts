#!/bin/bash

fastaFolder="/data1/MAPT_HCC1500_MiSeq_November2018"
dataFolder="/home/jkumar12/Projects/TauRegulation/data"
tmpFolder="/home/jkumar12/Projects/TauRegulation/tmp"

shapemapper --target ${dataFolder}/MAPT_primerAmplified.fa --out ${tmpFolder}/MAPT_amplicon_HCC1500_Take2 --modified --R1 ${fastaFolder}/LLPooledPlus_R1_001.fastq.gz --R2 ${fastaFolder}/LLPooledPlus_R2_001.fastq.gz --untreated --R1 ${fastaFolder}/LLPooledMinus_R1_001.fastq.gz --R2 ${fastaFolder}/LLPooledMinus_R2_001.fastq.gz --amplicon --primers ${dataFolder}/MAPT_primers.fa --overwrite --min-depth 100 --nproc 8 --output-processed-reads