#!/bin/bash
#SBATCH --job-name=dupcaller_trim
#SBATCH --output=dupcaller_trim_%A_%a.log
#SBATCH --error=dupcaller_trim_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-89%10  # Process 100 files, 20 at a time

cd fastq

files=(*_1.fq.gz)
fq1=${files[$SLURM_ARRAY_TASK_ID - 1]}
sample_name=${fq1/_1*/}
fq2="${fq1/1.fq.gz/2.fq.gz}"

# Run DupCallerTrim.py
DupCallerTrim.py -i "$fq1" -i2 "$fq2" -p NNNNNNNN -o "$sample_name" &&

# Compress output FASTQ files
gzip "${sample_name}_1.fq" &&
gzip "${sample_name}_2.fq"
