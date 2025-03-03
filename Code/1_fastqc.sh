#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --ntasks=10
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=qc.out
#SBATCH --error=qc.err

module load fastqc/0.12.1
module load multiqc/1.10.1

in="/share/cdfwlab/maherr/18S_sub/data/cutadapt"

cd ${in}

mkdir fastqc

fastqc *.fq -o ./fastqc -t 10

cd fastqc

mkdir multiqc

multiqc . --dirs -o multiqc --interactive
