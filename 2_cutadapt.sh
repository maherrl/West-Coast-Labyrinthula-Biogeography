#!/bin/bash

#SBATCH --job-name=cutadapt
#SBATCH --ntasks=10
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=cut.out
#SBATCH --error=cut.err

module load cutadapt/4.1

source activate cutadapt-4.1

in="/share/cdfwlab/maherr/18S_sub/data"

cd ${in}

mkdir cutadapt

for f in $(ls -1 *_L001_R1_001.fastq | sed 's/\_L001_R1_001.fastq//'); do cutadapt -g CACCTCTGACATACTCATACG -G CAATCCTGATACGGGGAGG --discard-untrimmed --no-indels --cores=10 -o cutadapt/${f}_1.fq -p cutadapt/${f}_2.fq ./${f}_L001_R1_001.fastq ./${f}_L001_R2_001.fastq; done
