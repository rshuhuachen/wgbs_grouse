#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=polya
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/polya.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/polya.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/3_polyA.sh

cd /vol/cluster-data/rchen/wgbs

cutadapt --poly-a -o data/processed/trim/polya/E197_R1.fastq -p data/processed/trim/polya/E197_R2.fastq ./data/processed/trim/E197_R1_val_1.fq.gz ./data/processed/trim/E197_R2_val_2.fq.gz

#trim_galore --polyA ./data/processed/trim/E197_R2_val_2.fq.gz -o data/processed/trim/E197_R2_val_2_trimmed.fq.gz 2>&1 | tee log/polya.log
#zcat ./data/processed/trim/E197_R2_val_2_trimmed.fq.gz | grep -A 3 PolyA | grep -v ^-- > ./data/processed/trim/E197_R2_val_2_PolyA_trimmed.fastq
