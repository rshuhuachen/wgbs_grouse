#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=trim
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/trim.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/trim.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/2_trim.sh

cd /vol/cluster-data/rchen/wgbs

trim_galore --paired data/E197_R1.fq.gz data/E197_R2.fq.gz --quality 20 --fastqc --illumina --stringency 3 --gzip --length 20 --trim-n --cores 8 -o ./data/processed/trim/ 2>&1 | tee log/trim.log
