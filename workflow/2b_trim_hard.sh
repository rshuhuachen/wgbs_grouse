#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=trim
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/trim_hard.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/trim_hard.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/2b_trim_hard.sh

cd /vol/cluster-data/rchen/wgbs

trim_galore --paired data/E197_R1.fq.gz data/E197_R2.fq.gz --quality 20 --fastqc --clip_R1 3 --clip_R2 3 --gzip --cores 8 -o ./data/processed/hardtrim/ 2>&1 | tee log/trim.log
