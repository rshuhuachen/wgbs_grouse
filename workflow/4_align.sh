#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=alignment
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/alignment.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/alignment.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/align.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd),/vol/cluster-data/rchen/geneticload/snpEff/data/genomes/ \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/bismark_0_24.sif \
  bismark -p 8 --parallel 4 --un --ambiguous --genome data/ref/ --ignore 3 --ignore_r2 3 \
  -1 data/processed/trim/E197_R1_val_1.fq.gz -2 data/processed/trim/E197_R2_val_2.fq.gz \
  -o data/processed/alignment --rg_tag 
