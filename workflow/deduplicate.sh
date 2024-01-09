#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=deduplicate
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/deduplicate.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/deduplicate.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/deduplicate.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd) \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/bismark_0_24.sif \
  deduplicate_bismark E197_R1_val_1_bismark_bt2_pe.bam 
