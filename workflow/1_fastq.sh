#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=fastqc
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/fastqc.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/fastqc.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/1_fastq.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd) \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/fastqc_0.11.9--hdfd78af_1.sif \
  fastqc -t 2 --quiet --outdir /vol/cluster-data/rchen/wgbs/data/processed/fastqc/ data/E197_R1.fq.gz data/E197_R2.fq.gz
