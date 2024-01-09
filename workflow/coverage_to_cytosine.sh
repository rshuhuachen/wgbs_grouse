#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=10000
#SBATCH --job-name=report
#SBATCH --nodes=2
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/cov_to_cytosine.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/cov_to_cytosine.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/coverage_to_cytosine.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd),/vol/cluster-data/rchen/geneticload/snpEff/data/genomes/ \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/bismark_0_24.sif \
  coverage2cytosine data/processed/methcall/E197_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
  --genome_folder data/ref/ --CX_context --coverage_threshold 10 --gzip \
  --dir data/processed/methcall -o methcall_report
