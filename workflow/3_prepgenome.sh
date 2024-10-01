#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=14000
#SBATCH --job-name=prepgenome
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/prepgenome.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/prepgenome.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/prepgenome.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd),/vol/cluster-data/rchen/geneticload/snpEff/data/genomes/ \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/bismark_0_24.sif \
  bismark_genome_preparation /vol/cluster-data/rchen/geneticload/snpEff/data/genomes/

 #mv /vol/cluster-data/rchen/geneticload/snpEff/data/genomes/Bisulfite_Genome/ ./data/ref

