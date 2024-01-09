#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=10000
#SBATCH --job-name=methcall
#SBATCH --nodes=2
#SBATCH --error=/vol/cluster-data/rchen/wgbs/log/methcall.err
#SBATCH --output=/vol/cluster-data/rchen/wgbs/log/methcall.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgbs/workflow/methcall.sh

cd /vol/cluster-data/rchen/wgbs
WGBS_SCRATCH=$(pwd)/scratch/

#enter container
/vol/apptainer/bin/apptainer run --cleanenv \
  --fakeroot --overlay ${WGBS_SCRATCH} \
  --bind ${WGBS_SCRATCH}/tmp,$(pwd),/vol/cluster-data/rchen/geneticload/snpEff/data/genomes/ \
  --env PYTHONNOUSERSITE=1 \
  /vol/cluster-data/rchen/epi/epigbs2/singularity/bismark_0_24.sif \
  bismark_methylation_extractor --multicore 4 --ignore 4 --ignore_r2 4 -p --CX --no_overlap --scaffolds --report \
  --cutoff 10 \
  --bedGraph --cytosine_report --comprehensive --gzip --merge_non_CpG --genome_folder data/ref/ data/processed/deduplicated/E197_R1_val_1_bismark_bt2_pe.deduplicated.bam \
  -o data/processed/methcall
