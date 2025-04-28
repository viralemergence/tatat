#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=make_vertebrata_core_nt_db
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G

ENV_FILE=$1
. $ENV_FILE

module load singularity

# To extract vertebrata sequences from core nt database
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $SCRATCH_DIR/core_nt:/src/blastdb \
    --bind $SCRATCH_DIR/blast:/src/tatat_blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    extract-taxon-nucleotide-sequences \
    7742 /src/tatat_blastdb/7742.fna

# To build vertebrata blast database from nucleotide sequences
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $SCRATCH_DIR/blast:/src/blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    make_nucleotide_blast_db \
    /src/blastdb/7742.fna "vertebrata_core_nt" 7742 /src/blastdb/vertebrata_core_nt