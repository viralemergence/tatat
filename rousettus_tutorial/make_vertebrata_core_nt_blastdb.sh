#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=make_vertebrata_core_nt_blastdb
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=15G

ENV_FILE=$1
. $ENV_FILE

mkdir $CORE_NT_BLASTDB_DIR

module load singularity

# To download core_nt blast database
singularity exec \
    --pwd /src \
    --no-home \
    --env LC_ALL=C,NCBI_API_KEY=$NCBI_API_KEY \
    --bind $APP_DIR:/src/app \
    --bind $CORE_NT_BLASTDB_DIR:/src/blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh pull-core-nt-database

# To extract vertebrata nucleotide sequences from core_nt blast database
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $CORE_NT_BLASTDB_DIR:/src/blastdb \
    --bind $TATAT_BLASTDB_DIR:/src/tatat_blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    extract-taxon-nucleotide-sequences \
    7742 /src/tatat_blastdb/7742.fna

# To build vertebrata blast database from nucleotide sequences
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    make_nucleotide_blast_db \
    /src/blastdb/7742.fna "vertebrata_core_nt" 7742 /src/blastdb/vertebrata_core_nt

# To remove nucleotide sequences to save space (optional)
rm $TATAT_BLASTDB_DIR/7742.fna