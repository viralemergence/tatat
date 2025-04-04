#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=ncbi_nr_pull
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --env LC_ALL=C,NCBI_API_KEY=$NCBI_API_KEY \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_BLASTDB_DIR:/src/blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh pull-nr-database