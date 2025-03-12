#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=evigene
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=150G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --env LC_ALL=C \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    --bind $EVIGENE_OUTPUT_DIR:/src/data/evigene_output \
    $SINGULARITY_IMAGE \
    bash /src/app/thinning/evigene_orchestration.sh \
    run_evigene /src/data/evigene_output \
    /src/data/transcripts/raw_transcriptome.fasta \
    /src/data/evigene_output/raw_transcriptome.fasta \
    15 149000