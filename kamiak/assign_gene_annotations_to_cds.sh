#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=gene_annotations
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
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $HOST_APP_DIR:/src/app \
    --bind $DIAMOND_HITS_DIR:/src/diamond_hits \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/assign_gene_annotations_to_cds.py \
    -diamond_results /src/diamond_hits/aa_hits.tsv \
    -accesion_numbers_out /src/diamond_hits/accession_numbers.csv \
    -cds_metadata /src/data/transcripts/cds_metadata.csv