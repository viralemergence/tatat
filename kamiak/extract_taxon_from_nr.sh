#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=extract_taxon_from_nr
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
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_BLASTDB_DIR:/src/blastdb \
    --bind $PROTEIN_ACCESSION_DIR:/src/protein_accessions \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    extract-taxon-protein-sequences \
    7742 /src/protein_accessions/7742.fsa