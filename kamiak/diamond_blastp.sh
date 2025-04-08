#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=diamond_blastp
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G

ENV_FILE=$1
. $ENV_FILE

module load singularity

singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    --bind $DIAMOND_HITS_DIR:/src/diamond_hits \
    $SINGULARITY_IMAGE \
    diamond blastp -d /src/blastdb/vertebrata \
    -q /src/data/transcripts/aa.fasta \
    -o /src/diamond_hits/aa_hits.tsv \
    --outfmt 6 qseqid sseqid qlen \
    -p 20 -b 8 -c 1 -k 1