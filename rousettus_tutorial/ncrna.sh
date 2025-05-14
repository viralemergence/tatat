#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=ncrna
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

module load singularity

# Generate empty ncrna table in sqlite db
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_ncrna_table

# Populate ncrna table with initial filtered transcripts ids
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/ncrna_initial_filtering.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcripts_fasta /src/transcriptome_data/raw_transcriptome.fna