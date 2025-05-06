#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=prep_sqlite_db
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

# Generate sqlite db with sample metadata table
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite3_prep.py \
    -sample_metadata /src/sqlite_db/sample_metadata.csv \
    -sqlite_db_dir /src/sqlite_db