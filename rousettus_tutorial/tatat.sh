#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=tatat
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

ENV_FILE=$1
. $ENV_FILE

mkdir $DATA_DIR
mkdir $SQLITE_DB_DIR
mkdir $TRANSCRIPTOME_DATA_DIR
mkdir $SCRATCH_DIR

module load singularity

# Generate sqlite db with sample metadata table
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -sample_metadata /src/sqlite_db/sample_metadata.csv \
    -create_transcripts_table \
    -create_cds_table \
    -create_acc_num_table

# Run assembly phase of TATAT
ASSEMBLY_ARRAY_SCRIPT=$ROUSETTUS_TUTORIAL_DIR/tatat_main/assembly_array.sh
sbatch -W $ASSEMBLY_ARRAY_SCRIPT $ENV_FILE

# Run thinning phase of TATAT
THINNING_SCRIPT=$ROUSETTUS_TUTORIAL_DIR/tatat_main/thinning.sh
sbatch -W $THINNING_SCRIPT $ENV_FILE

# Run annotation phase of TATAT
ANNOTATION_SCRIPT=$ROUSETTUS_TUTORIAL_DIR/tatat_main/annotation.sh
sbatch -W $ANNOTATION_SCRIPT $ENV_FILE