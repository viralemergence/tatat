#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=thinning
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G

ENV_FILE=$1
. $ENV_FILE

mkdir $EVIGENE_OUTPUT_DIR

module load singularity

# Generate sqlite db tables
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_transcripts_table \
    -create_cds_table

# To merge assemblies into single file
# and update corresponding sqlite metadata table "transcripts"
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/merge_fastas_and_set_metadata.py \
    -assembly_fasta_dir /src/data/collated \
    -merged_path /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db

# Use evigene to calculate candidate cds regions in assemblies,
# classify them as coding, noncoding, etc.,
# and update the sqlite tables with this information
singularity exec \
    --pwd /src \
    --no-home \
    --env LC_ALL=C \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $EVIGENE_OUTPUT_DIR:/src/evigene_output \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/evigene_orchestration.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -outdir /src/evigene_output \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcriptome rousettus -prefix_column sample_uid \
    -run_evigene -cpus $SLURM_CPUS_PER_TASK -mem $SLURM_MEM_PER_NODE -phetero 2 -minaa 99 \
    -run_transcript_metadata_appender -run_cds_and_metadata

singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/evigene_orchestration.py \
    -assembly_fasta holder -outdir holder \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcriptome holder -prefix_column holder \
    -update_transcript_cds_ids