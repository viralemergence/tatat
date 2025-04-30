#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=salmon_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-19:1%20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G

ENV_FILE=$1
. $ENV_FILE

mkdir $SALMON_COUNTS_DIR
mkdir $SALMON_COUNTS_COLLATED_DIR

module load singularity

SRA_NUMBER=$(singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $METADATA_DIR:/src/metadata \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sample_metadata_extraction.py \
    -sample_metadata /src/metadata/sample_metadata.csv \
    -array_index $SLURM_ARRAY_TASK_ID -encoding utf-8-sig)

singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR/salmon_index:/src/salmon_index \
    --bind $FASTQ_TRIMMED_DIR:/src/fastq_trimmed \
    --bind $SALMON_COUNTS_DIR:/src/salmon_counts \
    --bind $SALMON_COUNTS_COLLATED_DIR:/src/salmon_counts_collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_orchestration.py \
    -fastq_dir /src/fastq_trimmed \
    -sra $SRA_NUMBER -cpus 10 \
    -salmon_index /src/salmon_index/rousettus_core \
    -outdir /src/salmon_counts \
    -collated_dir /src/salmon_counts_collated