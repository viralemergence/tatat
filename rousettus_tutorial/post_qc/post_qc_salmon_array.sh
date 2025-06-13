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
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sample_metadata_extraction.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -array_index $SLURM_ARRAY_TASK_ID \
    -return_uid)

singularity exec \
    --pwd /src \
    --no-home \
    --bind $POST_QC_DIR/salmon_index:/src/salmon_index \
    --bind $FASTQ_TRIMMED_DIR:/src/fastq_trimmed \
    --bind $SALMON_COUNTS_DIR:/src/salmon_counts \
    --bind $SALMON_COUNTS_COLLATED_DIR:/src/salmon_counts_collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_orchestration.py \
    -fastq_dir /src/fastq_trimmed \
    -sra $SRA_NUMBER -cpus $SLURM_CPUS_PER_TASK \
    -salmon_index /src/salmon_index/rousettus_core \
    -outdir /src/salmon_counts \
    -collated_dir /src/salmon_counts_collated