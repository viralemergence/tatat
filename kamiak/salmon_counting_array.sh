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

mapfile -t SRA_NUMBERS < $SRA_LIST_FILE
SRA_NUMBER=${SRA_NUMBERS[$SLURM_ARRAY_TASK_ID]}

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR/salmon_index:/src/data/salmon_index \
    --bind $FASTQ_TRIMMED_DIR:/src/data/fastq_trimmed \
    --bind $SALMON_COUNTS_DIR:/src/data/salmon_counts \
    --bind $SALMON_COUNTS_COLLATED_DIR:/src/data/salmon_counts_collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_orchestration.py \
    -fastq_dir /src/data/fastq_trimmed \
    -sra $SRA_NUMBER -cpus 10 \
    -salmon_index /src/data/salmon_index/rousettus_core \
    -outdir /src/data/salmon_counts \
    -collated_dir /src/data/salmon_counts_collated