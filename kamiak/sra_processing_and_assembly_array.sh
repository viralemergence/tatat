#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=sra_paa_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --array=0-19:1%5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

ENV_FILE=$1
. $ENV_FILE

mapfile -t SRA_NUMBERS < $SRA_LIST_FILE
SRA_NUMBER=${SRA_NUMBERS[$SLURM_ARRAY_TASK_ID]}

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $SRA_DOWNLOAD_DIR:/src/data/download_dir \
    --bind $SRA_COLLATE_DIR:/src/data/collate_dir \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sra_read_download.py \
    -sra_number $SRA_NUMBER \
    -download_dir /src/data/download_dir \
    -collate_dir /src/data/collate_dir

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $SRA_COLLATE_DIR:/src/data/fastq_dir \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/fastp_orchestration.py \
    -i /src/data/fastq_dir \
    -o /src/data/trimmed \
    -u $SRA_NUMBER -r1a $R1_ADAPTER -r2a $R2_ADAPTER --cpus 10

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    --bind $RNASPADES_ASSEMBLY_DIR:/src/data/assembly \
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/rnaspades_orchestration.py \
    -fastq_dir /src/data/trimmed \
    -assembly_dir /src/data/assembly \
    -collated_dir /src/data/collated \
    -unique_identifier $SRA_NUMBER -cpus 10 -memory 49