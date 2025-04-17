#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=sra_processing_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-19:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G

ENV_FILE=$1
. $ENV_FILE

mkdir $SRA_DOWNLOAD_DIR
mkdir $SRA_COLLATE_DIR
mkdir $FASTQ_TRIMMED_DIR
mkdir $RNASPADES_ASSEMBLY_DIR
mkdir $RNASPADES_COLLATED_ASSEMBLY_DIR

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
    python3 -u /src/app/assembly/sra_read_download.py \
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
    python3 -u /src/app/assembly/fastp_orchestration.py \
    -i /src/data/fastq_dir \
    -o /src/data/trimmed \
    -u $SRA_NUMBER -r1a $R1_ADAPTER -r2a $R2_ADAPTER --cpus 10