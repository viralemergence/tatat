#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=sra_paa_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-19:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

ENV_FILE=$1
. $ENV_FILE

mkdir $SRA_DOWNLOAD_DIR
mkdir $SRA_COLLATE_DIR
mkdir $FASTQ_TRIMMED_DIR
mkdir $RNASPADES_ASSEMBLY_DIR
mkdir $RNASPADES_COLLATED_ASSEMBLY_DIR

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
    --bind $APP_DIR:/src/app \
    --bind $SRA_COLLATE_DIR:/src/data/fastq_dir \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/fastp_orchestration.py \
    -i /src/data/fastq_dir \
    -o /src/data/trimmed \
    -u $SRA_NUMBER -r1a $R1_ADAPTER -r2a $R2_ADAPTER --cpus 10

singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    --bind $RNASPADES_ASSEMBLY_DIR:/src/data/assembly \
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/rnaspades_orchestration.py \
    -fastq_dir /src/data/trimmed \
    -assembly_dir /src/data/assembly \
    -collated_dir /src/data/collated \
    -unique_identifier $SRA_NUMBER -cpus 10 -memory 49