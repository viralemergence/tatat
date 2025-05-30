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

MEM_IN_GB=$((SLURM_MEM_PER_NODE / 1024))

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
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sample_metadata_extraction.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -array_index $SLURM_ARRAY_TASK_ID)

singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SRA_DOWNLOAD_DIR:/src/data/download_dir \
    --bind $SRA_COLLATE_DIR:/src/data/collate_dir \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/sra_read_download.py \
    -sra_number $SRA_NUMBER \
    -download_dir /src/data/download_dir \
    -collate_dir /src/data/collate_dir \
    -sqlite_db /src/sqlite_db/tatat.db

singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SRA_COLLATE_DIR:/src/data/fastq_dir \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/fastp_orchestration.py \
    -fastq_dir /src/data/fastq_dir \
    -outdir /src/data/trimmed \
    -sqlite_db /src/sqlite_db/tatat.db \
    -uid $SRA_NUMBER -r1_adapter $R1_ADAPTER -r2_adapter $R2_ADAPTER -cpus $SLURM_CPUS_PER_TASK

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
    -unique_identifier $SRA_NUMBER -cpus $SLURM_CPUS_PER_TASK -memory $MEM_IN_GB