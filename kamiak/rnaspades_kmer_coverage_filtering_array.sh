#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=kmer_cov_filter_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-19:1%20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
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
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    --bind $KMER_COVERAGE_FILTERED_DIR:/src/data/filtered \
    --bind $HOST_DATA_DIR/kmer_qc:/src/data/kmer_qc \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/kmer_coverage_filtering.py \
    -assembly_path /src/data/collated/$SRA_NUMBER.fasta \
    -coverage_cutoff 1.5 \
    -outdir /src/data/filtered \
    -qc_dir /src/data/kmer_qc