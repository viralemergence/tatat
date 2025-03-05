#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=assembling_merging
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $KMER_COVERAGE_FILTERED_DIR:/src/data/filtered \
    --bind $RNASPADES_ASSEMBLIES_MERGED:/src/data/merged \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/de_novo_assembly_merging.py \
    -assembly_fasta_dir /src/data/filtered \
    -outdir /src/data/merged \
    -cpus 20 -mem 40000 -concat_size 6