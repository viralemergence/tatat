#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=accession_gene_mapping
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $HOST_APP_DIR:/src/app \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $ACCESSION_GENE_MAPPING_DIR:/src/accession_gene_mapping \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/make_accession_gene_symbol_mapping.py \
    -diamond_results /src/blast_hits/cds_hits.tsv \
    -datasets_mapping /src/accession_gene_mapping/datasets_mapping.csv \
    -accepted_accession_numbers /src/accession_gene_mapping/accepted_accession_numbers.csv \
    -accession_gene_id_mapping /src/accession_gene_mapping/accession_gene_id_mapping.csv