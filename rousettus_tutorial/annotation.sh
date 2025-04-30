#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=annotation
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G

ENV_FILE=$1
. $ENV_FILE

mkdir $BLAST_HITS_DIR
mkdir $ACCESSION_GENE_MAPPING_DIR

module load singularity

# Extract candidate cds for blastn search
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $METADATA_DIR:/src/metadata \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -transcript_metadata /src/metadata/transcriptome_metadata.csv \
    -cds_metadata /src/metadata/cds_metadata.csv \
    -extraction_fields /src/app/example_extraction_fields/cds_for_blastn_extraction_fields.json \
    -cds_fasta /src/transcriptome_data/cds.fna

# Perform blastn search with candidate cds as queries
# and vertebrata core nt database for subject matches
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    $SINGULARITY_IMAGE \
    blastn -db /src/blastdb/vertebrata_core_nt \
    -query /src/transcriptome_data/cds.fna \
    -out /src/blast_hits/cds_hits.tsv \
    -evalue 0.0001 -num_threads 19 -mt_mode 1 \
    -outfmt "6 qseqid sacc qlen" -max_target_seqs 10

# Generate accession number to gene symbol mapping from blastn results
singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $APP_DIR:/src/app \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $ACCESSION_GENE_MAPPING_DIR:/src/accession_gene_mapping \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/make_accession_gene_symbol_mapping.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -datasets_mapping /src/accession_gene_mapping/datasets_mapping.csv

# Append accession numbers and genes to cds metadata,
# and pick "best" cds per genes, i.e. the "core" genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $ACCESSION_GENE_MAPPING_DIR:/src/accession_gene_mapping \
    --bind $METADATA_DIR:/src/metadata \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/assign_gene_annotations_to_cds.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -accession_gene_mapping /src/accession_gene_mapping/datasets_mapping.csv \
    -cds_metadata /src/metadata/cds_metadata.csv

# Extract core cds as final "core" transcriptome
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $METADATA_DIR:/src/metadata \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -transcript_metadata /src/metadata/transcriptome_metadata.csv \
    -cds_metadata /src/metadata/cds_metadata.csv \
    -extraction_fields /src/app/example_extraction_fields/core_gene_extraction_fields.json \
    -cds_fasta /src/transcriptome_data/cds_core.fna \
    -add_gene_name