#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=annotate_evigene
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G

ENV_FILE=$1
. $ENV_FILE

module load singularity

singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    --bind $DIAMOND_HITS_DIR:/src/diamond_hits \
    $SINGULARITY_IMAGE \
    diamond blastp -d /src/blastdb/vertebrata \
    -q /src/data/transcripts/aa.fasta \
    -o /src/diamond_hits/aa_hits.tsv \
    --outfmt 6 qseqid sseqid qlen \
    -p 20 -b 8 -c 1 -k 5

singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $HOST_APP_DIR:/src/app \
    --bind $DIAMOND_HITS_DIR:/src/diamond_hits \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/assign_gene_annotations_to_cds.py \
    -diamond_results /src/diamond_hits/aa_hits.tsv \
    -accesion_numbers_out /src/diamond_hits/accession_numbers.csv \
    -cds_metadata /src/data/transcripts/cds_metadata.csv

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/intersect_core_genes.py \
    -cds_metadata /src/data/assemblies_and_metadata/cds_metadata.csv \
    -ncbi_genes_path /src/data/ncbi_gene_comparison/ncbi_gene_results_9407.txt

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/process_human_tissue_expression.py \
    -ncbi_genes_path /src/data/ncbi_gene_comparison/ncbi_gene_results_9407.txt \
    -tissue_expression_path /src/data/ncbi_gene_comparison/human_tissue_expression.tsv \
    -missing_genes_path /src/data/ncbi_gene_comparison/missing_ncbi_genes.txt

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/data/transcripts/raw_transcriptome.fasta \
    -transcript_metadata /src/data/transcripts/transcriptome_metadata.csv \
    -cds_metadata /src/data/transcripts/cds_metadata.csv \
    -extraction_fields /src/app/example_extraction_fields/core_gene_extraction_fields.json \
    -aa_fasta /src/data/transcripts/aa_core.fasta -add_gene_name

singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    --bind $HOST_DATA_DIR/ncbi_datasets_download/ncbi_dataset/data/GCF_014176215.1:/src/data/9407_data \
    --bind $HOST_DATA_DIR/ncbi_gene_comparison:/src/data/ncbi_gene_comparison \
    $SINGULARITY_IMAGE \
    python3 /src/app/post_annotation_qc/pairwise_align_core_genes.py \
    -tatat_aa_fasta /src/data/transcripts/aa_core.fasta \
    -ncbi_cds_fasta /src/data/9407_data/cds_from_genomic.fna \
    -ncbi_aa_fasta /src/data/9407_data/protein.faa