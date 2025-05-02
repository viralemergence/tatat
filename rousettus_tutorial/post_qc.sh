#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=post_qc
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G

ENV_FILE=$1
. $ENV_FILE

mkdir $BUSCO_DATA_DIR

module load singularity

# To calculate intersection of TATAT core genes with NCBI genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/intersect_core_genes.py \
    -cds_metadata /src/data/metadata/cds_metadata.csv \
    -ncbi_genes_path /src/data/ncbi_gene_comparison/ncbi_gene_results_9407.txt

# To examine TATAT core genes in relation to human gene expression
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/process_human_tissue_expression.py \
    -ncbi_genes_path /src/data/ncbi_gene_comparison/ncbi_gene_results_9407.txt \
    -tissue_expression_path /src/data/ncbi_gene_comparison/human_tissue_expression.tsv \
    -missing_genes_path /src/data/ncbi_gene_comparison/missing_ncbi_genes.txt

# Extract TATAT core aa sequences for pairwise alignments to NCBI aa seqs
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
    -aa_fasta /src/transcriptome_data/aa_core.faa -add_gene_name

# Perform pairwise alignments of TATAT core genes to NCBI genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $DATA_DIR/ncbi_datasets_download/ncbi_dataset/data/GCF_014176215.1:/src/data/9407_data \
    --bind $DATA_DIR/ncbi_gene_comparison:/src/data/ncbi_gene_comparison \
    $SINGULARITY_IMAGE \
    python3 /src/app/post_annotation_qc/pairwise_align_core_genes.py \
    -tatat_cds_fasta /src/transcriptome_data/cds_core.fna \
    -tatat_aa_fasta /src/transcriptome_data/aa_core.faa \
    -ncbi_cds_fasta /src/data/9407_data/cds_from_genomic.fna \
    -ncbi_aa_fasta /src/data/9407_data/protein.faa \
    -outdir /src/data/ncbi_gene_comparison -cpus 10

# Index core genes with salmon
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $DATA_DIR/salmon_index:/src/salmon_index \
    $SINGULARITY_IMAGE \
    salmon index -t /src/transcriptome_data/cds_core.fna \
    -i /src/salmon_index/rousettus_core -p 10

# Calculate salmon counts as array (depending on sample count)
SALMON_ARRAY_SCRIPT=$ROUSETTUS_TUTORIAL_DIR/post_qc_salmon_array.sh
sbatch -W $SALMON_ARRAY_SCRIPT $ENV_FILE

# Collate salmon counts
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    --bind $SALMON_COUNTS_COLLATED_DIR:/src/salmon_counts_collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_count_collation.py \
    -counts_dir /src/salmon_counts_collated \
    -outdir /src/data

# MDS salmon counts
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_count_mds.py \
    -counts /src/data/collated_salmon_counts.csv \
    -metadata /src/data/metadata/sample_metadata.csv \
    -outdir /src/data

# Extract TATAT core aa sequences for BUSCO
# (adding gene names causes BUSCO to fail, so they are ignored here)
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
    -aa_fasta /src/transcriptome_data/aa_core_busco.faa

# Run BUSCO on aa_core_busco.faa
singularity exec \
    --pwd /src/data/busco \
    --no-home \
    --bind /etc:/etc \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $BUSCO_DATA_DIR:/src/data/busco \
    $SINGULARITY_IMAGE \
    busco -i /src/transcriptome_data/aa_core_busco.faa \
    -m proteins -l vertebrata_odb12 -c 10 \
    --opt-out-run-stats