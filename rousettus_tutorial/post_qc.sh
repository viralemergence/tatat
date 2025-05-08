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

mkdir $POST_QC_DIR
mkdir $POST_QC_DIR/ncbi
mkdir $POST_QC_DIR/gene_comparison
mkdir $POST_QC_DIR/salmon_index
mkdir $BUSCO_DATA_DIR

module load singularity

# Download Rousettus cds and aa sequences for comparisons to TATAT output
singularity exec \
    --pwd /src \
    --no-home \
    --bind /etc:/etc \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind $POST_QC_DIR/ncbi:/src/ncbi \
    $SINGULARITY_IMAGE \
    datasets download genome taxon 9407 --include cds,protein \
    --assembly-source "RefSeq" --assembly-version "latest" \
    --filename /src/ncbi/9407_datasets.zip

# Unzip NCBI download
singularity exec \
    --pwd /src \
    --no-home \
    --bind $POST_QC_DIR/ncbi:/src/ncbi \
    $SINGULARITY_IMAGE \
    unzip /src/ncbi/9407_datasets.zip -d /src/ncbi/

# Extract NCBI genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $POST_QC_DIR/ncbi:/src/ncbi \
    --bind $POST_QC_DIR/ncbi/ncbi_dataset/data/GCF_014176215.1:/src/cds \
    $SINGULARITY_IMAGE \
    python3 /src/app/post_annotation_qc/extract_ncbi_genes.py \
    -ncbi_cds_fasta /src/cds/cds_from_genomic.fna \
    -outpath /src/ncbi/9407_ncbi_genes.txt

# To calculate intersection of TATAT core genes with NCBI genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $POST_QC_DIR/ncbi:/src/ncbi \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/intersect_core_genes.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -ncbi_genes_path /src/ncbi/9407_ncbi_genes.txt \
    -outdir /src/gene_comparison

# Download human expression data
singularity exec \
    --pwd /src \
    --no-home \
    --bind /etc:/etc \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    wget https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz \
    -O /src/gene_comparison/human_tissue_expression.tsv.gz

# Decompress human expression data
singularity exec \
    --pwd /src \
    --no-home \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    gzip -d /src/gene_comparison/human_tissue_expression.tsv.gz

# Remove human expression data extraneous header
singularity exec \
    --pwd /src \
    --no-home \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    sed -i "1,2d" /src/gene_comparison/human_tissue_expression.tsv

# To examine TATAT core genes in relation to human gene expression
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $POST_QC_DIR/ncbi:/src/ncbi \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/process_human_tissue_expression.py \
    -ncbi_genes_path /src/ncbi/9407_ncbi_genes.txt \
    -tissue_expression_path /src/gene_comparison/human_tissue_expression.tsv \
    -missing_genes_path /src/gene_comparison/missing_ncbi_genes.txt \
    -outdir /src/gene_comparison

# Extract TATAT core cds sequences for pairwise alignments to NCBI cds seqs
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/core_gene_sql_queries.json \
    -cds_fasta /src/transcriptome_data/cds_core.fna \
    -add_gene_name

# Perform pairwise alignments of TATAT core genes to NCBI genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $POST_QC_DIR/ncbi/ncbi_dataset/data/GCF_014176215.1:/src/cds \
    --bind $POST_QC_DIR/gene_comparison:/src/gene_comparison \
    $SINGULARITY_IMAGE \
    python3 /src/app/post_annotation_qc/pairwise_align_core_genes.py \
    -tatat_cds_fasta /src/transcriptome_data/cds_core.fna \
    -ncbi_cds_fasta /src/cds/cds_from_genomic.fna \
    -outdir /src/gene_comparison -cpus 10

# Index core genes with salmon
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $POST_QC_DIR/salmon_index:/src/salmon_index \
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
    --bind $POST_QC_DIR:/src/data \
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
    --bind $POST_QC_DIR:/src/data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/post_annotation_qc/salmon_count_mds.py \
    -counts /src/data/collated_salmon_counts.csv \
    -sqlite_db /src/sqlite_db/tatat.db \
    -outdir /src/data

# Extract TATAT core aa sequences for BUSCO
# (adding gene names causes BUSCO to fail, so they are ignored here)
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/core_gene_sql_queries.json \
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