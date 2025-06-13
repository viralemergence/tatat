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

module load singularity

# Extract candidate cds for blastn search
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/cds_for_blastn_sql_queries.json \
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
    -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -mt_mode 1 \
    -outfmt "6 qseqid sacc qlen" -max_target_seqs 10

# Generate sqlite db "accession_numbers" table
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_acc_num_table

# Generate accession number to gene symbol mapping from blastn results
# and add to sqlite db as "accession_numbers" table
singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/make_accession_gene_symbol_mapping.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db \
    -table_name accession_numbers -rna_type coding

# Append accession numbers and genes to cds metadata sqlite table,
# and pick "best" cds per gene, i.e. the "core" genes
singularity exec \
    --pwd /src \
    --no-home \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/assign_gene_annotations_to_cds.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db -transcriptome rousettus

# Extract core cds as final "core" transcriptome
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/core_gene_sql_queries.json \
    -cds_fasta /src/transcriptome_data/rousettus_cds_core.fna \
    -add_gene_name -transcriptome rousettus