#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=ncrna
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G

ENV_FILE=$1
. $ENV_FILE

mkdir $NCRNA_DIR
mkdir $BLAST_HITS_DIR

module load singularity

# Generate empty ncrna table in sqlite db
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_ncrna_table

# Populate ncrna table with initial filtered transcripts ids
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/ncrna_initial_filtering.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcripts_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -transcriptome rousettus

# Remove ncrna candidates that map to core cds (via CD-HIT-EST-2D)
# and then remove ncrna sequence duplicates/fragments (via CD-HIT-EST)
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $NCRNA_DIR:/src/ncrna \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/cd_hit_orchestration.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcriptome rousettus \
    -transcripts_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -ncrna /src/ncrna \
    -cds_fasta /src/transcriptome_data/rousettus_cds_core.fna \
    -cpus $SLURM_CPUS_PER_TASK -memory $SLURM_MEM_PER_NODE

# Extract ncRNA for BLAST search
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $NCRNA_DIR:/src/ncrna \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/blast_ncrna_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -ncrna_fasta /src/ncrna/blast_ncrna.fna

# Perform blastn search with candidate ncrna as queries
# and vertebrata core nt database for subject matches
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    --bind $NCRNA_DIR:/src/ncrna \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    $SINGULARITY_IMAGE \
    blastn -db /src/blastdb/vertebrata_core_nt \
    -query /src/ncrna/blast_ncrna.fna \
    -out /src/blast_hits/ncrna_hits.tsv \
    -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -mt_mode 1 \
    -outfmt "6 qseqid sacc qlen" -max_target_seqs 1

# Generate sqlite db "nc_accession_numbers" table
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_nc_acc_num_table

# Generate accession number to gene symbol mapping from blastn results
# and add to sqlite db as "nc_accession_numbers" table
singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $APP_DIR:/src/app \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/make_accession_gene_symbol_mapping.py \
    -blast_results /src/blast_hits/ncrna_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db \
    -table_name nc_accession_numbers -rna_type non_coding

# Append accession numbers and genes to ncrna metadata sqlite table,
# and pick "best" sequence per gene, i.e. the "core" ncRNA
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/assign_gene_annotations_to_ncrna.py \
    -blast_results /src/blast_hits/ncrna_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db -transcriptome rousettus

# Extract core ncRNA
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/ncrna/core_ncrna_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -ncrna_fasta /src/transcriptome_data/rousettus_ncrna_core.fna \
    -add_gene_name -transcriptome rousettus