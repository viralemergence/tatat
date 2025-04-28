#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=blastn
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=4-00:00:00
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
    --bind $SCRATCH_DIR/blast:/src/blastdb \
    --bind $HOST_DATA_DIR/assemblies_and_metadata:/src/data/transcripts \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    $SINGULARITY_IMAGE \
    blastn -db /src/blastdb/vertebrata_core_nt \
    -query /src/data/transcripts/cds.fasta \
    -out /src/blast_hits/cds_hits.tsv \
    -evalue 0.0001 -num_threads 19 -mt_mode 1 \
    -outfmt "6 qseqid sacc qlen" -max_target_seqs 5