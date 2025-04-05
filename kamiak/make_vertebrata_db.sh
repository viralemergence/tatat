#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=make_vertebrata_db
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

# To extract vertebrata sequences from nr database (~20 min on HPC)
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $NR_BLASTDB_DIR:/src/blastdb \
    --bind $TATAT_BLASTDB_DIR:/src/tatat_blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    extract-taxon-protein-sequences \
    7742 /src/tatat_blastdb/7742.fsa

# To build vertebrata blast database from sequences (~15 min on HPC)
singularity exec \
    --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    $SINGULARITY_IMAGE \
    bash /src/app/annotation/database_prep.sh \
    make_protein_blast_db \
    /src/blastdb/7742.fsa "vertebrata" 7742 /src/blastdb/vertebrata

# Remove original fsa
rm $TATAT_BLASTDB_DIR/7742.fsa

# To prepare vertebrata database for diamond usage (~2 min on HPC)
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    $SINGULARITY_IMAGE \
    diamond prepdb -d /src/blastdb/vertebrata