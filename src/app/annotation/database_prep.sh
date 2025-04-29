#!/bin/bash

pull-nr-database() {
	cd $BLASTDB \
	&& perl "/src/tools/ncbi-blast-2.16.0+/bin/update_blastdb.pl" --passive --decompress nr
}

extract-taxon-protein-sequences() {
	blastdbcmd -db nr -taxids $1 -out $2
}

make_protein_blast_db() {
	makeblastdb -in $1 -dbtype prot -title $2 -parse_seqids -blastdb_version 5 -taxid $3 -out $4
}

# Core NT blast database to vertebrata blast database commands
pull-core-nt-database() {
	cd $BLASTDB \
	&& perl "/src/tools/ncbi-blast-2.16.0+/bin/update_blastdb.pl" --passive --decompress core_nt
}

extract-taxon-nucleotide-sequences() {
	blastdbcmd -db core_nt -dbtype nucl -taxids $1 -out $2
}

make_nucleotide_blast_db() {
	makeblastdb -in $1 -dbtype nucl -title $2 -parse_seqids -blastdb_version 5 -taxid $3 -out $4
}

"$@"