#!/bin/bash

pull-nr-database() {
	cd $BLASTDB \
	&& perl "/src/tools/ncbi-blast-2.16.0+/bin/update_blastdb.pl" --passive --decompress nr
}

extract-taxon-protein-sequences() {
	blastdbcmd -db nr -taxids $1 -out $2
}

"$@"