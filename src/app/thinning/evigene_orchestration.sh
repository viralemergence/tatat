#!/bin/bash

run_evigene() {
	cd $1 \
	&& ln -s $2 $3 || true \
	&& echo "Starting evigene" \
    && $EVIGENE/scripts/prot/tr2aacds4.pl -NCPU $4 -MAXMEM $5 -log -cdna $3
}

"$@"