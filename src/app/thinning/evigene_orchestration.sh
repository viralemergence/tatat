#!/bin/bash

run_evigene() {
	cd $1 \
	&& ln -s $2 $3 \
    && head $3
}

"$@"