#!/bin/bash

cd-hit \
    -i ../../data/seqs/unknown_nterm.fa \
    -o ../../data/seqs/unknown_nterm_clustered_70.fa \
    -c 0.7 \
    -sf 1 \
    -sc 1 \
    -d 0 \
    -g 1  # Sensitive clustering

mv ../../data/seqs/unknown_nterm_clustered_70.fa.clstr ../../data
