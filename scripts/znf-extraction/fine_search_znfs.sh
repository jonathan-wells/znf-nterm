#!/usr/bin/env bash

DATADIR="../../data"
ZNF="PF00096.hmm"

while read -r line; do

    species=$(awk -F '\t' '{ OFS="\t"; print $3 }' <<< $line | sed 's/ /_/g')
   
    echo "Searching for ZNF sequences in ${species}..."
    hmmsearch \
        -o tmp.out \
        --domtblout "${DATADIR}/hmmer-out/${species}_znf_domains.out" \
        --noali \
        -E 0.1 \
        --cpu 8 \
        "${DATADIR}/phmms/${ZNF}" \
        "${DATADIR}/seqs/${species}_znfs.fa"

done < "${DATADIR}/refseq_metazoans.tsv"

./extract_terminal_domains.py
