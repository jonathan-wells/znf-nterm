#!/bin/bash

cat ../../data/seqs/*znfs_nterm.fa > ../../data/seqs/all_nterm.fa
cat ../../data/seqs/*znfs_cterm.fa > ../../data/seqs/all_cterm.fa

echo "Searching for known domains in N' sequence"
hmmsearch \
    -o tmp.out \
    --tblout ../../data/hmmer-out/all_nterm.out \
    --domtblout ../../data/hmmer-out/all_nterm_domains.out \
    --noali \
    -E 0.01 \
    --domE 0.01 \
    --incE 0.01 \
    --incdomE 0.01 \
    --cpu 14 \
    ~/Genomes/Pfam/Pfam-A.hmm \
    ../../data/seqs/all_nterm.fa

rg -v '#' ../../data/hmmer-out/all_nterm.out | awk '{ print $1 }' > known_nterm_domains.names

seqkit grep \
    -v \
    -f known_nterm_domains.names \
    ../../data/seqs/all_nterm.fa > ../../data/seqs/unknown_nterm.fa

rm known_nterm_domains.names

echo "Searching for known domains in C' sequence"
hmmsearch \
    -o tmp.out \
    --tblout ../../data/hmmer-out/all_cterm.out \
    --domtblout ../../data/hmmer-out/all_cterm_domains.out \
    --noali \
    -E 0.01 \
    --domE 0.01 \
    --incE 0.01 \
    --incdomE 0.01 \
    --cpu 14 \
    ~/Genomes/Pfam/Pfam-A.hmm \
    ../../data/seqs/all_cterm.fa

rg -v '#' ../../data/hmmer-out/all_cterm.out | awk '{ print $1 }' > known_cterm_domains.names

seqkit grep \
    -v \
    -f known_cterm_domains.names \
    ../../data/seqs/all_cterm.fa > ../../data/seqs/unknown_cterm.fa

rm known_cterm_domains.names
