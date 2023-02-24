#!/usr/bin/env python3

import sys
# from numpy import random
import random
from Bio import SeqIO
from collections import Counter

def extract_kmers(seq, k=3, window=None):
    if window == None:
        seq = list(seq)
    else:
        seq = [c for c in seq[window[0]:window[-1]]]
    seq_kmers = [''.join(seq[i:i+k]) for i in range(len(seq)-k)]
    return seq_kmers

def main():
    seq_kmers = []
    background_kmers = []
    seqs = SeqIO.parse(sys.argv[1], 'fasta')
    for seq in seqs:
        skmer = extract_kmers(str(seq.seq), k=int(sys.argv[2]))
        seq_kmers += skmer
    seq_counter = Counter(seq_kmers)
    for (k, v) in sorted(seq_counter.items(), key=lambda x: x[1]):
        print(k, v)

if __name__ == "__main__":
    main()
