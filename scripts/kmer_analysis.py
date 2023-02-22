#!/usr/bin/env python3

import sys
# from numpy import random
import random
from Bio import SeqIO
from collections import Counter

def extract_kmers(seq, k=2, window=(0, 80)):
    seq = [c for c in seq[window[0]:window[-1]]]
    seq_kmers = [''.join(seq[i:i+k]) for i in range(len(seq)-k)]
    background = random.choices(seq, k=5000)
    background_kmers = [''.join(background[i:i+k]) for i in range(len(background)-k)]
    return seq_kmers, background_kmers 

def main():
    seq_kmers = []
    background_kmers = []
    seqs = SeqIO.parse(sys.argv[1], 'fasta')
    for seq in seqs:
        skmer, bgkmer = extract_kmers(str(seq.seq))
        seq_kmers += skmer
        background_kmers += bgkmer
    seq_counter = Counter(seq_kmers)
    bg_counter = Counter(background_kmers)
    for (k, v) in sorted(seq_counter.items(), key=lambda x: x[1]):
        print(k, v, bg_counter[k]/5000, v/(bg_counter.get(k, 1)/5000))

if __name__ == "__main__":
    main()
