#!/usr/bin/env python3

from Bio import SeqIO
import markov_clustering as mc
import networkx as nx
from collections import namedtuple
import random

def load_seqs(filename, samplesize=None):
    SimpleRecord = namedtuple('SimpleRecord', ['id', 'seq'])
    seqlist = []
    for record in SeqIO.parse(filename, 'fasta'):
        seqlist.append(SimpleRecord(record.id, str(record.seq)))
    if samplesize:
        return random.sample(seqlist, samplesize)
    else:
        return seqlist

def jaccard(seq1, seq2, kmersize=5):
    set1 = set(seq1[i:i+kmersize] for i in range(len(seq1) - kmersize + 1))
    set2 = set(seq2[i:i+kmersize] for i in range(len(seq2) - kmersize + 1))
    return len(set1.intersection(set2))/len(set1.union(set2))

def calculate_distances(seqlist, outfilename):
    pairs = 1
    with open(outfilename, 'w') as outfile:
        kmersize = 5
        for i, seq1 in enumerate(seqlist):
            for j, seq2 in enumerate(seqlist):
                if j >= i:
                    break
                elif len(seq1.seq) < kmersize:
                    continue
                elif len(seq2.seq) < kmersize:
                    continue
                jaccard_dist = jaccard(seq1.seq, seq2.seq, kmersize)
                if jaccard_dist:
                    outfile.write(f'{seq1.id}\t{seq2.id}\t{jaccard_dist:.3f}\n')
                if pairs % 1e5 == 0:
                    print(f'{pairs} checked')
                pairs += 1
    print(pairs)

def markov_clustering(weighted_edgelist, outfilename):
    graph = nx.read_weighted_edgelist(path=weighted_edgelist, delimiter='\t')
    matrix = nx.to_scipy_sparse_array(graph)
    result = mc.run_mcl(matrix)           # run MCL with default parameters
    clusters = mc.get_clusters(result)
    nodelist = list(graph.nodes)
    with open(outfilename, 'w') as outfile:
        for i, cluster in enumerate(clusters):
            for node in cluster:
                label = nodelist[node]
                outfile.write(f'{label}\t{i}\n')
        
def main():
    seqlist = load_seqs('../../data/seqs/unknown_nterm.fa', 30000)
    calculate_distances(seqlist, '../../data/unknown_nterm_adjacency_list.txt')
    markov_clustering('../../data/unknown_nterm_adjacency_list.txt',
                      '../../data/unknown_nterm_clusters.txt')
    
if __name__ == "__main__":
    main()
