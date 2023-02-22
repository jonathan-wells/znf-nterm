#!/usr/bin/env python3

import re
import sys
from Bio import SeqIO

def filter_znfs(filename, min_znf_domains=5, max_znf_domains=40):
    """Return list of ZNF proteins with between min and and max domain copies."""
    filtered_znf_accessions = []
    with open(filename) as infile:
        for line in infile:
            if re.match('#', line):
                continue
            line = re.findall(r'(.+?)\s+', line)
            line = line[:4] + \
                   [float(x) for x in line[4:11]] + \
                   [int(x) for x in line[11:18]] + \
                   [' '.join(line[18:])]
            if line[16] >= min_znf_domains and line[16] <= max_znf_domains:
                filtered_znf_accessions.append(line[0])
    return set(filtered_znf_accessions)

def main(datadir, species):
    hmmfile = f'{datadir}/hmmer-out/{species}_znf.out'
    seqfile = f'{datadir}/seqs/{species}.longest_isoform.fa'
    outfile = f'{datadir}/seqs/{species}_znfs.fa'
    
    znf_accessions = filter_znfs(hmmfile)
    filtered_records = []
    for record in SeqIO.parse(seqfile, 'fasta'):
        if record.name in znf_accessions:
            filtered_records.append(record)
    SeqIO.write(filtered_records, outfile, 'fasta')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
