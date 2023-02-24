#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import re
from collections import defaultdict


def load_znf_doms(species):
    """Takes domtblout file and extracts start and end coordinates of first and last znf domains."""
    hmmfile = f'../../data/hmmer-out/{species}_znf_domains.out'
    znfdict = defaultdict(list)
    
    with open(hmmfile) as infile:
        for line in infile:
            if re.match('#', line):
                continue
            line = re.findall(r'(.+?)\s+', line)
            target = line[0]
            start, end = int(line[19]), int(line[20])
            znfdict[target].append(start)
            znfdict[target].append(end)
    
    for key, val in znfdict.items():
        sval = sorted(val)
        znfdict[key] = [sval[0], sval[-1]]
    
    return znfdict

def extract_terminal_sequence(species, znf_doms):
    seqfile = f'../../data/seqs/{species}_znfs.fa'
    nterm_records, cterm_records = [], []
    for record in SeqIO.parse(seqfile, 'fasta'):
        coords = znf_doms.get(record.name, None)
        if coords == None:
            continue
        
        nterm_id = record.name + '_nterm'
        nterm_seq = record.seq[:coords[0]]
        nterm_records.append(SeqRecord(id=nterm_id, seq=nterm_seq, description=record.description))

        cterm_id = record.name + '_cterm'
        cterm_seq = record.seq[coords[1]+1:]
        cterm_records.append(SeqRecord(id=cterm_id, seq=cterm_seq, description=record.description))
    
    with (open(f'../../data/seqs/{species}_znfs_nterm.fa', 'w') as nfile,
          open(f'../../data/seqs/{species}_znfs_cterm.fa', 'w') as cfile):
        SeqIO.write(nterm_records, nfile, 'fasta')
        SeqIO.write(cterm_records, cfile, 'fasta')

def main():
    with open('../../data/refseq_metazoans.tsv') as infile:
        for line in infile:
            line = line.split('\t')
            species = line[2].replace(' ', '_')
            try:
                znf_doms = load_znf_doms(species)
                extract_terminal_sequence(species, znf_doms)
            except:
                print(f'No data for {species}')

if __name__ == "__main__":
   main() 
