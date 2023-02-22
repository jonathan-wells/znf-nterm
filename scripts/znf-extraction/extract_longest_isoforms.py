#!/usr/bin/env python3

import sys
import re
from collections import defaultdict
from Bio import SeqIO

def get_gene_to_protein_dict(filename):
    """Uses GFF file to assign protein isoforms to their gene of origin."""
    rna_to_gene = {}
    protein_to_rna = {}
    with open(filename) as gff:
        for line in gff:
            if re.match(r'#', line):
                continue
            
            feature_type = re.search(r'\t(mRNA|CDS)\t', line)
            if feature_type:
                feature_type = feature_type.group(1)
            else:
                continue
            
            parent_data = re.search(r'\tID=.+?-(.+?);Parent=.+?-(.+?);', line)
            if parent_data:
                parent_data = parent_data.groups()
            else:
                raise ValueError(line)
            
            match feature_type:
                case "CDS":
                    protein_to_rna[parent_data[0]] = parent_data[1]
                case "mRNA":
                    rna_to_gene[parent_data[0]] = parent_data[1]

    gene_to_proteins = defaultdict(list)
    for protein, rna in protein_to_rna.items():
        gene = rna_to_gene.get(rna, None)
        gene_to_proteins[gene].append(protein)
    
    return gene_to_proteins

def get_protein_lengths(filename):
    """Returns dictionaries of protein lengths and original SeqIO records."""
    prot_lengths = {}
    records = {}
    for record in SeqIO.parse(filename, 'fasta'):
        prot_lengths[record.name] = len(record.seq)
        records[record.name] = record
    return prot_lengths, records

def extract_longest_isoform(seqfile, gfffile, outfile):
    """Extracts single longest isoform for set of protein isoforms associated with a gene"""
    prot_lengths, records = get_protein_lengths(seqfile)
    gene_to_proteins = get_gene_to_protein_dict(gfffile)
    
    longest_isoforms = [sorted(proteins, key=lambda x: prot_lengths.get(x, 0))[-1] for proteins in
                        gene_to_proteins.values()]  # Sorted by increasing length
    
    longest_isoform_records = [records[prot] for prot in longest_isoforms if records.get(prot, None)]

    SeqIO.write(longest_isoform_records, outfile, 'fasta')

def main(datadir, species):
    """Run if called by search_proteomes.sh script"""
    seqfile = f'{datadir}/seqs/{species}.aa.fa'
    gfffile = f'{datadir}/gffs/{species}.gff'
    outfile = f'{datadir}/seqs/{species}.longest_isoform.fa'
    extract_longest_isoform(seqfile, gfffile, outfile)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
