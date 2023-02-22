#!/usr/bin/env bash

DATADIR="../../data"
ZNF="PF00096.hmm"

c=1
while read -r line; do
    # if [[ $c -gt 3 ]]; then
    #     break
    # fi

    accession=$(awk -F '\t' '{ OFS="\t"; print $1 }' <<< $line)
    species=$(awk -F '\t' '{ OFS="\t"; print $3 }' <<< $line | sed 's/ /_/g')
   
    if [[ $accession == "Assembly Accession" ]]; then
        continue
    fi

    if [[ -s "${DATADIR}/seqs/${species}_znfs.fa" ]]; then
        echo "${species} data already exists"
        continue
    fi

    echo "Downloading ${species}..."
    datasets download genome accession \
        $accession \
        --include "protein,gff3" \
        --no-progressbar \
        --filename "${accession}.zip"
    
    if ! [ -s "${accession}.zip" ]; then
        echo "${species} genome not downloaded"
        continue
    fi

    echo "Unzipping genome..."
    unzip -qq "${accession}.zip"
    mv "./ncbi_dataset/data/${accession}/protein.faa" "${DATADIR}/seqs/${species}.aa.fa" 
    mv "./ncbi_dataset/data/${accession}/genomic.gff" "${DATADIR}/gffs/${species}.gff"
    
    echo "Extracting longest isoforms..."
    ./extract_longest_isoforms.py $DATADIR $species
    rm "${DATADIR}/seqs/${species}.aa.fa"  # Redundant fasta containing irrelevant isoforms 

    echo "Searching for ZNF sequences..."
    hmmsearch \
        -o tmp.out \
        --tblout "${DATADIR}/hmmer-out/${species}_znf.out" \
        --domtblout "${DATADIR}/hmmer-out/${species}_znf_domains.out" \
        --noali \
        -E 0.01 \
        --domE 0.01 \
        --incE 0.01 \
        --incdomE 0.01 \
        "${DATADIR}/phmms/${ZNF}" \
        "${DATADIR}/seqs/${species}.longest_isoform.fa"

    echo "Extracting ZNFs from HMMER output"
    ./parse_hmmer.py $DATADIR $species

    echo "Cleaning up..."
    rm -r ./README.md ./ncbi_dataset "./${accession}.zip"
    rm tmp.out
    rm "${DATADIR}/seqs/${species}.longest_isoform.fa"
    
    ((c=c+1))

done < "${DATADIR}/refseq_metazoans.tsv"
