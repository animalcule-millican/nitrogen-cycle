#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip


def parse_tax_info(input):
    tax = None
    tax_id = None
    parts = input.split(" ")
    for part in parts:
        if part.startswith('Tax='):
            tax = part[4:]
        elif part.startswith('TaxID='):
            tax_id = part[6:]
        if tax and tax_id is not None:
            return tax, tax_id

def parse_seq_info(input_file, tax, tax_id):
    with open(input_file, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            seq_id = row[2].split(" ")[0].split("_")[1]
            seq = row[3]
            sec_rec = SeqRecord(Seq(seq), id=seq_id, description=f"{tax}|{tax_id}")
            yield sec_rec

def process_file(input_file, output_file):
    seq_records = parse_seq_info(input_file)
    with gzip.open(output_file, 'wt') as output_file:
        SeqIO.write(seq_records, output_file, 'fasta')


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file)
    