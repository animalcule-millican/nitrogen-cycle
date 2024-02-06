#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

rec_list = []
input_file = sys.argv[1]
output_file = sys.argv[2]
with open(input_file, 'r') as f:
    with open(output_file, 'w') as o:
        for line in f:
            row = line.strip().split("\t")
            seq_record = SeqRecord(Seq(row[1]), id=row[0], description="")
            rec_list.append(seq_record)
        SeqIO.write(rec_list, o, "fasta")