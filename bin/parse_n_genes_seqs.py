import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import concurrent.futures

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

def parse_seq_info(input_file):
    with open(input_file, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            seq_id = row[2].split(" ")[0].split("_")[1]
            seq = row[3]
            sec_rec = SeqRecord(Seq(seq), id=seq_id, description="")
            yield sec_rec

def process_file(input_file):
    base = os.path.basename(input_file).split(".")[0]
    output_file = f"/home/glbrc.org/millican/ref_db/nitrogen-cycle/{base}.faa.gz"
    seq_records = parse_seq_info(input_file)
    with gzip.open(output_file, 'wt') as output_file:
        SeqIO.write(seq_records, output_file, 'fasta')


def main():
    file_list = glob.glob("/home/glbrc.org/millican/repos/trait-mapper/nitrogen-cycle/*.txt")
    with concurrent.futures.ThreadPoolExecutor(max_workers=24) as executor:
        executor.map(process_file, file_list)