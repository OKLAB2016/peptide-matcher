from Bio import SeqIO
from random import randint
import re

fasta_file = str(snakemake.input)
peptide_file = str(snakemake.output)

min_len = snakemake.params['min_len']
max_len = snakemake.params['max_len']


with open(fasta_file) as fasta_fh:
    with open(peptide_file, 'w') as peptide_fh:
        fasta = SeqIO.parse(fasta_fh, 'fasta')
        for record in fasta:
            seq_len = randint(min_len, max_len)
            if seq_len < len(record.seq):
                start = randint(0, len(record.seq) - seq_len)
                peptide = str(record.seq[start:start + seq_len])
                if re.match('[ARNDCQEGHILKMFPSTWYV]+$', peptide):
                    peptide_fh.write(peptide + '\n')
