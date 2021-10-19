from Bio import SeqIO
import re
from collections import Counter

class PeptideMatcher:

    def __init__(self, fasta, peptides, flanks, grid, progress_dialog):
        self.fasta = fasta
        self.peptides = peptides
        self.flanks = flanks
        self.grid = grid
        self.progress_dialog = progress_dialog

    def parse_fasta(self):
        with open(self.fasta) as (fp):
            self.fasta_sequences = list(SeqIO.parse(fp, 'fasta'))
            assert len(self.fasta_sequences) > 0, 'The fasta file seems to be empty or mal-formed'

    def run(self):
        self.parse_fasta()
        peptide_num = 0
        with open(self.peptides) as (fp):
            row = 0
            for line in fp:
                peptide_num += 1
                peptide = line.rstrip().upper()
                peplen = len(peptide)
                matches = 0
                all_n_terms = []
                all_c_terms = []
                for i in range(self.flanks):
                    all_n_terms.append([])
                    all_c_terms.append([])

                for record in self.fasta_sequences:
                    index = record.seq.find(peptide)
                    seq_len = len(record.seq)
                    while index > -1:
                        start = index + 1
                        end = index + peplen
                        to_c_term = seq_len - end
                        if index > self.flanks:
                            n_term = record.seq[index - self.flanks:index]
                        else:
                            n_term = '[' + record.seq[0:index]
                        offset = self.flanks - len(n_term)
                        for i in reversed(range(self.flanks)):
                            let = str(n_term[(i - offset)])
                            all_n_terms[i].append(let)
                            if let == '[':
                                break

                        if to_c_term > self.flanks:
                            c_term = record.seq[end:end + self.flanks]
                        else:
                            c_term = record.seq[end:] + ']'
                        for i in range(self.flanks):
                            let = str(c_term[i])
                            all_c_terms[i].append(let)
                            if let == ']':
                                break

                        self.grid.AppendRows()
                        self.grid.SetCellValue(row, 0, peptide)
                        self.grid.SetCellValue(row, 1, str(peplen))
                        self.grid.SetCellValue(row, 2, record.id)
                        self.grid.SetCellValue(row, 3, str(start))
                        self.grid.SetCellValue(row, 4, str(end))
                        self.grid.SetCellValue(row, 5, str(to_c_term + 1))
                        self.grid.SetCellValue(row, 6, str(n_term))
                        self.grid.SetCellValue(row, 7, str(c_term))
                        index = record.seq.find(peptide, start)
                        matches += 1
                        row += 1
                        self.progress_dialog.Pulse()

                if matches == 0:
                    assert re.match('[ACDEFGHIKLMNPQRSTVWY]+$', peptide), "Mal-formed peptide string '%s'" % peptide
                    self.grid.AppendRows()
                    self.grid.SetCellValue(row, 0, peptide)
                    self.grid.SetCellValue(row, 1, str(peplen))
                    self.grid.SetCellValue(row, 2, 'No match')
                    row += 1
                else:
                    n_terms = ''
                    for pos in all_n_terms:
                        counts = Counter(pos)
                        n_terms += '{'
                        for aa in counts:
                            n_terms += '%d%s|' % (counts[aa], aa)

                        n_terms += '}'

                    c_terms = ''
                    for pos in all_c_terms:
                        counts = Counter(pos)
                        c_terms += '{'
                        for aa in counts:
                            c_terms += '%d%s|' % (counts[aa], aa)

                        c_terms += '}'

                    for row_i in range(row - matches, row):
                        self.grid.SetCellValue(row_i, 8, n_terms)
                        self.grid.SetCellValue(row_i, 9, c_terms)

        assert peptide_num > 0, 'The list of the peptides seems to be empty'
        assert matches > 0, 'No matches found'
