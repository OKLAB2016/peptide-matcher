from Bio import SeqIO
import re
from time import time
from collections import Counter, defaultdict
from ahocorasick import Automaton

class PeptideMatcher:

    def __init__(self, fasta, secstruct_included, peptides, flanks, grid, progress_dialog):
        self.fasta = fasta
        self.secstruct_included = secstruct_included
        self.peptides = peptides
        self.flanks = flanks
        self.flanks_range = list(range(flanks))
        self.flanks_revrange = list(reversed(range(flanks)))
        self.grid = grid
        self.progress_dialog = progress_dialog
        self.secstruct_re = re.compile('secstruct:([^\s]+)')
        self.acc_re = re.compile('accessibility:([^\s]+)')

    def parse_peptides(self):
        self.automaton = Automaton()
        self.peptide_seqs = []
        with open(self.peptides) as fp:
            for line in fp:
               peptide = line.rstrip().upper()
               assert re.match('[ACDEFGHIKLMNPQRSTVWY]+$', peptide), "Malformed peptide string '%s'" % peptide
               self.automaton.add_word(peptide, peptide)
               self.peptide_seqs.append(peptide)
        assert len(self.automaton) > 0, "The peptide file seems to be empty"
        self.automaton.make_automaton()

    def run(self):
        time_start = time()
        self.parse_peptides()
        data = defaultdict(list)
        with open(self.fasta) as fp:
            for record in SeqIO.parse(fp, 'fasta'):
                record_seq = str(record.seq)
                seq_len = len(record_seq)
                if self.secstruct_included:
                    secstruct = self.secstruct_re.search(record.description)
                    acc = self.acc_re.search(record.description)
                    assert secstruct and acc, "Record description does not contain secstruct or accessibility information: " + record.description
                    secstruct = secstruct.group(1)
                    acc = acc.group(1).split(',')
                for end_index, peptide in self.automaton.iter(record_seq):
                    start_index = end_index - len(peptide) + 1
                    start = start_index + 1
                    end = end_index + 1
                    to_c_term = seq_len - end
                    if start_index > self.flanks:
                        n_term = record_seq[start_index - self.flanks:start_index]
                    else:
                        n_term = '[' + record_seq[0:start_index]
                    if to_c_term > self.flanks:
                        c_term = record_seq[end:end + self.flanks]
                    else:
                        c_term = record_seq[end:] + ']'
                    if self.secstruct_included:
                        peptide_secstruct = '%s^%s|%s^%s' % (
                            secstruct[start_index - 1] if start_index > 0 else '[',
                            secstruct[start_index],
                            secstruct[end_index],
                            secstruct[end_index + 1] if to_c_term > 0 else ']'
                        )
                        peptide_acc = '%s^%s|%s^%s' % (
                            acc[start_index - 1] if start_index > 0 else '[',
                            acc[start_index],
                            acc[end_index],
                            acc[end_index + 1] if to_c_term > 0 else ']'
                        )
                    else:
                        peptide_secstruct = ''
                        peptide_acc = ''
                    data[peptide].append((record.id, start, end, n_term, c_term, to_c_term, peptide_secstruct, peptide_acc))
        row = 0
        for peptide in self.peptide_seqs:
            peplen = len(peptide)
            if peptide in data:
                matches = 0
                all_n_terms = []
                all_c_terms = []
                for i in range(self.flanks):
                    all_n_terms.append([])
                    all_c_terms.append([])
                data_peptide = data[peptide]
                for record_id, start, end, n_term, c_term, to_c_term, peptide_secstruct, peptide_acc in data_peptide:

                    offset = self.flanks - len(n_term)

                    for i in self.flanks_revrange:
                        let = n_term[i - offset]
                        all_n_terms[i].append(let)
                        if let == '[': break
                    for i in self.flanks_range:
                        let = c_term[i]
                        all_c_terms[i].append(let)
                        if let == ']': break

                    self.grid.AppendRows()
                    self.grid.SetCellValue(row, 0, peptide)
                    self.grid.SetCellValue(row, 1, str(peplen))
                    self.grid.SetCellValue(row, 2, record_id)
                    self.grid.SetCellValue(row, 3, str(start))
                    self.grid.SetCellValue(row, 4, str(end))
                    self.grid.SetCellValue(row, 5, str(to_c_term + 1))
                    self.grid.SetCellValue(row, 6, str(n_term))
                    self.grid.SetCellValue(row, 7, str(c_term))
                    self.grid.SetCellValue(row, 10, peptide_secstruct)
                    self.grid.SetCellValue(row, 11, peptide_acc)
                    row += 1

                n_terms = ''
                for pos in all_n_terms:
                    counts = Counter(pos)
                    logo = []
                    for aa in counts:
                        logo.append(str(counts[aa]) + aa)
                    n_terms += '{' + '|'.join(logo) + '}'

                c_terms = ''
                for pos in all_c_terms:
                    counts = Counter(pos)
                    logo = []
                    for aa in counts:
                        logo.append(str(counts[aa]) + aa)
                    c_terms += '{' + '|'.join(logo) + '}'

                matches = len(data_peptide)
                for i in range(matches):
                    self.grid.SetCellValue(row - i - 1, 8, n_terms)
                    self.grid.SetCellValue(row - i - 1, 9, c_terms)
                self.progress_dialog.Pulse()
            else:
                self.grid.AppendRows()
                self.grid.SetCellValue(row, 0, peptide)
                self.grid.SetCellValue(row, 1, str(peplen))
                self.grid.SetCellValue(row, 2, 'No match')
                row += 1
        print("Run time: " + str(time() - time_start))
