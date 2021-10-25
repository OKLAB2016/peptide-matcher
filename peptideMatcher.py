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
        self.conf_re = re.compile('confidence:([^\s]+)')

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
                    secstruct_match = self.secstruct_re.search(record.description)
                    acc_match = self.acc_re.search(record.description)
                    conf_match = self.conf_re.search(record.description)
                    # assert secstruct, "Record description does not contain secondary structure: " + record.description
                    # assert acc, "Record description does not contain accessibility: " + record.description
                    # assert conf, "Record description does not contain confidence scores: " + record.description
                else:
                    secstruct_match = ''
                    acc_match = []
                    conf_match = []
                secstruct = secstruct_match.group(1) if secstruct_match else ''
                acc = acc_match.group(1).split(',') if acc_match else []
                conf = conf_match.group(1).split(',') if conf_match else []
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
                    if secstruct:
                        ss_pept = secstruct[start_index:end]
                        if start_index > self.flanks:
                            ss_n_term = secstruct[start_index - self.flanks:start_index]
                        else:
                            ss_n_term = '[' + secstruct[0:start_index]
                        if to_c_term > self.flanks:
                            ss_c_term = secstruct[end:end + self.flanks]
                        else:
                            ss_c_term = secstruct[end:] + ']'
                    else:
                        ss_pept = ss_n_term = ss_c_term = ''

                    if acc:
                        acc_pept = ','.join(acc[start_index:end])
                        if start_index > self.flanks:
                            acc_n_term = ','.join(acc[start_index - self.flanks:start_index])
                        else:
                            acc_n_term = '[' + ','.join(acc[0:start_index])
                        if to_c_term > self.flanks:
                            acc_c_term = ','.join(acc[end:end + self.flanks])
                        else:
                            acc_c_term = ','.join(acc[end:]) + ']'
                    else:
                        acc_pept = acc_n_term = acc_c_term = ''
                    if conf:
                        conf_pept = ','.join(conf[start_index:end])
                        if start_index > self.flanks:
                            conf_n_term = ','.join(conf[start_index - self.flanks:start_index])
                        else:
                            conf_n_term = '[' + ','.join(conf[0:start_index])
                        if to_c_term > self.flanks:
                            conf_c_term = ','.join(conf[end:end + self.flanks])
                        else:
                            conf_c_term = ','.join(conf[end:]) + ']'
                    else:
                        conf_pept = conf_n_term = conf_c_term = ''
                    data[peptide].append((record.id, start, end, n_term, c_term, to_c_term, ss_n_term, ss_pept, ss_c_term, conf_n_term, conf_pept, conf_c_term, acc_n_term, acc_pept, acc_c_term))
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
                for record_id, start, end, n_term, c_term, to_c_term, *struct_info in data_peptide:

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

                    if self.secstruct_included:
                        ss_n_term, ss_pept, ss_c_term, conf_n_term, conf_pept, conf_c_term, acc_n_term, acc_pept, acc_c_term = struct_info
                        self.grid.SetCellValue(row, 10, ss_n_term)
                        self.grid.SetCellValue(row, 11, ss_pept)
                        self.grid.SetCellValue(row, 12, ss_c_term)
                        self.grid.SetCellValue(row, 13, conf_n_term)
                        self.grid.SetCellValue(row, 14, conf_pept)
                        self.grid.SetCellValue(row, 15, conf_c_term)
                        self.grid.SetCellValue(row, 16, acc_n_term)
                        self.grid.SetCellValue(row, 17, acc_pept)
                        self.grid.SetCellValue(row, 18, acc_c_term)
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
