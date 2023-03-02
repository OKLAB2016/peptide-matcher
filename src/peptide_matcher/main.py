
def run_gui():

    from peptide_matcher import PMFrame
    import wx

    class PeptideMatcherApp(wx.App):
        def OnInit(self):
            self.frame = PMFrame(None, wx.ID_ANY, '')
            self.SetTopWindow(self.frame)
            self.frame.Show()
            return True

    peptide_matcher = PeptideMatcherApp(0)
    peptide_matcher.MainLoop()

def run_cli():

    from peptide_matcher import PeptideMatcher, wrap_logos, wrap_scores
    from argparse import ArgumentParser
    from sys import argv
    import json
    import csv
    import sys

    parser = ArgumentParser(description = 'Match peptides in a protein database.')
    parser.add_argument('--peptides', '-p', metavar = 'FILENAME', required = True, help = 'list of peptides to match')
    parser.add_argument('--database', '-d', metavar = 'FILENAME', required = True, help = 'protein database in fasta format')
    parser.add_argument('--secstruct', '-s', action = 'store_true', help = 'whether the database also contains structural information')
    parser.add_argument('--flanks', '-f', metavar = 'N', type = int, default = 4, help = 'length of the flanks to report (default: 4)')
    parser.add_argument('--format', '-F', default = 'json', choices = [ 'json', 'tsv', ], help = 'output format (default: json)')
    parser.add_argument('--output', '-o', default = '-', help = 'output file (default: output to stdout)')

    args = parser.parse_args()

    out = open(args.output, 'w') if args.output != '-' else sys.stdout

    if args.format == 'json':
        top = '['
        bottom = ']'
    elif args.format == 'tsv' or args.format == 'csv':
        delimiter = '\t' if args.format == 'tsv' else ','
        header = [ 'peptide', 'peplen', 'record_id', 'start', 'end', 'c_term', 'n_flank', 'c_flank', 'n_logos', 'c_logos', 'sst_n_term', 'sst_pept', 'sst_c_term', 'tm_n_term', 'tm_pept', 'tm_c_term', 'conf_n_term', 'conf_pept', 'conf_c_term', 'acc_n_term', 'acc_pept', 'acc_c_term' ]
        writer = csv.DictWriter(out, delimiter = delimiter, fieldnames = header)
        writer.writeheader()
        top = ''
        bottom = ''
    
    with open(args.peptides) as peptides:
        if top: out.write(top)
        pm = PeptideMatcher(peptides, args.database, args.secstruct, args.flanks)
        is_first = True
        for output in pm.run():
            if args.format == 'json':
                json_str = json.dumps(output)
                if not is_first: out.write(',')
                out.write(json_str)
            elif args.format == 'tsv' or args.format == 'csv':
                if output['matches']:
                    n_logos = wrap_logos(output['n_logos'])
                    c_logos = wrap_logos(output['c_logos'])
                    for match in output['matches']:
                        row = {
                            'peptide':   output['peptide'],
                            'peplen':    len(output['peptide']),
                            'record_id': match['record_id'],
                            'start':     match['start'],
                            'end':       match['end'],
                            'c_term':    match['c_term'],
                            'n_flank':   match['n_flank'],
                            'c_flank':   match['c_flank'],
                            'n_logos':   n_logos,
                            'c_logos':   c_logos
                        }
                        if args.secstruct:
                            row['sst_n_term']  = match['sst_n_term']
                            row['sst_pept']    = match['sst_pept']
                            row['sst_c_term']  = match['sst_c_term']
                            row['tm_n_term']   = match['tm_n_term']
                            row['tm_pept']     = match['tm_pept']
                            row['tm_c_term']   = match['tm_c_term']
                            row['conf_n_term'] = wrap_scores(match['conf_n_term'])
                            row['conf_pept']   = wrap_scores(match['conf_pept'])
                            row['conf_c_term'] = wrap_scores(match['conf_c_term'])
                            row['acc_n_term']  = wrap_scores(match['acc_n_term'])
                            row['acc_pept']    = wrap_scores(match['acc_pept'])
                            row['acc_c_term']  = wrap_scores(match['acc_c_term'])
                        writer.writerow(row)
                else:
                    row = {
                        'peptide': output['peptide'],
                        'peplen': len(output['peptide']),
                        'record_id': 'No match'
                    }
                    writer.writerow(row)
            is_first = False
        if bottom: print(bottom)
    if out is not sys.stdout:
        out.close()

