
import re
from pandas import read_csv
from os.path import basename
from Bio import SeqIO
from Bio.Seq import Seq
from sys import stderr
from itertools import groupby

parsed_files  = snakemake.input['parsed']
uniprot_file  = str(snakemake.input['uniprot'])
fasta_file    = str(snakemake.output)

af_regex = re.compile(r'AF-(\w+)')
def parse_name(fname):
    search = af_regex.match(basename(fname))
    assert search, 'Unexpected file name: %s' % fname
    return search.group(1)

fnames = {}
for fname in parsed_files:
    record_id = parse_name(fname)
    fnames[record_id] = fname

def compress(s):
    return ''.join(str(sum(1 for x in g)) + k for k, g in groupby(s))
def to_hex(decs):
    return ''.join('%02x' % x for x in decs)

feature_types = {
    'transmembrane region': 'T',
    'signal peptide':       'S'
}
with open(fasta_file, 'w') as fasta:
    uniprot = SeqIO.parse(uniprot_file, 'uniprot-xml')
    for record in uniprot:
        TMs = {}
        for feature in record.features:
            if feature.type in feature_types and type(feature.location.end).__name__ == 'ExactPosition':
                for i in range(feature.location.start - 1, feature.location.end):
                    TMs[i] = feature_types[feature.type]
        if record.id in fnames:
            data = read_csv(fnames[record.id], index_col = 0)
            acc = []
            sst = []
            conf = []
            tm = []
            match = True
            for i in range(len(record.seq)):
                if i in data['res']:
                    match = data['res'][i] == record.seq[i]
                    if not match:
                        print("Sequence mismatch in %s" % record.id, file = stderr)
                        break
                    acc.append(data['acc'][i])
                    sst.append(data['sst'][i])
                    conf.append(round(data['conf'][i]))
                else:
                    acc.append(-1)
                    sst.append('_')
                    conf.append(-1)
                tm.append(TMs[i] if i in TMs else '-')
            record.description += ' confidence:%s secstruct:%s accessibility:%s' % (to_hex(conf), compress(sst), to_hex(acc))
            if match:
                record.description += ' transmembrane:%s' % compress(tm)
            else:
                record.seq = Seq(''.join(data['res']))
        SeqIO.write(record, fasta, 'fasta')
