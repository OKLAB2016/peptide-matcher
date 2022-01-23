
from pandas import read_fwf
from collections import defaultdict 

dssp_files = snakemake.input['dssps']
pdb_files  = snakemake.input['pdbs']
output_file    = str(snakemake.output)

dssp_specs = {
    'resnumber': (0, 5),
    'resname': (5, 10),
    'chain': (10, 12),
    'aminoacid': (12, 14),
    'secstruct': (14, 17),
    'turns_helix_3': (17, 19),
    'turns_helix_4': (19, 20),
    'turns_helix_5': (20, 21),
    'geometrical_bend': (21, 22),
    'chirality': (22, 23),
    'beta_bridge_label_1': (23, 24),
    'beta_bridge_label_2': (24, 25),
    'beta_bridge_partner_resnum_1': (25, 29),
    'beta_bridge_partner_resnum_2': (29, 33),
    'beta_sheet_label': (33, 34),
    'accessibility': (34, 38),
    'N_H___O_00': (38, 45),
    'N_H___O_01': (45, 50),
    'O___N_H_00': (50, 56),
    'O___N_H_01': (57, 61),
    'N_H___O_10': (61, 67),
    'N_H___O_11': (68, 72),
    'O___N_H_10': (72, 78),
    'O___N_H_11': (79, 83),
    'tco': (83, 91),
    'kappa': (91, 97),
    'alpha': (97, 103),
    'phi': (103, 109),
    'psi': (109, 115),
    'x_ca': (115, 122),
    'y_ca': (122, 129),
    'z_ca': (129, 136)
}
maxasa_theoretical = {
    'A': 129,
    'C': 167,
    'D': 193,
    'E': 223,
    'F': 240,
    'G': 104,
    'H': 224,
    'I': 197,
    'K': 236,
    'L': 201,
    'M': 224,
    'N': 195,
    'P': 159,
    'Q': 225,
    'R': 274,
    'S': 155,
    'T': 172,
    'V': 174,
    'W': 285,
    'Y': 263
}
maxasa_empirical = {
    'A': 121,
    'C': 148,
    'D': 187,
    'E': 214,
    'F': 228,
    'G':  97,
    'H': 216,
    'I': 195,
    'K': 230,
    'L': 191,
    'M': 203,
    'N': 187,
    'P': 154,
    'Q': 214,
    'R': 265,
    'S': 143,
    'T': 163,
    'V': 165,
    'W': 264,
    'Y': 255
}
residues = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'XAA': 'X'
}

def parse_dssp(fname):
    colspecs = list(dssp_specs.values())
    names = list(dssp_specs.keys())
    with open(fname) as dssp:
        for line in dssp:
            if not line.rstrip().endswith('.'):
                break
        data = read_fwf(dssp, colspecs = colspecs, names = names)
    secstruct = []
    dssp_seq = ''
    acc = []
    for row in data.itertuples():
        if row.chain == 'A':
            aminoacid = str(row.aminoacid)
            struct = str(row.secstruct)
            asa = int(row.accessibility)
            rsa = round(asa / maxasa_theoretical[aminoacid] * 100)
            secstruct.append(struct if struct != 'nan' else '-')
            dssp_seq += aminoacid
            acc.append(str(rsa))
    return dssp_seq, secstruct, acc

def parse_pdb(pdb_file):
    start = stop = -1
    chain_seq = ''
    confs = []
    with open(pdb_file) as pdb_fh:
        for line in pdb_fh:
            if line.startswith('ATOM'):
                atom = line[13:17].strip()
                chain = line[20:22].strip()
                if atom == 'N' and chain == 'A':
                    res = line[17:20]
                    conf = float(line[60:66])
                    chain_seq += residues[res]
                    confs.append(conf)
            elif start < 0 and line.startswith('DBREF'):
                values = line.split()
                chain = values[2]
                if chain == 'A':
                    start = int(values[8]) - 1
                    stop  = int(values[9])
    return chain_seq, start, stop, confs
 

with open(output_file, 'w') as output:
    num_wins = len(dssp_files)
    accs  = {}
    ssts  = {}
    confs = {}
    seq   = {}
    count = defaultdict(int)
    for win in range(num_wins):
        dssp_seq, win_ssts, win_accs = parse_dssp(dssp_files[win])
        pdb_seq, start, stop, win_confs = parse_pdb(pdb_files[win])
        assert dssp_seq == pdb_seq, "Sequence mismatch between dssp and pdb"
        for i in range(stop - start):
            pos = i + start
            count[pos] += 1
            if pos in seq:
                assert seq[pos] == dssp_seq[i], "Sequence mismatch between windows"
            else:
                seq[pos] = dssp_seq[i]
            if pos not in confs or win_confs[i] > confs[pos]:
                confs[pos] = win_confs[i]
                accs[pos]  = win_accs[i]
                ssts[pos]  = win_ssts[i]
    output.write('pos,res,count,acc,sst,conf\n')
    for pos, num in sorted(count.items()):
        output.write('%d,%s,%d,%s,%s,%.2f\n' % (pos, seq[pos], num, accs[pos], ssts[pos], confs[pos]))
