# peptide-matcher

peptide-matcher is a piece of software that can be used for matching peptide sequences identified in proteomics experiments using a database-match or a de novo approach against a sequence database. The main purpose is to extract sequence context for the corresponding matches, but peptide-matcher can also provide structural context if provided with a database that includes structural information, see [peptide-matcher-data](https://github.com/OKLAB2016/peptide-matcher-data/releases/).

There are three ways of how to use peptide-matcher:

* the GUI `peptide\_matcher\_gui`
* the CLI `peptide\_matcher`
* the python class `PeptideMatcher` -- `from peptide_matcher import PeptideMatcher`

The gui is written with [wxWidgets](https://www.wxwidgets.org/). Other dependencies include [biopython](https://biopython.org/) and [pyahocorasick](https://pyahocorasick.readthedocs.io/).

## Installation

Install via pipy `pip install peptide_matcher`.

## How to use the GUI

![interface](doc/interface.png)

Two files are needed: the database in fasta format with optional structural annotations (see below) and a plain list of peptide sequences.

The optional structural annotations should follow a custom format. Databases generated based on alphafold's models for a couple of popular model organisms are distributed at [peptide-matcher-data](https://github.com/OKLAB2016/peptide-matcher-data/releases/).

The results of the peptide matching are returned to the GUI and can be saved as xlsx. For each peptide the following output is generated:

| Field          | Description                                    | Example                       | Values                                            |
|----------------|------------------------------------------------|-------------------------------|---------------------------------------------------|
| `Peptide`      | peptide sequence                               | QVHAVSFYSK                    | string of amino acid symbols                      |
| `Length`       | peptide length                                 | 10                            | integer                                           |
| `Protein`      | matching protein id                            | A6NL46                        | string                                            |
| `Start`        | start position (1-based)                       | 150                           | integer                                           |
| `End`          | end position (1-based)                         | 159                           | integer                                           |
| `C-term`       | distance to protein's C-terminus               | 182                           | integer                                           |
| `N-flank`      | N-flanking residues in this protein            | TDKA                          | string                                            |
| `C-flank`      | C-flanking residues in this protein            | GHGV                          | string                                            |
| `N-flank*`     | weblogo for each position of the N-flank       | 2T\|2D\|2K\|2A                | `\|` - separator between positions                |
| `C-flank*`     | weblogo for each position of the C-flank       | 1G1D\|2H\|1G1E\|2V            |                                                   |
| `N-flank SS`   | secondary structure for the N-flank            | HHHH                          | string of DSSP codes                              |
| `Peptide SS`   | same for the peptide itself                    | HH------EE                    |                                                   |
| `C-flank SS`   | same for the C-flank region                    | EEEE                          |                                                   |
| `N-flank TM`   | transmembrane region for the N-flank           | TTTT                          | string of: `T` - TM region, `S` - signal peptide  |
| `Peptide TM`   | same for the peptide itself                    | TT--------                    |                                                   |
| `C-flank TM`   | same for the C-flank region                    | ----                          |                                                   |
| `N-flank conf` | alphafold's pLDDT score for the N-flank        | 43,46,40,49                   | list of integers 0-100                            |
| `Peptide conf` | same for the peptide itself                    | 44,44,45,44,50,39,48,39,56,46 |                                                   |
| `C-flank conf` | same for the C-flank                           | 49,47,42,46                   |                                                   |
| `N-flank RSA`  | relative solvent accessibility for the N-flank | 81,79,84,71                   | list of integers 0-100                            |
| `Peptide RSA`  | same for the peptide itself                    | 90,78,75,78,54,62,73,84,73,81 |                                                   |
| `C-flank RSA`  | same for the C-flank                           | 67,78,71,80                   |                                                   |

In the provided database, the RSA values are calculated by dividing the absolute solvent accessibility (ASA) as produced by dssp (mkdssp v.3.0.0) by the theoretical maximum values for ASA from [Tien et al 2013](https://dx.doi.org/10.1371%2Fjournal.pone.0080635).

## How to use CLI

Check out `peptide_matcher -h`:

```
$ peptide_matcher -h
usage: peptide_matcher [-h] --peptides FILENAME --database FILENAME [--secstruct] [--flanks N] [--format {json,tsv,csv}] [--output OUTPUT]

Match peptides in a protein database.

optional arguments:
  -h, --help            show this help message and exit
  --peptides FILENAME, -p FILENAME
                        list of peptides to match
  --database FILENAME, -d FILENAME
                        protein database in fasta format
  --secstruct, -s       whether the database also contains structural information
  --flanks N, -f N      length of the flanks to report (default: 4)
  --format {json,tsv,csv}, -F {json,tsv,csv}
                        output format (default: json)
  --output OUTPUT, -o OUTPUT
                        output file (default: output to stdout)
```

The output is similar to that of the GUI. The header of the tabular output formats looks as follows: `[ 'peptide', 'peplen', 'record_id', 'start', 'end', 'c_term', 'n_flank', 'c_flank', 'n_logos', 'c_logos', 'sst_n_term', 'sst_pept', 'sst_c_term', 'tm_n_term', 'tm_pept', 'tm_c_term', 'conf_n_term', 'conf_pept', 'conf_c_term', 'acc_n_term', 'acc_pept', 'acc_c_term' ]`. The json output is a list of dictionaries with each one of the following format: `{"peptide": "IYGALAVGAP", "matches": [{"record_id": "P77549", "start": 157, "end": 166, "c_term": 227, "n_flank": "NGMA", "c_flank": "LGLL", "sst_n_term": "HHHH", "sst_pept": "HHHHHHHHHH", "sst_c_term": "HHHH", "tm_n_term": "----", "tm_pept": "----------", "tm_c_term": "----", "conf_n_term": [94, 89, 91, 94], "conf_pept": [93, 86, 88, 94, 89, 85, 90, 92, 86, 88], "conf_c_term": [93, 94, 91, 94], "acc_n_term": [3, 6, 25, 6], "acc_pept": [9, 24, 19, 0, 25, 50, 44, 0, 22, 45], "acc_c_term": [36, 0, 37, 47]}], "n_logos": [{"N": 1}, {"G": 1}, {"M": 1}, {"A": 1}], "c_logos": [{"L": 1}, {"G": 1}, {"L": 1}, {"L": 1}]}`.

## How to use the API

```
from peptide_matcher import PeptideMatcher, wrap_logos, wrap_scores
peptides = [ 'IYGALAVGAP', 'LTCDETPVFSGSVLN', 'KRFARESGMTLL', 'GAGFAELLSSLQTPEIK', 'RTGHKLV' ] # or a file handle
database = 'UP000000625_83333_ECOLI.fasta' # or a file handle
flanks = 4
secstruct = True
pm = PeptideMatcher(peptides, database, secstruct, flanks)
for output in pm.run():
    print(output)
```
