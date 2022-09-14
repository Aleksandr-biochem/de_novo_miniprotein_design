#! /usr/bin/env python
import glob
import argparse

### This script extracts protein sequences from pdb files

parser = argparse.ArgumentParser(description='Extract protein sequences from pdb files')
parser.add_argument('-i', nargs=1, help='A folder with pdb files to extract', type=str, required=False)
parser.add_argument('-f', nargs=1, help='A single pdb file to extract', type=str, required=False)
args = parser.parse_args()

if not args.i is None:
  pdb_files = glob.glob(f"{args.i[0]}/*pdb")
elif not args.f is None:
  pdb_files = [args.f[0]]
else:
  raise Exception('No input provided. Provide either folder or file to parse')

# a dictionary to convert canonical residue names
three_to_one_leter_name = {"GLY": "G",
                           "ALA": "A",
                           "VAL": "V",
                           "LEU": "L",
                           "ILE": "I",
                           "SER": "S",
                           "THR": "T",
                           "MET": "M",
                           "CYS": "C",
                           "ASN": "N",
                           "GLN": "Q",
                           "ASP": "D",
                           "GLU": "E",
                           "LYS": "K",
                           "ARG": "R",
                           "HIS": "H",
                           "PHE": "F",
                           "TRP": "W",
                           "TYR": "Y",
                           "PRO": "P"}

### extract AA sequences
for file in pdb_files:
    name = file.split('/')[-1][:-4]

    seq, res_ids = [], []
    with open(file, 'r') as f:
      for line in f:
          data = line.strip().split()
          if len(data) > 10:
              if data[0] == "ATOM" and data[4] == "A":
                  if data[5] not in res_ids:
                      seq.append(data[3])
                      res_ids.append(data[5])

    seq = ''.join([three_to_one_leter_name[res] for res in seq])
    print(f">{name}")
    print(seq)

    