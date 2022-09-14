#! /usr/bin/env python
import os
import glob
import argparse
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description='Prepare scaffolds for future use')
parser.add_argument('-i', nargs=1, help='Scaffold dir(s) to list and prepare, comma-separated if many', type=str, required=True)
args = parser.parse_args()

scaffolds_to_list = args.i[0].split(',')

# list all pdb scaffolds
print(f"Listing scaffolds from {', '.join(scaffolds_to_list)}")
scaffolds_counts = dict()
total = 0
cwd = os.getcwd()
with open(f"scaffolds_{'_'.join(scaffolds_to_list)}.list", 'w') as f:
	for batch in scaffolds_to_list:
		scaffolds = glob.glob(f"{batch}/*.pdb")
		total += len(scaffolds)
		scaffolds_counts[batch] = len(scaffolds)
		for scaf in scaffolds:
			f.write(f"{cwd}/{scaf}\n")

print(f"Found {total} scaffolds")
for batch in scaffolds_counts:
	print(f"{scaffolds_counts[batch]} scaffolds in {batch}")
print()

# convert all scaffolds to silent
print("Converting scaffolds to silent")
scaffolds_silent = open(f"scaffolds_{'_'.join(scaffolds_to_list)}.silent", 'w')
p1 = Popen(["cat", f"scaffolds_{'_'.join(scaffolds_to_list)}.list"], stdout=PIPE)
p2 = Popen(["silentfrompdbs"], stdin=p1.stdout, stdout=scaffolds_silent)
returncode = p2.wait()
scaffolds_silent.close()

# prepare scaffolds splits 
print("Preparing splits for grafting")
os.makedirs(f"{'_'.join(scaffolds_to_list)}_splits")
os.chdir(f"{'_'.join(scaffolds_to_list)}_splits")
process = Popen(["silentsplitshuf", f"../scaffolds_{'_'.join(scaffolds_to_list)}.silent", "1000"])
returncode = process.wait()
os.chdir(f"../")

# list all splits
print("Listing scaffolds splits")
scaffold_splits = glob.glob(f"{'_'.join(scaffolds_to_list)}_splits/*.silent")
print(f"Found {len(scaffold_splits)}")
with open(f"scaffolds_{'_'.join(scaffolds_to_list)}_split.list", 'w') as f:
	for split in scaffold_splits:
		f.write(f"{cwd}/{split}")
print(f"Wrote splits into scaffolds_{'_'.join(scaffolds_to_list)}_split.list")
