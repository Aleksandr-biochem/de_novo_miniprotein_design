#! /usr/bin/env python
import os
import glob
import argparse
from subprocess import Popen, PIPE

#####################
# This script prepares poly-Valine scaffolds for PatchDock docking
#####################

parser = argparse.ArgumentParser(description='Prepare pily-Valine scaffolds for future use')
parser.add_argument('-i', nargs=1, help='Scaffold dir(s) to list and prepare, comma-separated if many', type=str, required=True)
args = parser.parse_args()

scaffolds_to_list = args.i[0].split(',')

os.makedirs(f"{'_'.join(scaffolds_to_list)}_polyV", exist_ok=True)

## list all pdb scaffold files

print(f"Listing scaffolds from {', '.join(scaffolds_to_list)}")
scaffolds_counts = dict()
total = 0
cwd = os.getcwd()
with open(f"{'_'.join(scaffolds_to_list)}_polyV/original_pdb.list", 'w') as f:
	for batch in scaffolds_to_list:
		scaffolds = glob.glob(f"{batch}/*.pdb")
		total += len(scaffolds)
		scaffolds_counts[batch] = len(scaffolds)
		for scaf in scaffolds:
			f.write(f"{cwd}/{scaf}\n")
print(f"Found {total} scaffolds")
for batch in scaffolds_counts:
	print(f"{scaffolds_counts[batch]} scaffolds in {batch}")

## convert all scaffolds to poly-Valines

print("\nConverting scaffolds to poly-Valines")
os.chdir(f"./{'_'.join(scaffolds_to_list)}_polyV")
log = open("log.log", 'w')
rosetta_flags = [os.environ['ROSETTA'],
 				 "-parser:protocol", f"{os.environ['SCRIPTS']}/polyV.xml",
 				 "-beta_nov16", "-l", "original_pdb.list",
 				 "-mute", "protocols.rosetta_scripts.ParsedProtocol.REPORT",
 				 "-parser:script_vars", f"CAO_2021_PROTOCOL={os.environ['SCRIPTS']}"]
p = Popen(rosetta_flags, stdout=log, stderr=log)
returncode = p.wait()
log.close()

## list poly-Valine scaffolds

plyV_scaffolds = glob.glob('*.pdb')
cwd = os.getcwd()
with open(f"../{'_'.join(scaffolds_to_list)}_polyV.list", 'w') as f:
	for scaffold in plyV_scaffolds:
		f.write(f"{cwd}/{scaffold}\n")
