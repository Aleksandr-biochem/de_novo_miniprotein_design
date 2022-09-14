#! /usr/bin/env python
import os
import argparse
import glob
from os.path import exists
from subprocess import Popen, PIPE

##############################
# This script prepares interface design jobs
##############################

parser = argparse.ArgumentParser(description='Prepare interface design jobs')
parser.add_argument('-i', nargs=1, help='A silent file with structures for design', type=str, required=True)
parser.add_argument('-n', nargs=1, help='Node list to send jobs to, comma-separated', type=str, required=True)
parser.add_argument('-c', nargs=1, help='Number of CPUs per node to work with', type=int, required=True)
args = parser.parse_args()

## prepare splits

print("Preparing splits for design...")
os.makedirs("design_splits", exist_ok=True)
os.chdir("./design_splits")
process = Popen(["silentsplitshuf", f'../{args.i[0]}', "20"])
returncode = process.wait()
os.chdir(f"../")

## list splits

splits = glob.glob(f"design_splits/*.silent")
splits = [f"{os.getcwd()}/{split}" for split in splits]
print(f"There are {len(splits)} splits of structures for design")
with open(f"design_split.list", 'w') as f:
	for split in splits:
		f.write(f"{split}\n")
print("Wrote splits list into design_split.list")

## write design.flag

if exists("../1_target_preparation/patchdock_residues.txt"):
	with open("../1_target_preparation/patchdock_residues.txt", 'r') as inp:
		patchdock_residues = inp.readlines()[0].strip()
elif exists("../../1_target_preparation/patchdock_residues.txt"):
	with open("../../1_target_preparation/patchdock_residues.txt", 'r') as inp:
		patchdock_residues = inp.readlines()[0].strip()
else:
	raise Exception("Can not find the ../(../)1_target_preparation/patchdock_residues.txt file!")

with open("design.flags", 'w') as f:
	f.write(f"-script_vars patchdock_res={patchdock_residues}\n")
	f.write(f"-script_vars runpsipred_single={os.environ['PSIPRED']}\n")
	f.write(f"-dalphaball {os.environ['DALPAHBALL']}\n")
	f.write(f"-indexed_structure_store:fragment_store {os.environ['SS_GROUPED_VALL_ALL']}\n")
print("Wrote design.flags\n")

## write subdirs and commands

print("Creating design subdirs and setting up commands...")
commands = []
for i, split in enumerate(splits):
	os.makedirs(f"designs/{i:05d}", exist_ok=True)

	design_args = [os.environ['ROSETTA'],
				   "-parser:protocol", f"{os.environ['SCRIPTS']}/paper_interface_design.xml",
				   "-beta_nov16",
				   "-in:file:silent", split,
				   "-keep_input_scores", "False",
				   "-silent_read_through_errors",
				   "-out:file:silent", "out.silent",
				   "-out:file:silent_struct_type", "binary",
				   "-mute", "protocols.rosetta_scripts.ParsedProtocol.REPORT",
				   "-parser:script_vars", f"CAO_2021_PROTOCOL={os.environ['SCRIPTS']}",
				   f"@{os.getcwd()}/design.flags",
				   "-mute", "all"]

	commands.append([f"cd {os.getcwd()}/designs/{i:05d}",
					 ' '.join(design_args)])

## write sbatch job scripts

os.makedirs('jobs', exist_ok=True)
nodes = args.n[0].split(',')
for i, node in enumerate(nodes):
	for j in range(args.c[0]):
		with open(f"jobs/design_job{i}_{j}.sh", 'w') as f:
		    f.write('#!/bin/bash\n')
		    f.write('\n')
		    f.write(f'#SBATCH --job-name="design{i}_{j}"\n')
		    f.write('#SBATCH --ntasks=1\n')
		    f.write('#SBATCH --ntasks-per-core=1\n')
		    f.write('#SBATCH --partition=cpu\n')
		    f.write(f'#SBATCH --nodelist={node}\n')
		    f.write('\n')

k = 0
while k < len(commands):
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			if k == len(commands):
				break
			else:
				with open(f"jobs/design_job{i}_{j}.sh", 'a') as f:
					cd_com, rosetta_com = commands[k]
					f.write(f"{cd_com}\n")
					f.write("\n")
					f.write(f"{rosetta_com}\n")
					f.write("\n")
					k += 1

## prepare a script to launch all jobs

with open(f'run_design_jobs.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			f.write(f'sbatch jobs/design_job{i}_{j}.sh\n')

