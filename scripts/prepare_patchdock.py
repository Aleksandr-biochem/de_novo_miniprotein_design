#! /usr/bin/env python
import os
from os.path import exists
import argparse
from subprocess import Popen

####################
# This script prepares PatchDock runs to generate seeds for RifDock
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Run PatchDock')
parser.add_argument('-i', nargs=1, help='The list of input poly-Valine scaffolds', type=str, required=True)
parser.add_argument('-n', nargs=1, help='Node list to send jobs to, comma-separated', type=str, required=True)
parser.add_argument('-c', nargs=1, help='Number of CPUs per node to work with', type=int, required=True)
args = parser.parse_args()

## prepare pathdock residues list

os.makedirs("patchdock_xforms", exist_ok=True)
if exists("../1_target_preparation/patchdock_residues.txt"):
	with open('patchdock_xforms/patchdock_residues.list', 'w') as f:
		with open("../1_target_preparation/patchdock_residues.txt", 'r') as inp:
			residues = inp.readlines()[0].strip().split(',')
		if len(residues) > 9:
			print(f"Warning! The residues list for PatchDock contains {len(residues)} residues. You may want to reduce the list below 9-10.")
		for res in residues:
			f.write(f"{res} B\n")
else:
	raise Exception("Can not find the patchdock_residues.txt file!")

## read poly-Valine scaffolds to dock

with open(args.i[0], 'r') as f:
	scaffolds = [line.strip() for line in f]

## prepare PatchDock commands to launch

os.makedirs("patchdock_xforms/params", exist_ok=True)
cwd = os.getcwd()

commands = []
for scaffold in scaffolds:
	name = scaffold.split('/')[-1][:-4]

	# write params file
	with open(f"patchdock_xforms/params/{name}.params", "w") as f:
		f.write(f"receptorPdb {cwd}/receptor_centered_chainchanged.pdb\n")
		f.write(f"ligandPdb {scaffold}\n")
		f.write(f"protLib {'/'.join(os.environ['PATCHDOCK'].split('/')[:-1])}/chem.lib\n")
		f.write("log-file patchdock_log.log\n")
		f.write("log-level 0\n")
		f.write("receptorSeg 10.0 20.0 1.5 1 0 1 0\n")
		f.write("ligandSeg 10.0 20.0 1.5 1 0 1 0\n")
		f.write("scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0\n")
		f.write("desolvationParams 500.0 1.0\n")
		f.write("clusterParams 0.1 4 2.0 3.000000\n")
		f.write("baseParams 4.0 13.0 2\n")
		f.write("matchingParams 1.5 1.5 0.4 0.5 0.9\n")
		f.write("matchAlgorithm 1\n")
		f.write("receptorGrid 0.5 6.0 6.0\n")
		f.write("ligandGrid 0.5 6.0 6.0\n")
		f.write("receptorMs 10.0 1.8\n")
		f.write("ligandMs 10.0 1.8\n")
		f.write(f"receptorActiveSite {cwd}/patchdock_xforms/patchdock_residues.list\n")

	# appends command
	commands.append(f"{os.environ['PATCHDOCK']} {cwd}/patchdock_xforms/params/{name}.params {cwd}/patchdock_xforms/patchdock_out/{name}.out")

## write commands into a separate file

with open(f'patchdock_xforms/patchdock_commands.list', 'w') as f:
	for command in commands:
		f.write(f"{command}\n")

# write sbatch job scripts
os.makedirs("patchdock_xforms/jobs", exist_ok=True)
nodes = args.n[0].split(',')
for i, node in enumerate(nodes):
	for j in range(args.c[0]):
		with open(f"patchdock_xforms/jobs/patchdock_job{i}_{j}.sh", 'w') as f:
		    f.write('#!/bin/bash\n')
		    f.write('\n')
		    f.write(f'#SBATCH --job-name="patchdock{i}_{j}"\n')
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
				with open(f"patchdock_xforms/jobs/patchdock_job{i}_{j}.sh", 'a') as f:
					f.write(f"{commands[k]}\n")
					f.write("\n")
					k += 1

# prepare a script to launch all jobs
with open(f'patchdock_xforms/run_patchdock_jobs.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			f.write(f'sbatch jobs/patchdock_job{i}_{j}.sh\n')

# create a directory for out files
os.makedirs("patchdock_xforms/patchdock_out", exist_ok=True)
