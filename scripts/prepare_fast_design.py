#! /usr/bin/env python
import os
import argparse
import glob
from os.path import exists
from subprocess import Popen, PIPE

####################
# This script prepares fast interface design jobs
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Prepare fast inteface design jobs')
parser.add_argument('-i', nargs=1, help='Input silent file to sample from', type=str, required=True)
parser.add_argument('-n', nargs=1, help='Node list to send jobs to, comma-separated', type=str, required=True)
parser.add_argument('-c', nargs=1, help='Number of CPUs per node to work with', type=int, required=True)
args = parser.parse_args()

## prepare splits for predictor jobs

print("Preparign splits for predictor jobs")
os.makedirs('fast_design/predictor_splits', exist_ok=True)
os.chdir('./fast_design/predictor_splits')
process = Popen(["silentsplitshuf", f'../../{args.i[0]}', "1000"])
returncode = process.wait()
os.chdir('../../')

## list the generated splits

os.makedirs('fast_design/predictor', exist_ok=True)
splits = glob.glob('fast_design/predictor_splits/*.silent')
cwd = os.getcwd()
splits = [f"{cwd}/{split}" for split in splits]
with open('fast_design/predictor/predictor_splits.list', 'w') as f:
	for split in splits:
		f.write(f"{split}\n")
print("Wrote splits list into fast_design/predictor/predictor_splits.list")

## write predictor.flag

if exists("../1_target_preparation/patchdock_residues.txt"):
	with open("../1_target_preparation/patchdock_residues.txt", 'r') as inp:
		patchdock_residues = inp.readlines()[0].strip()
else:
	raise Exception("Can not find the ../1_target_preparation/patchdock_residues.txt file!")

with open('fast_design/predictor/predictor.flag', 'w') as f:
	f.write(f"-script_vars patchdock_res={patchdock_residues}\n")
	f.write(f"-script_vars runpsipred_single={os.environ['PSIPRED']}\n")
	f.write("-dunbrack_prob_buried 0.8\n")
	f.write("-dunbrack_prob_nonburied 0.8       # Use a reduced rotamer set to speed up calculations\n")
	f.write("-dunbrack_prob_buried_semi 0.8\n")
	f.write("-dunbrack_prob_nonburied_semi 0.8\n")
print("Wrote fast_design/predictor/predictor.flag")

## write predictor commands

commands = []
with open(f'fast_design/predictor/predictor_commands.list', 'w') as o:
	for i, split in enumerate(splits):
		os.makedirs(f'fast_design/predictor/{i:05d}', exist_ok=True)
		cd_comand = f"cd {cwd}/fast_design/predictor/{i:05d}"
		rosetta_args = [os.environ['ROSETTA'],
					   "-parser:protocol", f"{os.environ['SCRIPTS']}/paper_predictor.xml",
					   "-beta_nov16",
					   "-in:file:silent", split,
					   "-keep_input_scores", "False",
					   "-silent_read_through_errors",
					   "-out:file:score_only", "score.sc",
					   "-mute", "protocols.rosetta_scripts.ParsedProtocol.REPORT",
					   "-parser:script_vars", f"CAO_2021_PROTOCOL={os.environ['SCRIPTS']}",
					   f"@{cwd}/fast_design/predictor/predictor.flag",
					   "-mute", "all"]
		rosetta_args = ' '.join(rosetta_args)
		commands.append([cd_comand, rosetta_args])

		o.write(f"{cd_comand}\n")
		o.write(f"{rosetta_args}\n\n")
print("Wrote commands for predictor into fast_design/predictor/predictor_commands.list\n")

## prepare splits for pilot design jobs

print("Preparing splits for pilot design job")
os.makedirs('fast_design/design_splits', exist_ok=True)
os.chdir('./fast_design')

# slice out 1000 random docks
silent_for_pilot_design = open(f"for_pilot_design.silent", 'w')
p1 = Popen(["silentls", f'../{args.i[0]}'], stdout=PIPE)
p2 = Popen(["shuf"], stdin=p1.stdout, stdout=PIPE)
p3 = Popen(["head", "-n", "1000"], stdin=p2.stdout, stdout=PIPE)
p4 = Popen(["silentslice", f'../{args.i[0]}'], stdin=p3.stdout, stdout=silent_for_pilot_design)
returncode = p4.wait()
silent_for_pilot_design.close()

# prepare splits
os.chdir('./design_splits')
process = Popen(["silentsplitshuf", '../for_pilot_design.silent', "20"])
returncode = process.wait()
os.chdir('../../')

## list prepared design splits

os.makedirs('fast_design/pilot_design', exist_ok=True)
splits = glob.glob('fast_design/design_splits/*.silent')
splits = [f"{cwd}/{split}" for split in splits]
with open('fast_design/pilot_design/pilot_design_splits.list', 'w') as f:
	for split in splits:
		f.write(f"{split}\n")
print("Wrote splits list into fast_design/pilot_design/pilot_design_splits.list")

## write design.flag

with open(f"fast_design/pilot_design/design.flags", 'w') as f:
	f.write(f"-script_vars patchdock_res={patchdock_residues}\n")
	f.write(f"-script_vars runpsipred_single={os.environ['PSIPRED']}\n")
	f.write(f"-dalphaball {os.environ['DALPAHBALL']}\n")
	f.write(f"-indexed_structure_store:fragment_store {os.environ['SS_GROUPED_VALL_ALL']}\n")
print("Wrote fast_design/pilot_design/design.flags")

## write pilot design commands

with open(f'fast_design/pilot_design/design_commands.list', 'w') as o:
	for i, split in enumerate(splits):
		os.makedirs(f'fast_design/pilot_design/{i:05d}', exist_ok=True)
		cd_comand = f"cd {cwd}/fast_design/pilot_design/{i:05d}"
		rosetta_args = [os.environ['ROSETTA'],
					   "-parser:protocol", f"{os.environ['SCRIPTS']}/paper_interface_design.xml",
					   "-beta_nov16",
					   "-in:file:silent", split,
					   "-keep_input_scores", "False",
					   "-silent_read_through_errors",
					   "-out:file:silent", "out.silent",
					   "-out:file:silent_struct_type", "binary",
					   "-mute", "protocols.rosetta_scripts.ParsedProtocol.REPORT",
					   "-parser:script_vars", f"CAO_2021_PROTOCOL={os.environ['SCRIPTS']}",
					   f"@{cwd}/fast_design/pilot_design/design.flags",
					   "-mute", "all"]
		rosetta_args = ' '.join(rosetta_args)
		commands.append([cd_comand, rosetta_args])

		o.write(f"{cd_comand}\n")
		o.write(f"{rosetta_args}\n\n")
print("Wrote commands for pilot design into fast_design/pilot_design/design_commands.list\n")

## write predictor and pilot design jobs

os.makedirs('fast_design/jobs', exist_ok=True)
nodes = args.n[0].split(',')
for i, node in enumerate(nodes):
	for j in range(args.c[0]):
		with open(f"fast_design/jobs/fast_design_job{i}_{j}.sh", 'w') as f:
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
				with open(f"fast_design/jobs/fast_design_job{i}_{j}.sh", 'a') as f:
					cd_com, rosetta_com = commands[k]
					f.write(f"{cd_com}\n")
					f.write("\n")
					f.write(f"{rosetta_com}\n")
					f.write("\n")
					k += 1

## prepare a script to launch all jobs

with open(f'fast_design/run_fast_design_jobs.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			f.write(f'sbatch jobs/fast_design_job{i}_{j}.sh\n')