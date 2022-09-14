#! /usr/bin/env python
import os
import argparse

##########################
# This script prepares grafting jobs for obtained motifs
########################## 

# parse arguments for the job
parser = argparse.ArgumentParser(description='Prepare grafting jobs')
parser.add_argument('-i', nargs=1, help='Scaffolds splits list file', type=str, required=True)
parser.add_argument('-n', nargs=1, help='Node list to send jobs to, comma-separated', type=str, required=True)
parser.add_argument('-c', nargs=1, help='Number of CPUs per node to work with', type=int, required=True)
args = parser.parse_args()

## read scaffold splits

scaffold_splits = []
with open(args.i[0], 'r') as f:
	scaffold_splits =[line.strip() for line in f if len(line)>0]
print(f"Found {len(scaffold_splits)} scaffold splits")

## read motifs with hotspots

motifs_with_hotspots = []
with open('motif_extraction/motifs_with_hotspots.list', 'r') as f:
	motifs_with_hotspots =[line.strip().split(' ') for line in f if len(line)>0]
print(f"Found {len(motifs_with_hotspots)} motifs")

## create dirs and commands for grafting

os.makedirs(f'grafting', exist_ok=True)

commands = []
cwd = os.getcwd()
for motif in motifs_with_hotspots:
	motif_name = motif[0].split('/')[-1][:-4]
	os.makedirs(f'grafting/{motif_name}')

	for i, split in enumerate(scaffold_splits):
		os.makedirs(f'grafting/{motif_name}/{i:05d}')

		# create commands
		motif_graft_command = [os.environ['ROSETTA'],
			  				   "-parser:protocol", f"{os.environ['SCRIPTS']}/paper_motif_graft.xml",
							   "-beta_nov16",
							   "-in:file:silent", split,
							   "-keep_input_scores", "False",
							   "-silent_read_through_errors",
							   "-out:file:silent", "out.silent",
							   "-out:file:silent_struct_type", "binary",
							   "-mute", "protocols.rosetta_scripts.ParsedProtocol.REPORT",
							   "-parser:script_vars", f"CAO_2021_PROTOCOL={os.environ['SCRIPTS']}",
							   f"@{cwd}/grafting/grafting.flags",
							   "-parser:script_vars", f"motifpdb={motif[0]}",
							   f"hotspots={motif[1]}",
							   "-out:prefix", f"{motif_name}_",
							   "-mute", "all"]

		commands.append([f"cd {cwd}/grafting/{motif_name}/{i:05d}",
						 ' '.join(motif_graft_command)])

## write commands into separate file

with open(f'grafting/grafting_commands.list', 'w') as f:
	for command in commands:
		f.write(f"{command[0]}\n")
		f.write(f"{command[1]}\n")
		f.write("\n")

## write grafting.flags file

with open(f'grafting/grafting.flags', 'w') as f:
	f.write('-jd2::failed_job_exception 0\n')
	f.write('-parse_script_once_only\n')
	f.write(f"-parser:script_vars contextpdb={'/'.join(os.environ['SCRIPTS'].split('/')[:-1])}/2_docking/receptor_centered_chainchanged.pdb\n")
	f.write('-parser:script_vars initial_sasa_threshold=1400\n')
	f.write('-run::crash_to_console false\n')

#№ write sbatch job scripts

os.makedirs(f'grafting/jobs', exist_ok=True)
nodes = args.n[0].split(',')
for i, node in enumerate(nodes):
	for j in range(args.c[0]):
		with open(f"grafting/jobs/grafting_job{i}_{j}.sh", 'w') as f:
		    f.write('#!/bin/bash\n')
		    f.write('\n')
		    f.write(f'#SBATCH --job-name="grafting{i}_{j}"\n')
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
				with open(f"grafting/jobs/grafting_job{i}_{j}.sh", 'a') as f:
					cd_com, rosetta_com = commands[k]
					f.write(f"{cd_com}\n")
					f.write("\n")
					f.write(f"{rosetta_com}\n")
					f.write("\n")
					k += 1

## prepare a script to launch all jobs

with open(f'grafting/run_grafting_jobs.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			f.write(f'sbatch jobs/grafting_job{i}_{j}.sh\n')