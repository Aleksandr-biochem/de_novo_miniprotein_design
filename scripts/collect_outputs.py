#! /usr/bin/env python
import os
import argparse
import glob
from tqdm import tqdm
from subprocess import Popen, PIPE

####################
# This script collects outputs for predictor, design or grafting jobs
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Collect outputs for predictor, design or grafting jobs:')
parser.add_argument('--predictor', nargs=1, help='Input dir with predictor results', type=str, required=False)
parser.add_argument('--design', nargs=1, help='Input dir with design results', type=str, required=False)
parser.add_argument('--grafting', nargs=1, help='Input dir with grafting results', type=str, required=False)
args = parser.parse_args()

## collect score files from predictor

if not args.predictor is None:
	print("Collecting the predictor results...")
	score_files = glob.glob(f'{args.predictor[0]}/*/*.sc')
	with open('predictor_score_combined.sc', 'w') as f:
		got_header = False
		for file in score_files:
			with open(file, 'r') as sc:
				for line in sc:
					if 'total_score' in line:
                        if not got_header:
                            f.write(line)
                    elif not 'SEQUENCE' in line:
                        f.write(line)
            got_header = True
	print("Saved the predictor_score_combined.sc\n")

## collect silent files from design jobs

if not args.design is None:
	outs = glob.glob(f"{args.design[0]}/*/out.silent")
	with open(f"design_combined.silent", 'w') as combined:
		for out in tqdm(outs, desc="Collecting the design results"):
			with open(out, 'r') as f:
				for line in f:
					combined.write(line)
	print("Combined designs into design_combined.silent")
	process = Popen(["silentscorefile", 'design_combined.silent'])
	returncode = process.wait()
	print("Generated a score file design_combined.sc\n")


## collect silent files from grafting jobs

if not args.grafting is None:
	outs = glob.glob(f"{args.grafting[0]}/*/*/out.silent")
	with open(f"graftings_combined.silent", 'w') as combined:
		for out in tqdm(outs, desc="Collecting the grafting results"):
			with open(out, 'r') as f:
				for line in f:
					combined.write(line)
	# count graftings
	print(f"Combined graftings into graftings_combined.silent")
	print("The number of found graftings:") 
	process1 = Popen(["silentls", "graftings_combined.silent"], stdout=PIPE)
	process2 = Popen(["wc", "-l"], stdin=process1.stdout)
	returncode = process2.wait()
	
