#! /usr/bin/env python
import os
import argparse
import glob
from subprocess import Popen, PIPE

####################
# This script runs motif extraction and clustering
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Prepare raw motif and cluster them:')
parser.add_argument('-i', nargs=1, help='Input silent file to extract motifs from', type=str, required=True)
parser.add_argument('--dump_og', help='Dump pdb of the motif in presence of the target.', action='store_true')
parser.add_argument('--less_motifs', help='More stringent motif extraction criteria. Will produce less motifs. ', action='store_true')
args = parser.parse_args()

## prepare splits for motif extraction

print("Preparing splits for predictor jobs")
os.makedirs('motif_extraction/designs_splits', exist_ok=True)
os.chdir('./motif_extraction/designs_splits')
process = Popen(["silentsplitshuf", f'../../{args.i[0]}', "1000"])
returncode = process.wait()
os.chdir('../../')

## list the generated splits

splits = glob.glob('motif_extraction/designs_splits/*.silent')
cwd = os.getcwd()
splits = [f"{cwd}/{split}" for split in splits]
with open('motif_extraction/designs_splits.list', 'w') as f:
	for split in splits:
		f.write(f"{split}\n")
print("Wrote splits list into motif_extraction/designs_splits.list\n")

## run raw motif extraction

os.makedirs('motif_extraction/raw_motifs', exist_ok=True)
print("Starting the motif extraction...")
for i, split in enumerate(splits):
	os.makedirs(f'motif_extraction/raw_motifs/{i:05d}', exist_ok=True)
	os.chdir(f'motif_extraction/raw_motifs/{i:05d}')

	motif_extraction_args = [f"{os.environ['SCRIPTS']}/motif_tools/motif_extraction.py",
							 "-in:file:silent", split,
							 "-ref_pdb", f"{'/'.join(os.environ['SCRIPTS'].split('/')[:-1])}/2_docking/receptor_centered_chainchanged.pdb",
	    					 "-out_prefix", "mot_"]

	if args.dump_og:
		motif_extraction_args.append("-dump_og")

	if args.less_motifs:
		motif_extraction_args.append("-ddg_threshold")
		motif_extraction_args.append("-25")

	log = open('log.log', 'w')
	process = Popen(motif_extraction_args, stdout=log)
	returncode = process.wait()
	log.close()
	print(f"Processed split {i + 1} out of {len(splits)}")
	os.chdir(cwd)

## list extracted motifs

raw_motifs = glob.glob('motif_extraction/raw_motifs/*/*.gz')
raw_motifs = [motif for motif in raw_motifs if '_og.pdb.gz' not in motif]
print(f"\nFound {len(raw_motifs)} raw motifs")

with open('motif_extraction/all_motifs.list', 'w') as f:
	for motif in raw_motifs:
		f.write(f"{cwd}/{motif}\n")
print("Wrote the list of extracted motifs into motif_extraction/all_motifs.list\n")

## cluster motifs

os.makedirs('motif_extraction/motif_clusters', exist_ok=True)
os.chdir('motif_extraction/motif_clusters')

print("Launched motif clustering...")
process = Popen([os.environ['CLUSTER'], '../all_motifs.list', '0.7'])
returncode = process.wait()

with open('cluster_results.list', 'r') as f:
	num_clusters = sum(1 for line in f if line.startswith('Cluster:'))
os.chdir('../../')
print(f"Generated {num_clusters} motif clusters. Listed them at motif_extraction/motif_clusters/cluster_results.list")


