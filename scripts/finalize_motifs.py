#! /usr/bin/env python
import os
import argparse
import glob
from subprocess import Popen, PIPE

####################
# This script runs motif finalization before grafting
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Finalize motifs before grafting')
parser.add_argument('-n', nargs=1, help='The number of motifs to pick', type=str, required=True)
args = parser.parse_args()

## pick best motifs

print(f"Selecting {args.n[0]} best motifs...")
os.chdir('./motif_extraction')
selected_motifs = open('selected_motifs.list', 'w')
process = Popen([f"{os.environ['SCRIPTS']}/motif_tools/motif_selection.py",
				 "motif_clusters/cluster_results.list",
  				 "raw_motifs/*/*.json",
  				 "-num_motifs", args.n[0]], stdout=selected_motifs)
returncode = process.wait()
selected_motifs.close()
print("Listed selected motifs at motif_extraction/selected_motifs.list\n")

## Trim down the selected motifs

os.makedirs('selected_motifs', exist_ok=True)
os.makedirs('selected_motifs_trimmed', exist_ok=True)

print(f"Finalizing the selected motifs...")
process = Popen([f"{os.environ['SCRIPTS']}/motif_tools/motif_finalization.py",
				 "all_motifs.list",
  				 "selected_motifs.list",
  				 "selected_motifs_trimmed", "selected_motifs"])
returncode = process.wait()
os.chdir('../')
print("Listed final motifs at motif_extraction/motifs_with_hotspots.list")


