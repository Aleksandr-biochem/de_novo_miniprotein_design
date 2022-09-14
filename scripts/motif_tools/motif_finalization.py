#!/usr/bin/env python

import os
import sys

from pyrosetta import *
from pyrosetta.rosetta import *
init()


import numpy as np

all_motifs_list = sys.argv[1]
longxing_list = sys.argv[2]
motif_folder = sys.argv[3]
og_folder = sys.argv[4]

assert( os.path.exists( motif_folder ) )
assert( os.path.exists( og_folder ) )

tag_to_path = {}
with open(all_motifs_list) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0):
            continue
        tag = os.path.basename(line)
        if ( tag.endswith(".gz") ):
            tag = tag[:-len(".gz")]
        if ( tag.endswith(".pdb") ):
            tag = tag[:-len(".pdb")]
        if ( tag.endswith("_og")):
            continue

        tag_to_path[tag] = line

longxing_splits = []
with open(longxing_list) as f:
    for line in f:
        line = line.strip()
        if ( len(line) == 0 ):
            continue
        longxing_splits.append(line.split())


def delete_residues_smart(pose, start, end):

    if ( start > end ):
        return
    ft = pose.fold_tree()
    if ( start == 1 ):
        ft.reorder(pose.size())
    else:
        ft.reorder(1)
    pose.fold_tree(ft)

    pose.delete_residue_range_slow(start, end)
    pose.conformation().detect_disulfides()


def trim_motif(pose, start, end):

    motif_size = pose.conformation().chain_end(1)

    delete_residues_smart(pose, end+1, motif_size)
    delete_residues_smart(pose, 1, start-1)

    pose.pdb_info( core.pose.PDBInfo( pose ) )

with open("motifs_with_hotspots.list", "w") as f:

    for tag, _, _, _, start, end, hotspots in longxing_splits:
        start = int(start)
        end = int(end)

        path = tag_to_path[tag]
        og_path = path.replace(".pdb", "_og.pdb")

        assert(os.path.exists(path))
        #assert(os.path.exists(og_path))

        pose = pose_from_file(path)
        trim_motif(pose, start, end)
        motif_final_path = os.path.join(motif_folder, os.path.basename(path))
        pose.dump_pdb(motif_final_path)

        if ( os.path.exists(og_path) ):
            pose = pose_from_file(og_path)
            trim_motif(pose, start, end)
            pose.dump_pdb(os.path.join(og_folder, os.path.basename(og_path)))

        f.write("%s %s\n"%(os.path.abspath(motif_final_path), hotspots.replace(',', ':')))




