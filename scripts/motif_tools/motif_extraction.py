#!/usr/bin/env python

import os
import sys
import json
import argparse
import numpy as np
import itertools


from pyrosetta import *
from pyrosetta.rosetta import *

init("-corrections::beta_nov16 -use_truncated_termini -in:file:silent_struct_type binary")

scorefxn = get_fa_scorefxn()


def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose

abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))


def better_dssp_0(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = dssp.get_dssp_secstruct()

    my_dssp = ""

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos-1]
        if ( the_dssp[seqpos-1] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case
            if ( abego == "B" ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser()
parser.add_argument("-ref_pdb", type=str, default="", help="the reference pdb file")
parser.add_argument("-ddg_threshold", type=float, default=-20, help="the ddg cutoff value of the motif")
parser.add_argument("-out_prefix", type=str, default="", help="prefix on out files")
parser.add_argument("-dump_og", action="store_true", help="Output full motif complex")
parser.add_argument("pdbs", type=str, nargs="*", help="Inputs")
parser.add_argument("-pdb_list", type=str, default="", help="")
parser.add_argument("-in:file:silent", type=str, default="")

args = parser.parse_args(sys.argv[1:])

silent = args.__getattribute__("in:file:silent")



def get_my_interface_graph(pose):
    energy_graph = pose.energies().energy_graph()

    graph = np.zeros((pose.size()+1, pose.size()+1), np.float)

    weights = scorefxn.weights()

    it = energy_graph.const_edge_list_begin()
    while it.valid():
        edge = it.__mul__()
        it.plus_plus()
        assert( not edge is None )

        low_node = edge.get_first_node_ind()
        high_node = edge.get_second_node_ind()
        if ( pose.chain(low_node) == pose.chain(high_node ) ):
            continue

        assert( low_node < high_node )

        score = edge.dot( weights )

        graph[low_node, high_node] = score

    return graph



# align two things using sequence and accepting gaps
def pymol_align( move_pose, to_pose, sel_move=None, sel_to=None, atoms=["N", "CA", "C"], throw_away=0.1 ):

    if ( not sel_move is None ):
        move_res = np.array(list(core.select.get_residues_from_subset(sel_move)))
    else:
        move_res = np.array(list(range(1, move_pose.size()+1)))

    if ( not sel_to is None ):
        to_res = np.array(list(core.select.get_residues_from_subset(sel_to)))
    else:
        to_res = np.array(list(range(1, to_pose.size()+1)))

    seq_move = "x" + move_pose.sequence()
    seq_to = "x" + to_pose.sequence()

    seq_move = "".join(np.array(list(seq_move))[move_res])
    seq_to = "".join(np.array(list(seq_to))[to_res])

    from Bio import pairwise2
    alignment = align_move, align_to, idk1, idk2, idk3 = pairwise2.align.globalxs(seq_move,seq_to, -2, -1)[0]
    # print(align_move, align_to)

    move_to_pairs = []
    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()

    i_move = 0
    i_to = 0
    for i in range(len(align_move)):
        if ( align_move[i] == align_to[i] ):

            seqpos_move = move_res[i_move]
            seqpos_to = to_res[i_to]

            move_to_pairs.append((seqpos_move, seqpos_to))

            for atom in atoms:
                coords_move.append(move_pose.residue(seqpos_move).xyz(atom))
                coords_to.append(to_pose.residue(seqpos_to).xyz(atom))


        if ( align_move[i] != "-" ):
            i_move += 1
        if ( align_to[i] != "-" ):
            i_to += 1

    move_pose_copy = move_pose.clone()

    rmsd = 0

    distances = []

    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )
        
        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance
                distances.append(distance)

        rmsd /= len(move_to_pairs)*len(atoms)
        rmsd = np.sqrt(rmsd)

    move_pose = move_pose_copy

    distances = np.array(distances)

    cutoff = np.percentile(distances, 90)

    mask = distances <= cutoff

    print(mask.sum(), len(mask))

    coords_move_old = list(coords_move)
    coords_to_old = list(coords_to)

    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()

    for i in range(len(coords_move_old)):
        if ( not mask[i] ):
            continue
        coords_move.append(coords_move_old[i])
        coords_to.append(coords_to_old[i])

    rmsd = 0
    imask = -1
    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )
        
        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                imask += 1
                if ( not mask[imask] ):
                    continue
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance

        rmsd /= imask
        rmsd = np.sqrt(rmsd)

    return rmsd, move_to_pairs, move_pose


rmsd_sel = core.select.residue_selector.ChainSelector("B")

ref_fname = args.ref_pdb
save_ref_pose = pose_from_file(ref_fname)

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    fnames = sfd_in.tags()
elif (args.pdb_list != "" ):
    fnames = []
    with open(args.pdb_list) as f:
        for line in f:
            line = line.strip()
            if ( len(line) == 0 ):
                continue
            fnames.append(line)

else:
    fnames = args.pdbs


for ifname, fname in enumerate(fnames):

    # try:
    for k in [1]:
        basename = os.path.basename(fname)

        ftag = ""
        if (fname.endswith(".pdb")):
            ftag = basename[:-4]
        if (fname.endswith(".pdb.gz")):
            ftag = basename[:-7]
        if ( not ftag ):
            ftag = basename
        assert(ftag)
        ftag = args.out_prefix + ftag

        print("Processing pdb file %s"%fname )

        ddg_threshold = args.ddg_threshold

        if ( silent == "" ):
            pose = pose_from_file(fname)
        else:
            pose = Pose()
            sfd_in.get_structure(fname).fill_pose(pose)

        ref_pose = save_ref_pose.clone()

        rmsd, _, pose = pymol_align(pose, ref_pose, rmsd_sel.apply(pose), rmsd_sel.apply(ref_pose))

        print("RMSD: %.3f"%rmsd)


        binder_target = pose.split_by_chain()
        binder = binder_target[1]
        target = binder_target[2]
        chainA_len = binder.size()
        print("length of the binder %i"%chainA_len)


####################################################################################
# Per residue ddg

        score = scorefxn(pose)
        print("total score of the complex: %.3f"%score)

        graph = get_my_interface_graph(pose)


        separate_pose = move_chainA_far_away(pose)
        scorefxn(separate_pose)

        per_res_ddg = utility.vector1_double()

        for i in range(1, pose.size()+1):
            ddg = np.sum( graph[:,i]) + np.sum( graph[i,:] )
            per_res_ddg.append(ddg)

        per_res_ddg0 = list(per_res_ddg)
        per_res_ddg = [0] + list(per_res_ddg)

############################################
# Per residue ddg when compared to glycine

        all_positions = utility.vector1_unsigned_long()
        for i in range(1, binder.size()+1):
            all_positions.append(i)

        poly_ala_binder_pose = pose.clone()
        protocols.toolbox.pose_manipulation.construct_poly_XXX_pose("GLY", poly_ala_binder_pose, all_positions, True, True, True)
        scorefxn(poly_ala_binder_pose)

        ala_graph = get_my_interface_graph(poly_ala_binder_pose)

        poly_ala_binder_separate_pose = move_chainA_far_away(poly_ala_binder_pose)
        scorefxn(poly_ala_binder_separate_pose)

        per_res_ddg_ala = utility.vector1_double()
        per_res_ala_scan = utility.vector1_double()
        for i in range(1, pose.size()+1):
            ddg = np.sum(ala_graph[:,i]) + np.sum(ala_graph[i,:])
            per_res_ddg_ala.append(ddg)
            per_res_ala_scan.append(per_res_ddg[i] - ddg)

        per_res_ddg0_ala = list(per_res_ddg_ala)

######################################
# Extract H and E motifs

        dssp_str = better_dssp_0(binder)

        if ( dssp_str[0] == "L" and dssp_str[1] != "L"):
            tmp = list(dssp_str)
            tmp[0] = tmp[1]
            dssp_str = ''.join(tmp) 
        if (dssp_str[-1] == "L" and dssp_str[-2] != "L"):
            tmp = list(dssp_str)
            tmp[-1] = tmp[-2]
            dssp_str = ''.join(tmp) 

        segs = []
        seg_temp = {}

        for ii in range(len(dssp_str)):
            if (ii == 0):
                seg_temp["start"] = 1
                seg_temp["sec_type"] = dssp_str[ii]

            if ( dssp_str[ii] != seg_temp["sec_type"]):
                seg_temp["end"] = ii
                segs.append(seg_temp)
                seg_temp = {}
                seg_temp["start"] = ii + 1
                seg_temp["sec_type"] = dssp_str[ii]

            if ( ii == len(dssp_str)-1 ):
                seg_temp["end"] = ii + 1
                segs.append(seg_temp)

##########################################################
# Filter motifs by ddg and output the good ones

        for seg in segs:
            seg["ddg"] = 0
            for ires in range(seg["start"], seg["end"]+1):    
                seg["ddg"] += per_res_ddg[ires]

        good_motifs = []

        for seg in segs:
            if ( seg["ddg"] > ddg_threshold ):
                print("Reject ddg")
                continue

            motif_res = utility.vector1_unsigned_long()

            for i in range(seg["start"], seg["end"] + 1 ):
                motif_res.append(i)


            motif_alone = core.pose.Pose()
            core.pose.pdbslice(motif_alone, binder, motif_res)

            motif_alone.pdb_info(core.pose.PDBInfo(motif_alone))
            for i in range(1, motif_alone.size()+1):
                motif_alone.pdb_info().chain(i, "A")

            motif_complex = motif_alone.clone()
            motif_complex.append_pose_by_jump(target, 1)

            motif_complex.pdb_info(core.pose.PDBInfo(motif_complex))
            for i in range(1, motif_complex.size()+1):
                if ( i <= motif_alone.size() ):
                    motif_complex.pdb_info().chain(i, "A")
                else:
                    motif_complex.pdb_info().chain(i, "B")

            seg["per_res_ddg"] = per_res_ddg0[seg["start"]-1:seg["end"]+1-1]
            seg["per_res_ddg_ala"] = per_res_ddg0_ala[seg["start"]-1:seg["end"]+1-1]

            fout = ftag + "_%i_%i_%s.pdb.gz"%(seg["start"], seg["end"], seg["sec_type"])

            motif_alone.dump_pdb(fout)
            if ( args.dump_og ):
                motif_complex.dump_pdb(fout.replace(".pdb.gz", "_og.pdb.gz"))

            def default(o):
                if isinstance(o, np.int64): return int(o)  
                raise TypeError

            fname = fout[:-7]
            seg["name"] = fname
            f = open(fname + ".json", "w")
            f.write(json.dumps(seg, indent=4, default=default))
            f.close()






        print("Job done!")

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



