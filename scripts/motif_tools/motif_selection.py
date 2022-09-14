#!/usr/bin/env python
# coding: utf-8

# In[2]:

import sys

import numpy as np
import pandas as pd

import glob
import json
import re
import random
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("CLUSTER_RESULTS", type=str, help="the output from the cluster program")
parser.add_argument("JSON_PATH", type=str, help="a * path that captures all the .json from the motif extraction step")
parser.add_argument("-num_motifs", type=int, default=1200, help="how many motifs to output")


args = parser.parse_args(sys.argv[1:])


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



CLUSTER_RESULTS = args.CLUSTER_RESULTS 
JSON_PATH = args.JSON_PATH 



clusters = defaultdict(dict)
with open(CLUSTER_RESULTS) as f:
    lines = [line.strip() for line in f]
for index, line in enumerate( lines ):
    if index % 3 ==0:
        cluster = int( line.split()[1] )
        statistics = line
        alignments = lines[index+1]
        pdb_names  = lines[index+2]
        clusters[cluster]['size'] = int(line.split()[3])
        clusters[cluster]['alignments'] = [ int(ii) for ii in alignments.split() ]
        clusters[cluster]['pdb_names'] = pdb_names.split()



motifs = glob.glob(JSON_PATH)
js = {}
for i, m in enumerate(motifs):
    if ( i % 20000 == 0 ):
        eprint(i, "/", len(motifs))
    try:
        with open(m) as f:
            j = json.load(f)
    except:
        eprint("Failed to load %s"%m)
        continue
    
    if ( j['sec_type'] == "P" ):
        continue
    j['len'] = j['end'] - j['start'] + 1
    js[j['name']] = j



def select_best_motif( cluster, min_motif_len = 5, verbose = False, bonus_key = None, bonus_value = -5.0 ):
    average_ddg = [0.0] * 100
    min_index = 999
    max_index = -999
    
    #lets set 50 to be the first residue
    for icluster in range(cluster['size']):
        alignment = cluster['alignments'][icluster]
        pdb = cluster['pdb_names'][icluster]
        tag = pdb.split('/')[-1][:-7]
        info = js[tag]
        length = info['len']
        # indexed from 0
        a_shift = alignment // 1000
        b_shift = alignment % 1000
        # indexed from 0
        for ires in range(length):
            new_index = ires - (b_shift-a_shift) + 50
            min_index = min( min_index, new_index-50 )
            max_index = max( max_index, new_index-50 )
            average_ddg[new_index] += info['per_res_ddg'][ires]
            
    for jj in range(min_index, max_index+1):
        if verbose: eprint( "{:<5d}".format(jj), end='')
        average_ddg[jj+50] /= cluster['size']
        if verbose: eprint( '{:<10.2f}'.format(average_ddg[jj+50]), end='' )
        for icluster in range(cluster['size']):
            alignment = cluster['alignments'][icluster]
            pdb = cluster['pdb_names'][icluster]
            tag = pdb.split('/')[-1][:-7]
            info = js[tag]
            length = info['len']
            a_shift = alignment // 1000
            b_shift = alignment % 1000
            cursor = jj - a_shift+b_shift
            if verbose:
                if cursor < 0 or cursor >= length:
                    eprint( "{:3s}".format("*"), end='' )
                else:
                    eprint( "{:3s}".format(info['sequence'][cursor]), end='' )
        if verbose: eprint('')
        
    hotspots_global = []
    motif_range = []
    # this loop will add some pre-required hotspots, like hotspots, cation-pi residues
    for icluster in range(cluster['size']):
        alignment = cluster['alignments'][icluster]
        pdb = cluster['pdb_names'][icluster]
        tag = pdb.split('/')[-1][:-7]
        info = js[tag]
        length = info['len']
        # indexed from 0
        a_shift = alignment // 1000
        b_shift = alignment % 1000
        # indexed from 0
        # OK, I need to make it clear.
        # 
        if bonus_key != None:
            for ires in info[bonus_key]:
                motif_range.append( ires - 1 - ( b_shift-a_shift ) )
            
    
    # define the left_cut and right_cut of this cluster
    ddg_resi_pair = []
    for jj in range(min_index, max_index+1):
        ddg_resi_pair.append( (average_ddg[jj+50], jj) )
    ddg_resi_pair.sort(key=lambda x:x[0])
    motif_range.append(ddg_resi_pair[0][1])
    hotspots_global.append( ddg_resi_pair[0][1] )
    best_index = 1
    left_cut  = min( motif_range )
    right_cut = max( motif_range )
    while right_cut - left_cut + 1 < min_motif_len and best_index < len(ddg_resi_pair):
        motif_range.append(ddg_resi_pair[best_index][1])
        hotspots_global.append( ddg_resi_pair[best_index][1] )
        left_cut  = min( motif_range )
        right_cut = max( motif_range )
        best_index += 1
        
    
    # the weights for each position, to downweight some uncommon hotspots, like wierd trp
    weights = average_ddg[:]
    s = sum( average_ddg[left_cut+50: right_cut+1+50] )
    for ii in range(len(weights)):
        weights[ii] /= s
        
    
    final_results = []
    
    for icluster in range(cluster['size']):
        alignment = cluster['alignments'][icluster]
        pdb = cluster['pdb_names'][icluster]
        tag = pdb.split('/')[-1][:-7]
        info = js[tag]
        length = info['len']
        # indexed from 0
        a_shift = alignment // 1000
        b_shift = alignment % 1000
        
        hotspots_specific = hotspots_global[:]
        if bonus_key != None:
            for igood in info[bonus_key]:
                hotspots_specific.append( igood - 1 - ( b_shift-a_shift ) )
        hotspots_specific = list(set(hotspots_specific))
        # this_should be the index for the new motif, so
        # final hotspots str
        hotspots_new_index = []
        
        # raw_ddg doesn't depends on the hbond_bonus
        raw_ddg = 0.0
        weighted_ddg = 0.0
        unweighted_ddg = 0.0
        if verbose: eprint(info['name'])
        #print(info['sequence'])
        left_resi = 999
        right_resi = -999
        for ii in range(left_cut, right_cut+1):
            resi = ii + (b_shift-a_shift)
            if resi >=0 and resi <length:
                left_resi = min( left_resi, resi )
                right_resi = max( right_resi, resi )
                raw_ddg += info['per_res_ddg'][resi]
                bonus_for_key_res = 0.0 if bonus_key == None else info[bonus_key].count(resi+1) * bonus_value
                weighted_ddg += ( bonus_for_key_res + info['per_res_ddg'][resi] ) * weights[ii+50]
                unweighted_ddg += bonus_for_key_res + info['per_res_ddg'][resi]
        for ii in hotspots_specific:
            resi = ii + (b_shift-a_shift)
            if resi >= left_resi and resi <= right_resi:
                hotspots_new_index.append( resi - left_resi + 1 )
        hotspots_new_index.sort()
        hotspots_str = ','.join( [ str(ii) for ii in hotspots_new_index ] )
        if verbose: eprint( "{:.2f}".format(weighted_ddg) )
            
        final_results.append( (tag, raw_ddg, unweighted_ddg, weighted_ddg, left_resi+1, right_resi+1, hotspots_str ) )
    
    return final_results



def get_top_motifs( motifs, num_motifs = 1000 ):
    cursor = num_motifs
    size = len(motifs)
    if size <= cursor: return motifs[:]
    raw = []
    unweighted = []
    weighted = []
    for r in motifs:
        raw.append(r[1])
        unweighted.append(r[2])
        weighted.append(r[3])
    raw.sort()
    unweighted.sort()
    weighted.sort()
    while True:
        cutoff_raw = raw[cursor]
        cutoff_unweighted = unweighted[cursor]
        cutoff_weighted = weighted[cursor]
        
        c = 0
        for m in motifs:
            if m[1] < cutoff_raw and m[2] < cutoff_unweighted and m[3] < cutoff_weighted:
                c += 1
        if c >= num_motifs:
            break
            
        cursor += 1
        
    cutoff_raw = raw[cursor]
    cutoff_unweighted = unweighted[cursor]
    cutoff_weighted = weighted[cursor]
    
    selected_motifs = []
    for m in motifs:
            if m[1] < cutoff_raw and m[2] < cutoff_unweighted and m[3] < cutoff_weighted:
                selected_motifs.append(m)
                
    return selected_motifs



standard_good_motifs = []
for k,cluster in clusters.items():
    any_P = False
    for name in cluster['pdb_names']:
        if ( "_P.pdb" in name):
            any_P = True
    if ( any_P ):
        continue
    final_results = select_best_motif( cluster, min_motif_len = 7, verbose=False, bonus_key=None)
    final_results.sort(key=lambda x:x[3])
    standard_good_motifs.append(final_results[0])



standard_good_motifs_top1000 = get_top_motifs( standard_good_motifs, num_motifs = args.num_motifs)



all_motifs = []
all_motifs.extend(standard_good_motifs_top1000)
uniq_motifs = {}


for m in all_motifs:
    uniq_motifs[ m[0] ] = m
    
    
for k,v in uniq_motifs.items():
    print('{} {:.2f} {:.2f} {:.2f} {} {} {}'.format(*v))




