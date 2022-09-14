#! /usr/bin/env python
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import argparse
from subprocess import Popen, PIPE

#############################
# This script preforms final design filtering
#############################

parser = argparse.ArgumentParser(description='Filter designs:')
parser.add_argument('-i', nargs=1, help='Silent files to combine and filter comma-separated', type=str, required=True)
parser.add_argument('-p', nargs=1, help="Output pdb structures of the filtered complexes. Number or 'all'", type=str, required=False)
parser.add_argument('-n', nargs=1, help="Optional. A desired number of designs to filter. By default it is set to 50%% of designs passing the hard cuts filters", type=int, required=False)
parser.add_argument('-hard_cut_only', help="In case you have only few designs passing hard cuts you may want to leave just these without best effort filtering", action='store_true', required=False)
args = parser.parse_args()

## All of the designs must pass these cutoffs

hard_cuts = {               # cutoff,  higher is better
    'ddg' :                     [ -30,  False],
    'contact_molecular_surface':[ 450,  True],
    'score_per_res' :           [-2.4, False],
    'mismatch_probability':     [ 0.1, False],
    'sap_score':                [  35, False],
    'binder_delta_sap':         [  12,  True],
#    'ss_sc':                    [0.77, True], # optional
}

## After the hard cuts find the best by the following metrics

best_effort = {             # higher is better
    'ddg':                         [  False],
    'contact_patch' :               [  True],
    'target_delta_sap':             [  True],
    'contact_molec_sq5_apap_target':[  True],
}

## supporting function for finding best effort designs after hard cuts

def top_x_by_multiple(data, desired_num):
    # creating lists where each design has it's rank within each term
    arg_sorted = np.argsort(data, axis=-1)
    ranked = np.zeros(arg_sorted.shape, np.int)
    for icol in range(len(ranked)):
        ranked[icol,arg_sorted[icol]] = np.arange(0, total, dtype=np.int)

    # first power of 2 where bigger_2 // 2 > total
    bigger_2 = 2**(np.floor(np.log(total)/np.log(2)).astype(int) + 2)
    assert( bigger_2 // 2 > total )
    percentiles = np.linspace(0, 1, bigger_2)
    space_size = bigger_2
    next_cut = bigger_2 // 2 - 1
    remaining = 0
    cutoff = total

    # binary search
    # We're trying to find a percentile where if we take everything better than
    #  that percentile from all categories, and look for members that exist in all
    #  lists, it's equal to the number of designs we want to keep
    while True:

        # this is the actual ranking process
        # take top X in each argsort and make sure they're in all top Xs
        eval_percentile = (percentiles[next_cut] + percentiles[next_cut+1]) / 2
        cutoff = eval_percentile * total
        mask = np.ones(total, np.bool)
        for icol in range(len(data)):
            mask &= ranked[icol] >= cutoff

        # end ranking
        remaining = mask.sum()
        # eprint(cutoff, next_cut, "/", bigger_2, remaining)

        if ( remaining == desired_num ):
            break
        space_size //= 2
        if ( space_size == 1 ):
            break
        if ( remaining < desired_num ):
            next_cut -= space_size // 2
        else:
            next_cut += space_size // 2
            
    print("Best effort == %ith percentile"%(eval_percentile*100))     
    return mask

designs = args.i[0].split(',')

## read all score files from specified dirs

score_dfs = []
for design_dir in dirs_to_parse:
	sc_file = glob.glob(f'{design_dir}/*.sc')[0]
	asdf = pd.read_csv(sc_file, sep='\s+')
	asdf['origin'] = design_dir
	print(f'In {design_dir} found {asdf.shape[0]} designs to filter')
	score_dfs.append(asdf)
score_df = pd.concat(score_dfs)
score_df.reset_index(drop=True, inplace=True)

## filter designs

# Print the pass rates for each term
print("")
print("=========================== Hard cuts: ===========================")
score_df['orderable'] = True
for pilot_term in hard_cuts:
    cut, higher_better = hard_cuts[pilot_term]
    ok_term = pilot_term.replace("_pilot", "") + "_ok"
    if ( higher_better ):
        score_df[ok_term] = score_df[pilot_term] >= cut
    else:
        score_df[ok_term] = score_df[pilot_term] <= cut
    score_df['orderable'] &= score_df[ok_term]
    print("%30s: %7.2f -- %5.0f%% pass-rate"%(pilot_term, cut, score_df[ok_term].sum() / len(score_df) * 100))


total = score_df['orderable'].sum()
print("")
print("                         %s (%.1f%%) designs remain"%(total, total/len(score_df)*100))

# check if any designs pass hard cuts
if total == 0:
	final_df = score_df
	print(f"Looks like no dessigns pass hard cuts. Try making mode designs or adjusting the filters after looking at the terms plots")
else:
	if args.n is None:
		number_to_order = total // 2
	else:
		number_to_order = args.n[0]

	if not args.hard_cut_only:
		after_hard_cuts = score_df[score_df['orderable']]
		data = np.zeros((len(best_effort), total))
		for iterm, term in enumerate(best_effort):
		    higher_better = best_effort[term][0]
		    if ( higher_better ):
		        data[iterm] = after_hard_cuts[term].values
		    else:
		        data[iterm] = -after_hard_cuts[term].values

		print("")
		is_in_the_top = top_x_by_multiple(data, number_to_order)
		final_df = after_hard_cuts[is_in_the_top]
		        
		print("")
		print("=========================== Best effort: ===========================")
		for term in best_effort:
		    higher_better = best_effort[term]
		    if ( higher_better ):
		        cut = final_df[term].min()
		        oks = (after_hard_cuts[term] >= cut).sum()
		    else:
		        cut = final_df[term].max()
		        oks = (after_hard_cuts[term] <= cut).sum()
		    
		    print("%30s: %7.2f -- %5.0f%% pass-rate"%(term, cut, oks / len(after_hard_cuts) * 100))

		print("")
		print("Final: %i designs"%(len(final_df)))
	else:
		final_df = score_df[score_df['orderable']]

	# save filtered designs
	final_df[['description']].to_csv(f"filtered_{'_'.join(dirs_to_parse)}/filtered_designs.list", index=None, header=None)
	filtered_designs = open(f"filtered_{'_'.join(dirs_to_parse)}/filtered_designs.silent", 'a')
	for design_dir in dirs_to_parse:
		p1 = Popen(["cat", f"filtered_{'_'.join(dirs_to_parse)}/filtered_designs.list"], stdout=PIPE)
		p2 = Popen(["silentslice", f"{design_dir}/designs_combined.silent"], stdin=p1.stdout, stdout=filtered_designs)
		returncode = p2.wait()
	filtered_designs.close()


## convert designs to pdb if needed

if not args.p is None:
	if args.p[0] == 'all':
		num_to_extract = final_df.shape[0]
	else:
		num_to_extract = int(args.p[0])	

	print(f"Converting {args.p[0]} of filtered designs into pdb")

	os.makedirs(f"filtered_{'_'.join(dirs_to_parse)}/designs_pdb", exist_ok=True)
	os.chdir(f"filtered_{'_'.join(dirs_to_parse)}/designs_pdb")
	p1 = Popen(["silentls", f'../filtered_designs.silent'], stdout=PIPE)
	p2 = Popen(["shuf"], stdin=p1.stdout, stdout=PIPE)
	p3 = Popen(["head", "-n", str(num_to_extract)], stdin=p2.stdout, stdout=PIPE)
	p4 = Popen(["silentextractspecific", f'../filtered_designs.silent'], stdin=p3.stdout, stdout=PIPE)
	returncode = p4.wait()
	os.chdir("../../")

## plot terms destributions

print("Plotting terms distributions for the provided designs")
fig, axes = plt.subplots(3, 2, figsize=(20, 30))
axes = axes.flatten()
for i, cut in enumerate(hard_cuts):
    # add cuts
    cutoff, higher_better = hard_cuts[cut]
    if higher_better:
        lim = score_df[cut].max()
        axes[i].axvspan(cutoff, lim, facecolor='#D6E1C9', alpha=0.5, label='cutoff')
    else:
        lim = score_df[cut].min()
        axes[i].axvspan(lim, cutoff, facecolor='#D6E1C9', alpha=0.5, label='cutoff')
    
    sns.histplot(score_df[cut], color='#3461d2', element="step", alpha=0.3,  ax=axes[i])

    axes[i].set_title(cut, fontsize=20, pad=15)
    axes[i].set_xlabel(xlabel=cut, fontsize=15, labelpad=15)
    axes[i].set_ylabel(ylabel='Count', fontsize=15, labelpad=15)
    
plt.savefig(f"filtered_{'_'.join(dirs_to_parse)}/terms_distributions.png", dpi=300)




