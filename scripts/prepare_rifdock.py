#! /usr/bin/env python
import os
from os.path import exists
import argparse
from subprocess import Popen

####################
# This script prepares prepares RifDock runs
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='Run RifDock')
parser.add_argument('--cache', nargs=1, help='Cache dir for RifDock.', type=str, required=True)
parser.add_argument('-i', nargs=1, help='The list of input scaffolds', type=str, required=True)
parser.add_argument('-n', nargs=1, help='Node list to send jobs to, comma-separated', type=str, required=True)
parser.add_argument('-c', nargs=1, help='Number of CPUs per node to work with', type=int, required=True)
parser.add_argument('--use_polar', help='Pass this flag if there are polar residues among PatchDock residues list.', action='store_true')
args = parser.parse_args()

## write rifdock.flag file

if exists("../1_target_preparation/patchdock_residues.txt"):
	with open("../1_target_preparation/patchdock_residues.txt", 'r') as inp:
		patchdock_residues = inp.readlines()[0].strip()
else:
	raise Exception("Can not find the patchdock_residues.txt file!")

with open('rifdock.flag', 'w') as f:

	# first extract necessary information from rifgen.log

	lines_from_rifgen_log = []
	with open('rifgen.log', 'r') as log:
		for line in log:
			if 'what you need for docking' in line or len(lines_from_rifgen_log) > 0:
				lines_from_rifgen_log.append(line)

	for line in lines_from_rifgen_log:
		f.write(line)

	# now other flags		

	f.write("################################ Constant paths ##########################################\n")
	f.write(f"-database {os.environ['ROSETTA_RIF_DATABASE']}\n\n")

	f.write(f"# A cache directory. Populated on first-run and then never changes.\n")
	f.write(f"-rif_dock:rotrf_cache_dir {args.cache[0]}\n")

	f.writelines("""############################### Last minute options ######################################

-outputsilent  # Make sure you enable this for your production run. Dealing with 6M pdbs is less than ideal

-n_pdb_out_global 300                     # n_pdb_out controls how many per patchdock output. This is how many total
-rif_dock:redundancy_filter_mag 0.5       # The "RMSD" cluster threshold for the output. Smaller numbers give more, but redundant output

#################################### Flags that control output ##########################################

-rif_dock:outdir  rifdock_out             # the output folder for this run
-rif_dock:dokfile all.dok                 # the "score file" for this run

-rif_dock:n_pdb_out 30                    # max number of output pdbs

#-rif_dock:target_tag conf01              # optional tag to add to all outputs
-rif_dock:align_output_to_scaffold false  # If this is false, the output is aligned to the target

# Pick either one of the following or none
# (None)                                  # Output target + scaffold. But scaffold may be poly ALA with rifres based on scaffold_to_ala
#-output_scaffold_only                    # Output just the scaffold. But scaffold may be poly ALA with rifres based on scaffold_to_ala
-output_full_scaffold                     # Output target + scaffold. Scaffold retains input sequence plus rifres
#-output_full_scaffold_only               # Output just the scaffold. Scaffold retains input sequence plus rifres


############################ Flags that affect runtime/search space ####################################

-beam_size_M 5                            # The number of search points to using during HSearch
-hsearch_scale_factor 1.2                 # The default search resolution gets multiplied by this. People don't usually change this.

#-rif_dock:tether_to_input_position 3     # Only allow results within this "RMSD" of the input scaffold

-rif_dock:global_score_cut -8.0          # After HSearch and after HackPack, anything worse than this gets thrown out


##################### Flags that only affect the PatchDock/RifDock runs ################################
# Uncomment everything here except seeding_pos if running PatchDock/RifDock

#-rif_dock:seeding_pos ""                  # Either a single file or a list of seeding position files
-rif_dock:seeding_by_patchdock true       # If true, seeding_pos is literally the PatchDock .out file
                                           # If false, seeding_pos file is list of transforms. 
                                           #   (Each row is 12 numbers. First 9 are rotation matrix and last 3 are translation.)
-rif_dock:patchdock_min_sasa  1000        # Only take patchdock outputs with more than this sasa
-rif_dock:patchdock_top_ranks 2000        # Only take the first N patchdock outputs\n""")

	f.write(f"-rif_dock:xform_pos {os.environ['SCRIPTS']}/small_sampling_10deg_by1.7_1.5A_by_0.6.x\n")
	f.writelines("""                                       # Which xform file do you want to use. Difference is how many degrees do you want to 
                                           #   deviate from the PatchDock outputs. Pick one from here:
                                           #                 /home1/05571/bcov/sc/scaffold_comparison/data/xform_pos_ang*


-rif_dock:cluster_score_cut -6.0          # After HackPack, what results should be thrown out before applying -keep_top_clusters_frac
-rif_dock:keep_top_clusters_frac 1.0      # After applying the cluster_score_cut, what fraction of remaining seeding positions should survive?
                                         
-rosetta_score_each_seeding_at_least 1    # When cutting down the results by rosetta_score_fraction, make sure at least this many from each 
                                           #   seeding position survive

-only_load_highest_resl                   # This will make rifdock use less ram. Highly recommended for the patchdock protocol.


##################### Advanced seeding position flags ##################################################

#-rif_dock:seed_with_these_pdbs *.pdb      # List of scaffolds floating in space above the target that you would like to use instead.
                                           #   of numeric seeding positions. The target shouldn't be present and the scaffold must match exactly
                                           #   Use this instead of -seeding_pos
#-rif_dock:seed_include_input true         # Include the input pdb as one of the pdbs for -seed_with_these_pdbs

#-rif_dock:write_seed_to_output true       # Use this if you want to know which output came from which seeding position


##################### Flags that affect how things are scored ##########################################

#-use_rosetta_grid_energies true/false     # Your choice. If True, uses Frank's grid energies during Hackpack

-hbond_weight 2.0                          # max score per hbond (Rosetta's max is 2.1)
-upweight_multi_hbond 0.0                  # extra score factor for bidentate hbonds (BrianC recommends don't do this)
-min_hb_quality_for_satisfaction -0.25     # If using require_satisfaction (or buried unsats). How good does a hydrogen bond need to be to "count"?
                                           #   The scale is from -1.0 to 0 where -1.0 is a perfect hydrogen bond.
-scaff_bb_hbond_weight 2.0                 # max score per hbond on the scaffold backbone 

-favorable_1body_multiplier 0.2            # Anything with a one-body energy less than favorable_1body_cutoff gets multiplied by this
-favorable_1body_multiplier_cutoff 4       # Anything with a one-body energy less than this gets multiplied by favorable_1body_multiplier
-favorable_2body_multiplier 5              # Anything with a two-body energy less than 0 gets multiplied by this

-user_rotamer_bonus_constant 0             # Anything that makes a hydrogen-bond, is a hotspot, or is a "requirement" gets this bonus
-user_rotamer_bonus_per_chi 0              # Anything that makes a hydrogen-bond, is a hotspot, or is a "requirement" gets this bonus * number of chis

-rif_dock:upweight_iface 2.0               # During RifDock and HackPack. rifres-target interactions are multiplied by this number


################ stuff related to picking designable and fixed positions #################

#### if you DO NOT supply scaffold_res files, this will attempt to pick which residues on the scaffold
#### can be mutated based on sasa, internal energy, and bb-sc hbonds
-scaffold_res_use_best_guess true

#### if scaffold_res is NOT used, this option will cause loop residues to be ignored
#### scaffold_res overrides this
-rif_dock::dont_use_scaffold_loops false

#### these cause the non-designable scaffold residues to still contribute sterically
#### and to the 1 body rotamer energies. use these flags if you have a fully-designed scaffold
#-rif_dock:scaffold_to_ala false
#-rif_dock:scaffold_to_ala_selonly true
#-rif_dock:replace_all_with_ala_1bre false
#### if you don't have a fully designed scaffold, treat non-designable positions as alanine
-rif_dock:scaffold_to_ala true            # Brian thinks that converting the whole scaffold to alanine works better during rosetta min
-rif_dock:scaffold_to_ala_selonly false
-rif_dock:replace_all_with_ala_1bre true



#################################### HackPack options #####################################
-hack_pack true                            # Do you want to do HackPack? (Probably a good idea)
-rif_dock:hack_pack_frac  1.0              # What fraction of your HSearch results (that passed global_score_cut) do you want to HackPack?


############################# rosetta re-scoring / min stuff ###################################

-rif_dock:rosetta_score_cut -10.0                    # After RosettaScore, anything with a score worse than this gets thrown out

-rif_dock:rosetta_score_fraction 0                   # These two flags greaty affect runtime!!!!!
-rif_dock:rosetta_min_fraction 1                     # Choose wisely, higher fractions give more, better output at the cost of runtime

-rif_dock:rosetta_min_at_least 30                    # Make sure at least this many survive the rosetta_min_fraction
-rif_dock:rosetta_min_at_most 300                    # Make sure no more than this get minned
-rif_dock:rosetta_score_at_most  3000              # Make sure that no more than this many go to rosetta score

-rif_dock:replace_orig_scaffold_res false            # If you converted to poly ALA with scaffold_to_ala, this puts the original residues
                                                     #   back before you do rosetta min.
-rif_dock:override_rosetta_pose true                 # Brian highly recommends this flag. This prevents the minimized pose from being output
-rif_dock:rosetta_min_scaffoldbb false               # Set BB movemap of scaffold to True
-rif_dock:rosetta_min_targetbb   false               # Set BB movemap of target to true
-rif_dock:rosetta_hard_min false                     # Minimize with the "hard" score function (alternative is "soft" score function)

-rif_dock:rosetta_score_rifres_rifres_weight   0.6   # When evaluating the final score, multiply rifres-rifres interactions by this
-rif_dock:rosetta_score_rifres_scaffold_weight 0.4   # When evaluating the final score, multiply rifres-scaffold interactions by this
                                                     #  These two flags only get used if the interaction is good. Bad interactions are
                                                     #    full weight.


######################### Special flags that do special things #################################

#-hack_pack_during_hsearch False           # Run HackPack during the HSearch. Doesn't usually help, but who knows.
#-require_satisfaction 4                   # Require at least this many hbonds, hotspots, or "requirements"
#-require_n_rifres  3                      # Require at least this may rifres (not perfect)

#-requirements 0,1,2,8                     # Require that certain satisfactions be required in all outputs
                                           # If one runs a standard RifDock, these will be individual hydrogen bonds to specific atoms
                                           # If one uses hotspots during rifgen, these will correspond the the hotspots groups
                                           #   However, due to some conflicts, these will also overlap with hydrogen bonds to specific atoms
                                           # Finally, if one uses a -tuning_file, these will correspond to the "requirements" set there


######################### Hydrophobic Filters ##################################################
# These are rather experimental flags. You'll have to play with the values.
# Hydrophobic ddG is roughly fa_atr + fa_rep + fa_sol for hydrophobic residues.

#-hydrophobic_ddg_cut -12                  # All outputs must have hydrophobic ddG at least this value\n""")

	f.write(f"-require_hydrophobic_residue_contacts {len(patchdock_residues.split(',')) - 2}  # All outputs must make contact with at least this many target hydrophobics\n")

	f.writelines("""#-one_hydrophobic_better_than -2           # Require that at least one rifres have a hydrophobic ddG better than this
#-two_hydrophobics_better_than -2          # Require that at least two rifres have a hydrophobic ddG better than this
#-three_hydrophobics_better_than -1        # Require that at least three rifres have a hydrophobic ddG better than this

# This next flag affects the *_hydrophobics_better_than flags. A rifres can only be counted towards those flags if it passes this one.
#-hydrophobic_ddg_per_atom_cut -0.3        # Require that hydrophobics for the *_hydrophobics_better_than flags have at least this much 
                                           #  ddG per side-chain heavy atoms.\n""")

	f.write(f"-hydrophobic_target_res {patchdock_residues}    # If you want your selection of hydrophobic residues to include only a subset of the ones\n")
	f.write("                                                 #  you selected for the target_res, place that selection here with commas.\n")

	if args.use_polar:
		f.write("-count_all_contacts_as_hydrophobic        # Use this if some of your hydrophobic_target_res aren't actually hydrophobic amino acids\n")

	f.writelines("""-hydrophobic_ddg_weight 3
######################### options to favor existing scaffold residues ##########################
-add_native_scaffold_rots_when_packing 0 # 1
-bonus_to_native_scaffold_res          0 # -0.5


################################# Twobody table caching ####################################

# RifDock caches the twobody tables so that you can save time later. If you use the same scaffolds
#  in the same directory mulitple times. This is a good idea. Otherwise, these take up quite
#  a bit of space and it might be smart to turn the caching off.

-rif_dock:cache_scaffold_data False
-rif_dock:data_cache_dir  /home/me/rifdock/scaff_cache




############################################################################################
############################################################################################
############################ END OF USER ADJUSTABLE SETTINGS ###############################
############################################################################################
############################################################################################


#### to use -beta, ask will if you don't want to use -beta
-beta
-score:weights beta_soft
-add_orbitals false

#### HackPack options you probably shouldn't change
-rif_dock:pack_n_iters    2
-rif_dock:pack_iter_mult  2.0
-rif_dock:packing_use_rif_rotamers        true
-rif_dock:extra_rotamers                  false
-rif_dock:always_available_rotamers_level 0


#### details for how twobody rotamer energies are computed and stored, don't change
-rif_dock:rotrf_resl   0.25
-rif_dock:rotrf_spread 0.0
-rif_dock:rotrf_scale_atr 1.0

### Brian doesn't know what these flags do
-rif_dock::rf_resl 0.5
-rif_dock::rf_oversample 2
-rif_dock:use_scaffold_bounding_grids 0
-rif_dock:target_rf_oversample 2


 # disulfides seem to cause problems... ignoring them isn't really an issue unless
 # you do bbmin where there should be disulfides
-detect_disulf 0


-mute core.scoring.ScoreFunctionFactory
-mute core.io.pose_from_sfr.PoseFromSFRBuilder)""")


## write commands into a separate file and list

commands = []
cwd = os.getcwd()
with open(f'rifdock_commands.list', 'w') as o:
	with open(args.i[0], 'r') as f:
		for line in f:
			scaffold = line.strip()
			rifdock_args = [os.environ['RIFDOCK'],
						    f"@{cwd}/rifdock.flag",
						    "-scaffolds", scaffold,
						    "-seeding_pos",
						    f"{cwd}/patchdock_xforms/patchdock_out/{scaffold.split('/')[-1][:-4]}_0001.out",
						    "> rifdock.log"]
			commands.append(' '.join(rifdock_args))
			o.write(f"{' '.join(rifdock_args)}\n")

## write sbatch job scripts

os.makedirs("rifdock_jobs", exist_ok=True)
nodes = args.n[0].split(',')
for i, node in enumerate(nodes):
	for j in range(args.c[0]):
		with open(f"rifdock_jobs/rifdock_job{i}_{j}.sh", 'w') as f:
		    f.write('#!/bin/bash\n')
		    f.write('\n')
		    f.write(f'#SBATCH --job-name="rifdock{i}_{j}"\n')
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
				with open(f"rifdock_jobs/rifdock_job{i}_{j}.sh", 'a') as f:
					f.write(f"{commands[k]}\n")
					f.write("\n")
					k += 1

# prepare a script to launch all jobs
with open(f'run_rifdock_jobs.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	for i in range(len(nodes)):
		for j in range(args.c[0]):
			f.write(f'sbatch rifdock_jobs/rifdock_job{i}_{j}.sh\n')

# prepare a script to run test job for tuning
with open(f'run_test_rifdock.sh', 'w') as f:
	f.write('#!/bin/bash\n')
	f.write('\n')
	f.write(f"{commands[0]}\n")

