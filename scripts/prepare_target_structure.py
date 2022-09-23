#! /usr/bin/env python
import os
import argparse
import subprocess
from pyrosetta import *
from pyrosetta.rosetta import *

####################
# This script performs PDB complex structure preparation for de novo design workflow
####################

# parse arguments for the job
parser = argparse.ArgumentParser(description='PDB structure preparation for binder design')
parser.add_argument('-i', nargs=1, help='Input structure file name', type=str, required=True)
parser.add_argument('-r', nargs=1, help='Receptor chain name in the input file. For example: A', type=str, required=True)
parser.add_argument('--just_sap', help='Only compute per residue SAP scores for the protein', action='store_true')
parser.add_argument('--just_relax', help='Relax the protein and compute SAP prior to manual trimming', action='store_true')
parser.add_argument('--trimmed', help='Finish preparation of manually trimmed target', action='store_true')
args = parser.parse_args()

## create dirs for the upcomming workflow
os.makedirs('../2_docking', exist_ok=True)
os.makedirs('../3_interface_design_on_docks', exist_ok=True)
os.makedirs('../4_motif_grafting', exist_ok=True)
os.makedirs('../5_interface_design_on_grafts', exist_ok=True)
os.makedirs('../6_final_filtering', exist_ok=True)

## if we are working with the structure requiring cleaning and relaxation

if not args.trimmed:

	## extract the desired receptor protein and strip miscellaneous molecules from pdb file

	print(f"Extracting the receptor chain from {args.i[0]}")
	struct_name = args.i[0][:-4]
	protein_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'MET',
					    'CYS', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 
					    'PHE', 'TRP', 'TYR', 'PRO', 'CYX', 'PYL', 'SEC', 'ASX', 'GLX']
	receptor_atoms = []
	with open(args.i[0], 'r') as f:
		for line in f:
			dat = line.split()
			if (dat[0] == 'ATOM') and (len(dat) > 10) and (dat[3] in protein_residues):
				if dat[4] == args.r[0]:
					receptor_atoms.append(line)

	with open("receptor.pdb", 'w') as o:
		for line in receptor_atoms:
			o.write(line)
		o.write("TER\n")

	print(f"Saved the target chain as receptor.pdb\n")

	## compure per residue SAP for the receptor protein

	print("Computing SAP scores per residue for receptor.pdb file\n")
	scoring_args = [os.environ['ROSETTA'],
				    "-parser:protocol", f"{os.environ['SCRIPTS']}/per_res_sap.xml",
				    "-beta_nov16",
				    "-renumber_pdb",
				    "-s", "receptor.pdb"]

	process = subprocess.Popen(scoring_args)
	returncode = process.wait()

	os.rename("receptor_0001.pdb",
			  "receptor_with_sap_scores.pdb")

	if not args.just_sap:
		## relax the protein in Rosetta

		print(f"Relaxing the protein in Rosetta\n")
		relax_args = [os.environ['ROSETTA'],
					  "-parser:protocol", f"{os.environ['SCRIPTS']}/relax_structure.xml",
					  "-beta_nov16",
					  "-s", "receptor_with_sap_scores.pdb"]

		process = subprocess.Popen(relax_args)
		returncode = process.wait()

		os.rename("receptor_with_sap_scores_0001.pdb",
				  "receptor_relaxed.pdb")


# if we are working with the structure requiring cleaning and relaxation

if not args.just_relax and not args.just_sap:

	## renumber structure continuously with Rosetta

	print("Renumbering the receptor protein\n")
	init()
	pose = pose_from_pdb("receptor_relaxed.pdb")
	pose.pdb_info(core.pose.PDBInfo(pose))
	pose.dump_pdb("receptor_relaxed_renumbered.pdb")

	## change chain names: ligand to A, receptor to B

	print(f"\nVerifying chain naming for the protein\n")
	table = {'A':'A', 'B':'B', 'C':'C', 'D':'D', 'E':'E', 'F':'F', 'G':'G', 'H':'H', 'I':'I', 'J':'J', 'K':'K', 'L':'L', 'M':'M', 'N':'N', 'O':'O', 'P':'P', 'Q':'Q', 'R':'R', 'S':'S', 'T':'T', 'U':'U', 'V':'V', 'W':'W', 'X':'X', 'Y':'Y', 'Z':'Z',
	        'a':'a', 'b':'b', 'c':'c', 'd':'d', 'e':'e', 'f':'f', 'g':'g', 'h':'h', 'i':'i', 'j':'j', 'k':'k', 'l':'l', 'm':'m', 'n':'n', 'o':'o', 'p':'p', 'q':'q', 'r':'r', 's':'s', 't':'t', 'u':'u', 'v':'v', 'w':'w', 'x':'x', 'y':'y', 'z':'z',
	        '1':'1', '2':'2', '3':'3', '4':'4', '5':'5', '6':'6', '7':'7', '8':'8', '9':'9'}
	table[args.r[0]] = 'A' # rename receptor chain to A
	results = []
	with open("receptor_relaxed_renumbered.pdb", 'r') as f:
		for line in f:
		    if line.startswith('ATOM'):
		        temp = list(line)
		        try:
		            temp[ 21 ] = table[ temp[21] ]
		            results.append( ''.join(temp) )
		        except KeyError:
		            print("Unknown chain identifier in the ATOM line!!")
		            exit(0)
		    elif line.startswith('SSBOND'):
		        temp = list(line)
		        try:
		            temp[ 15 ] = table[ temp[15] ]
		            temp[ 29 ] = table[ temp[29] ]
		            results.append( ''.join(temp) )
		        except KeyError:
		            print("Unknown chain identifier in the SSBOND line!!")
		            exit(0)
		    else:
		        results.append ( line )

	with open("receptor_relaxed_renumbered_chainchanged.pdb", 'w') as fout:
	    for line in results:
	        fout.write(line)

	## strip tables from the end of the file

	with open(f"receptor_relaxed_renumbered_chainchanged.pdb", 'r') as f:
		with open(f"receptor_relaxed_renumbered_chainchanged_striped.pdb", 'w') as o:
			for line in f:
				if "# All" in line:
					break
				else:
					o.write(line)






