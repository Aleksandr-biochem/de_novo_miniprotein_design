# Mini-protein *de novo* design pipeline

This is an automated remake of the [*de novo* mini-protein design protocol by Longxing Cao and Brian Coventry](https://www.nature.com/articles/s41586-022-04654-9#data-availability). Clone this repository to start a new design workflow.  
  
The original scripts and the corresponding supplementary materials are available by [download link](files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/scripts_and_main_pdbs.tar.gz). Before working with this protocol, I would suggest reading through the `cao_2021_protocol_guide.txt` in aformentioned supplement to get a full idea about mini-binder *de novo* design. The full supplement is available [online](files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/download_supplement.txt)

## Before running any steps

Install the software and dependencies listed in **Appendix** at the end of this guide. Then inspect `export_paths_and_vars.sh` to set up your paths to the dependencies and the python environmet. Sourcing this file prior to running the scripts will make the workflow plesantly automated:

```
source export_paths_and_vars.sh
```

## 0. Prepare the scaffold library

### List scaffolds, convert to silent format 

Here is an example of HHH_bc 3-helical scaffold subset preparation from Longxing Cao and Brian Coventry scaffolds collection. You can use any other subset if you wish. Download and unpack `scaffolds` folder. There is no need to download it into the main design folder. These scaffolds can be used for any other job:

```
wget files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/scaffolds.tar.gz
tar -xvf scaffolds.tar.gz
cd scaffolds
```

Use `prepare_scaffolds.py` to list and prepare silent files from scaffolds pdbs. Scaffolds conversion to silent files may take several hours so it is better to launch it with slurm:

```
# the script has help section 
$SCRIPTS/prepare_scaffolds.py -h

# several folders for listing and preparation may be provided comma-separated, no spaces (e.x. HHH_bc,HHH_gr)
sbatch --job-name="scaffold_prep" --partition=cpu --nodelist=node_name $SCRIPTS/prepare_scaffolds.py -i HHH_bc
```

As a result, `scaffolds_HHH_bc.list` and `scaffold_HHH_bc_splits.list` are produced. Copy these files to the folder, where the binder design process will be performed.

### Convert scaffolds to poly-Valines for PatchDock

Use `prepare_polyV_scaffolds.py` to prepare poly-Valine scaffolds pdbs for docking. Provide a dir with scaffolds pdb to convert. It may be better to launch with slurm:

```
$SCRIPTS/prepare_polyV_scaffolds.py -h
sbatch --job-name="polyV_prep" --partition=cpu --nodelist=node_name $SCRIPTS/prepare_polyV_scaffolds.py -i HHH_bc 
```

As a result `HHH_bc_polyV.list` is produced. Copy this file into the folder, where the binder design process will be performed. **Note:** poly-Valine scaffolds should originate from exactly same scaffold files as in `.list` from previous step. 

## 1. Target preparation

### 1.1 Extract and relax the target protein

Place the pdb file containing the target of interest into `1_target_preparation` dir (there is already a sample file to try) and execute `prepare_target_structure.py`. This will perform:

- Extraction and continuous renumbering of the desired receptor protein chain
- Calculation of target per residue SAP
- Target relaxation in Rosetta
- Chain renaming (at this step the receptor chain is named A)

**Important note:** You may benefit from manual target extraction and trimming. This is optional but will speed up all subsquent calculations. Rosetta considers all residues of your target protein even if they are very far from the interface. Removing residues 20Å away from the interface should be safe. If you want to do this launch `prepare_target_structure.py` with `--just_relax` flag. Open the relaxed structure `receptor_relaxed.pdb` in a molecular viewer and delete some residues. Try aiming for the target to be around 200 residues and **don't change the name of the file**. Then run `prepare_target_structure.py` on trimmed file with `--trimmed` flag. 

```
cd 1_target_preparation

# the script has the help section
$SCRIPTS/prepare_target_structure.py -h

# below is an example of the preparation with manual trimming
$SCRIPTS/prepare_target_structure.py -i 7jzm.pdb -r B --just_relax
# trim receptor_relaxed.pdb, then
$SCRIPTS/prepare_target_structure.py -i receptor_relaxed.pdb -r B --trimmed

# or prepare without trimming in a single run
$SCRIPTS/prepare_target_structure.py -i 7jzm.pdb -r B
```

If you just want to assess per residue SAP scores for a target, use `--just_sap` flag. This will only perform the extraction of the desired protein chain and SAP calculation:

```
$SCRIPTS/prepare_target_structure.py -i 7jzm.pdb -r B --just_sap
```

### 1.2 Select target residues for PatchDock and RifDock

We are aiming to design binders making contacts with the hydrophobic residues on selected interface. The file `receptor_with_sap_scores.pdb` is a receptor structure with per residue SAP scores, which will help to make this choice. This structure can be conveniently visualized in PyMol using `per_res_sap.py` script provided in `scripts` dir.  

Start PyMol at the directory where the `receptor_with_sap_scores.pdb` and then:

```
# execute the script in PyMol terminal after loading the `receptor_with_sap_scores.pdb`
run path/to/per_res_sap.py
```

Select 3-9 residues for PatchDock to aim for. **Consider these when choosing:**
- At the end of the protocol, you will filter on how much contact your binders are making with these residues. Pick your residues such that a binder will be able to interact with all of them at once.
- Clumps of residues are better than sparse residues. 
- The bigger residue SAP score the better. However, any hydrophobic residues with a SAP score > 1 are probably good. But again, focus on clumps. If there are two distinct clumps, run this protocol twice. 
- Ideally, one would only pick hydrophobic residues (ALA, VAL, THR, ILE, PRO, LEU, MET, PHE, TYR, TRP). However, if your site is so hard that polar residues must be picked, this is acceptable (there's a flag to give to RifDock if you do this). Know though that by needing to pick polar residues, your success rate may be lowered.

Write selected comma-separated residues numbers into a text file `patchdock_residues.txt`. **Note:** In case you have trimmed the target protein, you have to carefully adjust the numbers of residues to match the trimmed structure. It is better to check the selection in a molecular viewer. 

```
# for example
echo '116,122,123,153,156,157,159,172' > patchdock_residues.txt
```

Next, write residues for RifGen as `rifgen_residues.txt`. A typical RifGen residue selection will include all of the PathDock residues and their immediate neighbors, usually about 10-25 residues. Fewer residues selected is often best because we want RifDock to focus on your PathDock residues.

```
# for example
echo '115,116,117,121,122,123,124,152,153,154,155,156,157,158,159,160,171,172,173' > rifgen_residues.txt
```

## 2. Docking

### 2.1 Run RifGen

Create a cache dir for RifGen somewhere on your system. This dir is populated on first run and then never changes. Execute `run_rifgen.py` providing this cache dir. RifGen run will take about 6-9 hours, so it is better to launch with slurm:

```
cd 2_docking

sbatch --job-name="rifgen" --partition=cpu --nodelist=node_name ../scripts/run_rifgen.py --cache /home/user/tmp/rifgen
```

Apart from the log file and `rifgen_out` dir, a file `receptor_centered_chainchanged.pdb` is generated. This target protein chain is centered, renamed to B and used for the further dockings, so don't move or rename it. **Note:** It is also critical to keep the output `rifgen.log` file as the information from it is used by further scripts.

### 2.2 Run PatchDock

PatchDock runs fast, however if you are dealing with thousands of scaffolds you may want to distribute tasks over several cpus. Use `prepare_patchdock.py` to set up PatchDock jobs providing the list of nodes and the number of cpus per node you want to use:

```
# the script has help section 
$SCRIPTS/prepare_patchdock.py -h

# distribute the jobs between 2 nodes, 10 cpus each
$SCRIPTS/prepare_patchdock.py -i ../HHH_bc_polyV.list -n node_name1,node_name2 -c 10
```

This script will create `patchdock_xforms` subdir with all the necessary files and scripts. Launch all jobs by executing:

```
cd patchdock_xforms
chmod +x run_patchdock_jobs.sh
./run_patchdock_jobs.sh
```

### 2.3 Run RifDock

Use `prepare_rifdock.py` to set up RifDock jobs as earlier at `2_docking` dir:

```
# the script has help section 
$SCRIPTS/prepare_rifdock.py -h

# distribute jobs between 2 nodes, 10 cpus each
$SCRIPTS/prepare_rifdock.py --cache /home/user/tmp/rifdock -i ../HHH_bc.list -n node_name1,node_name2 -c 10
```

RifDock requires manual `rifdock.flag` tuning. For that purpose, first launch the trial run, which will perform the docking for one scaffold only in about 4-5 mins:

```
./run_test_rifdock.sh
```

The `rifdock_out` dir will be generated. Inspect the output dockings:

```
# count output dockings using silent-tools
silentls rifdock_out/scaffold_name.silent | wc -l

# in addition you can inspect dockings visually. Exctract 10 random dockings:
cd rifdock_out
silentls scaffold_name.silent | shuf | head -n 10 | silentextractspecific
```

The tuning goal is to achieve a particular number of docking outputs per scaffold. By defaukt we aim for 300 (as specified in `-n_pdb_out_global` at `rifdock.flag`), but you may want to produce fewer or more docking outputs. If you have less outputs then expected after the test run, delete `rifdock_out` dir and change some `rifdock.flag` parameters manually. Here is a brief advice on what flags to tackle:
- `-require_hydrophobic_residue_contacts` by default is set to the number of PatchDock residues - 2. The less hydrophobic residue interactions you require, the more outputs you will get. Decrease by 1 and run test again.
- `-rif_dock:redundancy_filter_mag` controls the RMSD between outputs after clustering. The default is 0.5Å, you can lower it to get more outputs. A goal would be to ensure that your outputted interfaces are at least 0.3Å RMSD different from each other.
- `-xform_pos` passes the file from `scripts` dir defining the sampling rate. By default you have a `small_sampling_10deg_by1.7_1.5A_by_0.6.x` defining 10K pertrubations. You can provide `large_sampling_10deg_by1.1_2A_by0.35.x` to sample 100K perturbations. Going from the small file to the large file will increase your runtime by 10-fold, however, you will probably get more output.

See the full tuning discussion in original supplementary manual. You can look throgh `rifdock.flag` to investigate more options for tuning.

After tuning, make sure to clear `rifdock_out` and launch all dockings:

```
./run_rifdock_jobs.sh
```

Once the jobs are finished, collect all the outputs into a single silent file:

```
cat rifdock_out/*.silent > rifdock_out.silent
```

**Note:** the aim of 300 docking outputs per scaffold is provided by the authors and seemingly allows to achieve reasonably large sampling when having \~20K scaffolds. You may want to work with a smaller sampling but remember that at this version of the protocol large sampling appears crucial for success.

## 3. Interface design on dockings

### 3.1 Predictor and filtering

This step is relevant, when you have too many dockings for design. Estimate how much time your designs will take: a design of 20 complexes requires 1 cpu and about 4 hours. Also consider that you need to run full design for a sample of 1000 dockings to set up the filtering, so if you have 2-10K dockings it may be easier to design them all.

In case you want to apply filtering, execute `prepare_fast_design.py` providing the combined silent file with dockings and the computational resources for jobs set up. Then launch prepared jobs:

```
cd 3_interface_design_on_docks
$SCRIPTS/prepare_fast_design.py -i ../2_docking/rifdock_out.silent -n node_name1,node_name2 -c 10

# launch prepared jobs at fast_design dir
cd fast_design
./run_fast_design_jobs.sh
```

When jobs finish, collect the outputs using the python script providing the predictor and design dirs:

```
# at fast_design dir
$SCRIPTS/collect_outputs.py --predictor predictor --design pilot_design
```

This will produce `predictor_score_combined.sc`, `design_combined.sc` and `design_combined.silent`. Use the `.sc` score files and `predictor_filter.ipynb` code to obtain the list of filtered docks `predicted_fastdesign_tags.list`. Then slice these docks from combined silent file and use them for interface design:

```
cat predicted_fastdesign_tags.list | silentslice path/to/rifdock_out.silent > rifdock_out_filtered.silent
```

### 3.2 Interface design

Setup and run design using the `prepare_design.py` script, providing the desired docks set. If you have not conducted filtering, then it will be convinient to setup design directly at `3_interface_design_on_docks`. Otherwise you may create a neighboor dir `full_design` to setup design from there.

```
# optional
cd interface_design

# setup the jobs
$SCRIPTS/prepare_design.py -h
$SCRIPTS/prepare_design.py -i path/to/docks.silent -n node_name1,node_name2 -c 10

# launch
./run_design_jobs.sh

# collect results
$SCRIPTS/collect_outputs.py --design .
```

## 4. Motif grafting

### 4.1 Motif extraction

Extract and cluster motifs from designed complexes, using the `prepare_motifs.py`. Motif extraction from 1000 designs takes about 8 mins, so in case of large designs amount it's better to launch with slurm (check `slurm.out` to see info on the number of extracted motifs and clusters).

```
cd 4_motif_grafting

sbatch --job-name="motif_extraction" --partition=cpu --nodelist=node_name $SCRIPTS/prepare_motifs.py -i ../3_interface_design_on_docks/design_combined.silent
```

**Note:**
- This step can produce as much as 4 files for every structure you give it in very short order. If that number of files is going to overload your system, you can pass `--less_motifs` flag, which will change the ddg threshold for motif extraction to -25. This will probably reduce the number of outputs to 0.2 files per structure, but you'll still have the really good ones.
- You can pass `--dump_og` to the above command for better visualization of the motifs (this dumps the motif in the presense of the target molecule). The downside is that we will get 50% more files and probably 10x more data.

Pick the final motif list, providing the number of motifs you want to use. A 1000 motifs can be a good choice. Consider that it takes about 1h of calculations per scaffold batch–motif combination (1 cpu, 1-2GB ram).  

```
$SCRIPTS/finalize_motifs.py -n 1000
```

`motifs_with_hotspots.list` contains the motifs with their hotspots for grafting. Don't move or rename this file!

### 4.2 Grafting

Set up, run and summarize grafting with the following commands:

```
cd 4_motif_grafting

$SCRIPTS/prepare_grafting.py -i ../scaffold_splits.list -n node_name1,node_name2 -c 10

cd grafting
./run_grafting_jobs.sh

cd ../
$SCRIPTS/collect_outputs.py --grafting grafting
```

As a result you will get `graftings_combined.silent` containing the combined graftings.

## 5. Interface design on graftings

Inteface design on graftings is prepared, launched and summarized the same way as described in section 3.2. If you have too many graftings to design, filter them as described in section 3.1.

## 6. Final filtering

This step can be done automatically with `filter_designs.py`. 

```
cd 6_final_filtering
$SCRIPTS/designs_filter.py -h

# you can provide several silent files to filter. See `-h` for more options
scripts/designs_filter.py -i ../3_interface_design_on_docks/design_combined.silent,../5_interface_design_on_grafts/design_combined.silent
```

The script will print stats on filtering with hard cuts and best effort percentile, and produce silent with filtered designs and a plot `terms_distributios.png` with terms distributions among designs subjected to filtering.

If you pass `-p` flag with number/`all` velue, some of the filtered designs will be extracted as pdb into `designs_pdb` subdir. In order to extract binder sequence use `extract_protein_sequence.py` providing a `-i` directory of pdb files or `-f` a single pdb file. This script can be applied to extract sequence from pdbs of motifs, graftings and unfiltered designs too.

```
$SCRIPTS/extract_protein_sequence.py -i filtered_design/designs_pdb
```

**Note:** The designs filtration procedure may require manual filter adjustment, therefore you may want to download scores `.sc` files of the corresponding designs and the `final_filtering.ipynb` to experiment in hand. 

Here are coments on filtering step from the original *de novo* design pipeline description:

The authors have been using these cutoffs at the time of writing:
- contact_molecular_surface > 450
- ddg < -30
- binder_delta_sap > 12         # totally polar binders don't seem to work. This looks like the cutoff
- mismatch_probability < 0.1    # you can go up to 0.3 for beta-sheet designs
- sap_score < 35   # not strictly necessary for yeast surface display. But your designs will aggregate in E. coli above this
- score_per_res < -2.4 # Super topology dependent
- optional, uncomment it in the script if you need it: ss_sc > 0.77 # If you have tons of outputs, this will increase your success rate. However, it's not a robust filter

And then trying to optimize the following metrics until they reach the number they wish to order:
- ddg
- contact_patch
- target_delta_sap
- contact_molec_sq5_apap_target

Ideally, you'd be taking the 90th percentile and above according to these metrics. If you're only taking the top 20th percentile, you may need to go back and make more designs.

## Programs and environment setup

In order to run this protocol from start to finish, you will need the following software:

- **Rosetta**. You will need a relatively new version (newer than December 2020). Also, you need to compile with HDF5 support. (Add extras=hdf5 to your scons command)
- **DAlphaBall**. Part of rosetta. Go to rosetta/external/DAlpahBall and type "make". DAlphaBall.gcc will be produced.
- **PyRosetta**. All python scripts in protocol assume PyRosetta is installed within your python environment. 
- **PatchDock**. Available [online](https://bioinfo3d.cs.tau.ac.il/PatchDock/php.php) as a pre-compiled binary
- **RifDock**. Requires a separate build of Rosetta and must be compiled with **gcc7**. Available [at GitHub](https://github.com/rifdock/rifdock)
- **silent_tools**. Available [at GitHub](https://github.com/bcov77/silent_tools)
- **PsiPred**. After installation, you'll need to identify the location of runpsipred_single which will be referred to as RUNPSIPRED_SINGLE
- **motif_clustering/cluster**. - Get this from [GitHub](https://github.com/LongxingCao/ppi_tools)

You'll also need this database (somewhat optional, but gives better results)
ss_grouped_vall_all.h5 -- files.ipd.uw.edu/pub/modular_repeat_protein_2020/ss_grouped_vall_all.h5
