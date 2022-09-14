#!/bin/bash

# activate your python environment
source /home/akovalenko/protein_design/venv/bin/activate

# export your path to scripts folder
export SCRIPTS=/home/akovalenko/protein_design/template_binder_de_novo_design/scripts

# export path to silent tools
export PATH=$PATH:/home/akovalenko/tools/silent_tools

# export path to mdutils
export PYTHONPATH=$PYTHONPATH:/home/akovalenko/cov2/github/md-utils

# export your path to Rosetta
# export ROSETTA=/home/olebedenko/tools/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/rosetta_scripts.mpi.linuxgccrelease
export ROSETTA=/home/olebedenko/tools/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/rosetta_scripts.hdf5.linuxgccrelease

# export your path to Rosetta database linked to RifDock
export ROSETTA_RIF_DATABASE=/home/olebedenko/tools/rosetta_src_2018.09.60072_bundle/main/database

# export your path to RifGen
export RIFGEN=/home/olebedenko/tools/rifdock/build_gcc7/apps/rosetta/rifgen

# export your path to RifDock
export RIFDOCK=/home/olebedenko/tools/rifdock/build_gcc7/apps/rosetta/rif_dock_test

# export your path to PatchDock
export PATCHDOCK=/home/akovalenko/tools/PatchDock/patch_dock.Linux

# export path to PsiPred
export PSIPRED=/home/akovalenko/tools/psipred/runpsipred_single

# export path to DAlpahBall
export DALPAHBALL=/home/olebedenko/tools/rosetta_bin_linux_2021.16.61629_bundle/main/source/external/DAlpahBall/DAlphaBall.gcc

# export path to ss_grouped_vall_all.h5 database
export SS_GROUPED_VALL_ALL=/home/akovalenko/tools/ss_grouped_vall_all.h5

# export path to motif clusterer
export CLUSTER=/home/akovalenko/tools/ppi_tools/motif_clustering/cluster