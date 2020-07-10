#!/bin/bash

pdbtag=$1
partnum=$2
num=$3
weights=$4

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s /home/dimaio/optE2/dualoptE/decoys/rotrecov/interface/$pdbtag'_clean.pdb' \
	-score:weights $weights \
	-set_weights cart_bonded 0.0 pro_close 1.25 \
	-mapfile /home/dimaio/optE2/dualoptE/decoys/rotrecov/interface/$pdbtag.map \
	-out::nooutput \
	-parser:protocol features.xml \
	-parser:script_vars outdir=opt_$num pdb=$pdbtag resfile=/home/dimaio/optE2/dualoptE/decoys/rotrecov/interface/$pdbtag"_cleanwat_0001_eval.resfile" database_partition=$partnum type=int \
	-mute all #-ignore_unrecognized_res

