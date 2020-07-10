#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.docking.out

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	-parser:protocol ddg.xml \
	-parser:script_vars wts=$weights \
	@flags \
	@flags_$num \
	-in:file:silent /home/dimaio/optE2/dualoptE/decoys/docking/$pdbtag.trim2 \
	-in:file:silent_struct_type binary \
	-force_silent_bitflip_on_read \
	-score:weights $weights \
	-set_weights cart_bonded 0.0 pro_close 1.25 \
	-in:file:native  /home/dimaio/optE2/dualoptE/decoys/docking/$pdbtag'_bound_native.pdb' \
	-evaluation:rmsd NATIVE 1 FULL \
	-out:file:score_only opt_$num/$pdbtag.docking.out  \
	-silent_read_through_errors \
	-mute all
