#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.sc.min

~/Rosetta/main/source/bin/min_test.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-ignore_unrecognized_res \
	-in:file:silent /home/dimaio/optE2/dualoptE/decoys/monomer/$pdbtag.relax.out.trim \
	-in:file:silent_struct_type binary \
	-force_silent_bitflip_on_read \
	-min:scoreonly \
	-out:file:scorefile opt_$num/$pdbtag.sc.min \
	-score:weights $weights \
	-in:file:native /home/dimaio/optE2/dualoptE/decoys/monomer/$pdbtag'_clean.pdb' \
	-evaluation:rmsd NATIVE 1 FULL \
	-min:minimizer lbfgs_armijo_nonmonotone \
	-min::pack \
	-default_max_cycles 100 \
	-min:cartesian \
	-keep_input_scores false \
	-mute all
