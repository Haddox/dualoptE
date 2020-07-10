#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.sc.static

~/Rosetta/main/source/bin/score_jd2.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-ignore_unrecognized_res \
	-in:file:silent /home/dimaio/optE2/dualoptE/decoys/monomer/$pdbtag.relax.out \
	-in:file:silent_struct_type binary \
	-force_silent_bitflip_on_read \
	-out:file:scorefile opt_$num/$pdbtag.sc.static \
	-score:weights $weights \
	-in:file:native /home/dimaio/optE2/dualoptE/decoys/monomer/$pdbtag'_clean.pdb' \
	-evaluation:rmsd NATIVE 1 FULL \
	-keep_input_scores false \
	-mute all
