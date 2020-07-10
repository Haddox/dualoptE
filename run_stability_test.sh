#!/bin/sh

silenttag=$1
num=$2
weights=$3

rm -f opt_$num/$silenttag.sc.stab

~/Rosetta/main/source/bin/min_test.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-in:file:silent /home/dimaio/optE2/dualoptE/decoys/stability_hugh/split_$silenttag \
	-in:file:silent_struct_type binary \
	-score:weights $weights \
    -min::cartesian \
    -default_max_cycles 50 \
	-out:file:score_only opt_$num/$silenttag.sc.stab \
	-mute all
