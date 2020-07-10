#!/bin/bash

pdbtag=$1
num=$2
weights=$3


~/Rosetta/main/source/bin/fasol_perres.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-score:weights $weights \
	-set_weights cart_bonded 0.0 pro_close 1.25 \
	-s /home/dimaio/optE2/dualoptE/decoys/fasol/$pdbtag.pdb \
	-mute all -unmute fasol_refit \
	-overwrite -nooutput \
	-ignore_unrecognized_res \
	-crystal_refine \
	-flip_HNQ -no_optH false > opt_$num/$pdbtag.EXPSOL
