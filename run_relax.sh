#!/bin/bash

pdbtag=$1
num=$2
weights=$3

~/Rosetta/main/source/bin/relax.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s /home/dimaio/optE2/dualoptE/decoys/xtal_refine_beta16/$pdbtag'_clean_0001.pdb' \
	-score:weights $weights \
	-crystal_refine \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-default_max_cycles 200 \
	-relax:script cartminpack.script \
	-relax:min_type lbfgs_armijo_nonmonotone \
	-out:prefix opt_$num/ \
	#-mute all
