#!/bin/bash

pdbtag=$1
num=$2
weights=$3

~/Rosetta/main/source/bin/calc_ssm_energies.default.linuxgccrelease \
		@flags \
		@flags_$num \
		-s /home/dimaio/optE2/dualoptE/decoys/seqrecov/$pdbtag'_0001.pdb' \
		-score:weights $weights \
    	-set_weights cart_bonded 0.0 pro_close 1.25 \
		-out::nooutput \
	    -mute all -unmute calc_ssm_energies > opt_$num/$pdbtag.ENERGIES

