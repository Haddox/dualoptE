#!/bin/sh

pdbtag=$1
num=$2
weights=$3

phenix.rosetta.run_phenix_interface \
	~/Rosetta/main/source/bin/rosetta_scripts.python.linuxgccrelease \
    -database ~/Rosetta/main/database \
	-parser:protocol refine.xml \
	-parser:script_vars symmdef=/home/dimaio/optE2/dualoptE/decoys/xtal_refine_beta16/$pdbtag.symm wts=$weights outfile=opt_$num/$pdbtag.grads  \
	@flags \
	@flags_$num \
	-s /home/dimaio/optE2/dualoptE/decoys/xtal_refine_beta16/$pdbtag'_clean_0001.pdb' \
	-cryst::mtzfile /home/dimaio/optE2/dualoptE/decoys/xtal_refine_beta16/$pdbtag'-sf.mtz' \
	-crystal_refine \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-renumber_pdb false \
	-out::nooutput \
	#-mute all -unmute protocols.cryst.cryst_movers &> opt_$num/$pdbtag.grad_detail

