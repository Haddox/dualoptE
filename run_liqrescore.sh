#!/bin/sh

silent=$1
tag=`echo $1 | awk -F. '{print $1}'`
num=$2
weights=$4

rm -f opt_$num/$tag.sc.liqsim

~/Rosetta/main/source/bin/score_jd2.default.linuxgccrelease \
	-in:file:silent  /usr/lusers/dimaio/dimaio/decoys/liquid_snapshots/$silent \
	-in:file:fullatom \
	-symmetry_definition dummy \
	@flags \
	@flags_$num \
	@flags_liq_$num \
	-score:weights $weights \
	-score:grpelec_fade_type grpsubtract \
	-score:elec_r_option true \
	-score:elec_max_dis 6.0 \
	-score:fa_max_dis 6.0 \
	-out:file:scorefile opt_$num/$silent.sc.liqsim \
	-overwrite \
	-score:elec_group_extrafile /usr/lusers/dimaio/dimaio/decoys/liquid_snapshots/params/$tag.groupdef \
	-extra_res_fa /usr/lusers/dimaio/dimaio/decoys/liquid_snapshots/params/$tag.params

