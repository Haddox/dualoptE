#!/bin/sh

# core RR
bash scripts/merge.sh $1/features.core.db3 $1/*core.db3_*
sqlite3 $1/features.core.db3  < scripts/average_rotamer_recovery_by_amino_acid.sql > $1/rr_core_by_amino_acid.csv
sqlite3 $1/features.core.db3 < scripts/average_rotamer_recovery.sql > $1/rr_core_result

# int RR
bash scripts/merge.sh $1/features.int.db3 $1/*int.db3_*
sqlite3 $1/features.int.db3  < scripts/average_rotamer_recovery_by_amino_acid.sql > $1/rr_interface_by_amino_acid.csv
sqlite3 $1/features.int.db3 < scripts/average_rotamer_recovery.sql > $1/rr_interface_result

# rescore
ls $1/*.sc.static > temp.scores
./calc1dboltzmann.pl score temp.scores > $1/score_decoy_result
rm $1/*.sc.static

# mintest
ls $1/*.sc.min > temp.scoresmin
./calc1dboltzmann.pl score temp.scoresmin > $1/score_min_result
rm $1/*.sc.min

# docking
ls $1/*.docking.out > temp.docking
./calc1dboltzmann.pl ddg temp.docking > $1/score_docking_result
rm $1/*.docking.out

# xtal grads
cat $1/*.grads | awk '{ sum += $1; n++ } END { if (n > 0) print sum/n, n; }' > $1/xtal_grad_result
rm -rf _ros* $1/*.grads


# seqrecov
echo  /home/dimaio/optE2/dualoptE/decoys/stability_hugh/out_hbnet.csv > stability_results
ls $1/*stab >> stability_results

ln -s /home/dimaio/optE2/dualoptE/decoys/seqrecov/*COUNTS $1/
ls $1/*ENERGIES | sed 's/\.ENERGIES//' > ssm_results

# ./fitref <native_tags> <stability_tags> <ddg_tags> wt_1atatime_recov wt_stability wt_ddg wt_fixbb_recov wt_kldiv
./fitref ssm_results stability_results ddg_results 1.0 0.000001 0.0 0.0 0.5 2000 10000 > $1/seqrecov_result
rm $1/*COUNTS $1/*ENERGIES $1/*stab $1/*ddg ssm_results stability_results ddg_results

# expsol
cat $1/*.EXPSOL > $1/ALL.EXPSOL
./exp_desolvation.pl $1/ALL.EXPSOL --median 5 > $1/expsol_result
rm $1/*EXPSOL

# distrs
python ./distdstr_0.3.py $1/*_0001.pdb > $1/distr_result
rm $1/*_0001.pdb
