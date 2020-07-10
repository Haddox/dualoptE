#!/bin/bash
#SBATCH -p long
#SBATCH -n 1 
#SBATCH -c 1
#SBATCH --mem=8g
#SBATCH --job-name=HH_run18D_master
#SBATCH -o opt.out
#SBATCH -e opt.err

source activate optE
python -u optimize.dask.py
