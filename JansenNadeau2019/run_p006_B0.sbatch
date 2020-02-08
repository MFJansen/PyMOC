#!/bin/bash
#SBATCH --job-name=p006_B0
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=3:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --db=0.006 --B=0.0 --pickup='save_ref.npz' --diagfile='diags_p006_B0'

