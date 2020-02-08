#!/bin/bash
#SBATCH --job-name=p006
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=3:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --db=0.006 --pickup='save_ref.npz' --diagfile='diags_p006'

