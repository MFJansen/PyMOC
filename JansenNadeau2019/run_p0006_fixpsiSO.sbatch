#!/bin/bash
#SBATCH --job-name=fixPsiSO
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=2:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --db=0.0006 --pickup='save_ref.npz' --diagfile='diags_p0006_fixPsiSO' --fixPsiSO

